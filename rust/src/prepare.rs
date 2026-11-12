//! Everything the test needs before it runs: the alignment as residue codes, the locus plateaus and rates, and the
//! pair structures (each taxon's candidate list with m, its floor and the ranking key; the closest match at every
//! locus; the classes of identical genomes; the full m up to FULL_M_MAX_TAXA).
//!
//! locusmasker.py walks pairs of taxon blocks with matrix products and reduces block results; its reductions are
//! order-independent by construction, so this port walks the rows of the pair matrix directly, one taxon against
//! every other, with the arithmetic of its compiled pair kernel (_pair_medians_py, _closest_rows_py).
use crate::input::Loci;
use crate::npmath::*;
use crate::pyfmt::commas;
use crate::{log, Args};
use rayon::prelude::*;
use std::collections::HashMap;

pub const MIN_LOCI_FOR_M: usize = 5;
pub const MIN_LOCI_PER_TAXON: i64 = 10;
pub const MIN_FIT_PAIRS: usize = 8;
pub const BS_PER_COMP: usize = 10;
pub const BS_POOL: usize = 100;
pub const BS_MIN: i64 = 5;
pub const BS_SEED: u64 = 0x5eed1517;
pub const SWEEP_BLOCK_MIN: usize = 1024;
pub const CANDIDATES: usize = 200;
pub const FULL_M_MAX_TAXA: usize = 20000;
pub const BAND_SLOTS: usize = 4;
const SWEEP_CHUNK_BYTES: usize = 512 << 10;       // the planes of a chunk of rows: about a core's share of the cache

/// the residues by locus (each taxon's row holds the loci one after another) and, up to SWEEP_BLOCK_MIN taxa, the
/// counts of every pair at every locus
pub struct Core {
    pub n: usize, pub nloc: usize, pub lt: usize,
    pub xl: Vec<i8>, pub off: Vec<usize>, pub keep: Vec<bool>, pub a_: usize,
    pub slab: Option<Slab>, pub bits: Bits,
}

/// the residues as bit planes, for counting: per taxon a presence plane and the planes of the residue code's bits, each
/// locus starting at a word. Two sequences share a site where both presence bits are set, and agree there where every
/// code plane agrees: two popcounts per 64 sites, the same integers as counting site by site
pub struct Bits { pub planes: usize, pub wpr: usize, pub woff: Vec<usize>, pub data: Vec<u64>, popcnt: bool, avx2: bool }
impl Bits {
    fn build(n: usize, xl: &[i8], lt: usize, off: &[usize], a_: usize) -> Bits {
        let planes = if a_ <= 4 { 2 } else { 5 }; let nloc = off.len() - 1;
        let mut woff = Vec::with_capacity(nloc + 1); let mut acc = 0; for li in 0..nloc { woff.push(acc); acc += (off[li + 1] - off[li]).div_ceil(64); } woff.push(acc);
        let wpr = acc; let stride = (planes + 1) * wpr; let mut data = vec![0u64; n * stride];
        if stride > 0 {
            data.par_chunks_mut(stride).enumerate().for_each(|(g, row)| {
                let x = &xl[g * lt..(g + 1) * lt];
                for li in 0..nloc {
                    for (k, &v) in x[off[li]..off[li + 1]].iter().enumerate() {
                        if v < 0 { continue; }
                        let (w, b) = (woff[li] + k / 64, k % 64);
                        row[w] |= 1u64 << b;
                        for p in 0..planes { if (v as u8 >> p) & 1 == 1 { row[(p + 1) * wpr + w] |= 1u64 << b; } }
                    }
                }
            });
        }
        #[cfg(target_arch = "x86_64")]
        let (popcnt, avx2) = (std::is_x86_feature_detected!("popcnt"), std::is_x86_feature_detected!("popcnt") && std::is_x86_feature_detected!("avx2") && std::env::var_os("LOCUSMASKER_NO_AVX2").is_none());
        #[cfg(not(target_arch = "x86_64"))]
        let (popcnt, avx2) = (false, false);
        Bits { planes, wpr, woff, data, popcnt, avx2 }
    }
    #[inline(always)]
    fn count_impl(&self, g: usize, c: usize, li: usize) -> (u32, u32) {
        let (w0, w1) = (self.woff[li], self.woff[li + 1]); let stride = (self.planes + 1) * self.wpr;
        let (rg, rc) = (&self.data[g * stride..(g + 1) * stride], &self.data[c * stride..(c + 1) * stride]);
        let (mut sh, mut sm) = (0u32, 0u32);
        for w in w0..w1 {
            let pa = rg[w] & rc[w]; let mut eq = pa;
            for p in 1..=self.planes { eq &= !(rg[p * self.wpr + w] ^ rc[p * self.wpr + w]); }
            sh += pa.count_ones(); sm += eq.count_ones();
        }
        (sh, sm)
    }
    #[inline(always)]
    fn all_impl(&self, g: usize, c: usize, w: &mut [u32], s: &mut [u32]) { for li in 0..w.len() { let (a, b) = self.count_impl(g, c, li); w[li] = a; s[li] = b; } }
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "popcnt")]
    unsafe fn count_popcnt(&self, g: usize, c: usize, li: usize) -> (u32, u32) { self.count_impl(g, c, li) }
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "popcnt")]
    unsafe fn all_popcnt(&self, g: usize, c: usize, w: &mut [u32], s: &mut [u32]) { self.all_impl(g, c, w, s) }
    /// every locus at once: the two masks of the whole row four words at a time, then the popcounts per locus
    /// (pa, eq: scratch of wpr words). The same integers as count_impl
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "popcnt,avx2")]
    unsafe fn all_avx2(&self, g: usize, c: usize, w: &mut [u32], s: &mut [u32], pa: &mut [u64], eq: &mut [u64]) {
        use std::arch::x86_64::*;
        let wpr = self.wpr; let stride = (self.planes + 1) * wpr; debug_assert!(pa.len() >= wpr && eq.len() >= wpr);
        let (rg, rc) = (self.data.as_ptr().add(g * stride), self.data.as_ptr().add(c * stride)); let mut x = 0;
        while x + 4 <= wpr {
            let p = _mm256_and_si256(_mm256_loadu_si256(rg.add(x) as *const __m256i), _mm256_loadu_si256(rc.add(x) as *const __m256i));
            let mut d = _mm256_setzero_si256();
            for pl in 1..=self.planes { d = _mm256_or_si256(d, _mm256_xor_si256(_mm256_loadu_si256(rg.add(pl * wpr + x) as *const __m256i), _mm256_loadu_si256(rc.add(pl * wpr + x) as *const __m256i))); }
            _mm256_storeu_si256(pa.as_mut_ptr().add(x) as *mut __m256i, p); _mm256_storeu_si256(eq.as_mut_ptr().add(x) as *mut __m256i, _mm256_andnot_si256(d, p)); x += 4;
        }
        while x < wpr { let p = *rg.add(x) & *rc.add(x); let mut d = 0u64; for pl in 1..=self.planes { d |= *rg.add(pl * wpr + x) ^ *rc.add(pl * wpr + x); } pa[x] = p; eq[x] = p & !d; x += 1; }
        for li in 0..w.len() {
            let (mut sh, mut sm) = (0u32, 0u32);
            for x in self.woff[li]..self.woff[li + 1] { sh += pa[x].count_ones(); sm += eq[x].count_ones(); }
            w[li] = sh; s[li] = sm;
        }
    }
    /// shared sites and same residues of two taxa at a locus
    #[inline] pub fn count(&self, g: usize, c: usize, li: usize) -> (u32, u32) {
        #[cfg(target_arch = "x86_64")]
        if self.popcnt { return unsafe { self.count_popcnt(g, c, li) }; }
        self.count_impl(g, c, li)
    }
    /// the same at every locus (pa, eq: scratch of wpr words each)
    #[inline] pub fn all(&self, g: usize, c: usize, w: &mut [u32], s: &mut [u32], pa: &mut [u64], eq: &mut [u64]) {
        #[cfg(target_arch = "x86_64")]
        {
            if self.avx2 { return unsafe { self.all_avx2(g, c, w, s, pa, eq) }; }
            if self.popcnt { return unsafe { self.all_popcnt(g, c, w, s) }; }
        }
        let _ = (pa, eq);
        self.all_impl(g, c, w, s)
    }
    /// bytes of one taxon's planes
    pub fn row_bytes(&self) -> usize { (self.planes + 1) * self.wpr * 8 }
}
pub struct Slab { pub w: Vec<u32>, pub s: Vec<u32> }                // [tri(g, c) * nloc + li], g <= c

impl Core {
    #[inline] pub fn seq(&self, g: usize, li: usize) -> &[i8] { &self.xl[g * self.lt + self.off[li]..g * self.lt + self.off[li + 1]] }
    #[inline] pub fn tri(&self, g: usize, c: usize) -> usize { let (a, b) = if g <= c { (g, c) } else { (c, g) }; a * self.n - a * a.saturating_sub(1) / 2 + (b - a) }
    /// shared sites and same residues of a pair at a locus; a pair with an excluded member counts 0 of each
    #[inline] pub fn at(&self, li: usize, g: usize, c: usize) -> (u32, u32) {
        if let Some(sl) = &self.slab { let k = self.tri(g, c) * self.nloc + li; return (sl.w[k], sl.s[k]); }
        if !(self.keep[g] && self.keep[c]) { return (0, 0); }
        self.bits.count(g, c, li)
    }
    /// the counts of a pair at every locus
    fn pair_counts(&self, g: usize, c: usize, w: &mut [u32], s: &mut [u32], pa: &mut [u64], eq: &mut [u64]) {
        if let Some(sl) = &self.slab { let k = self.tri(g, c) * self.nloc; w.copy_from_slice(&sl.w[k..k + self.nloc]); s.copy_from_slice(&sl.s[k..k + self.nloc]); return; }
        if !(self.keep[g] && self.keep[c]) { w.fill(0); s.fill(0); return; }
        self.bits.all(g, c, w, s, pa, eq);
    }
    fn build_slab(&mut self) {
        let n = self.n; let nloc = self.nloc; let npair = n * (n + 1) / 2;
        let mut w = vec![0u32; npair * nloc]; let mut s = vec![0u32; npair * nloc];
        // row g owns the pairs (g, g..n): a contiguous run of the triangle
        let this: &Core = self;
        let mut wrest: &mut [u32] = &mut w; let mut srest: &mut [u32] = &mut s; let mut parts = Vec::with_capacity(n);
        for g in 0..n { let len = (n - g) * nloc; let (a, b) = wrest.split_at_mut(len); let (c, d) = srest.split_at_mut(len); parts.push((g, a, c)); wrest = b; srest = d; }
        let wpr = this.bits.wpr;
        parts.into_par_iter().for_each_init(|| (vec![0u64; wpr], vec![0u64; wpr]), |(pa, eq), (g, wg, sg)| {
            for c in g..n {
                let k = (c - g) * nloc;
                if !(this.keep[g] && this.keep[c]) { continue; }
                this.bits.all(g, c, &mut wg[k..k + nloc], &mut sg[k..k + nloc], pa, eq);
            }
        });
        self.slab = Some(Slab { w, s });
    }
}

#[derive(Clone, Default)]
pub struct RowList { pub cand: Vec<i64>, pub key: Vec<f64>, pub m: Vec<f64>, pub mfl: Vec<f64>, pub nfin: i64 }

pub struct Prep {
    pub core: Core,
    pub thin: Vec<usize>, pub present: Vec<i64>, pub focal_sites: Vec<u32>,       // focal_sites: n x nloc
    pub plateaus: Vec<f64>, pub rate: Vec<f64>, pub rr: Vec<f64>, pub ok_r: Vec<bool>,
    pub cls: Vec<usize>, pub has_dup: bool,
    pub lists: Vec<RowList>, pub c_len: usize,
    pub cl_win: Vec<i64>, pub cl_dmin: Vec<f64>,                                  // n x nloc
    pub m_full: Option<Vec<f64>>, pub mfl_full: Option<Vec<f64>>,                 // n x n
    pub fin_pairs: i64, pub ext: HashMap<usize, RowList>, pub n_starved: usize,
    pub restart: Option<(String, String)>,
}

#[derive(Clone, Copy, PartialEq)]
pub enum Key { M0, Mj }

pub struct SweepOut {
    pub rows: Vec<RowList>, pub cl_win: Vec<i64>, pub cl_dmin: Vec<f64>, pub ident: Vec<(usize, usize)>,
    pub counters: (i64, i64, i64), pub full: Option<(Vec<f64>, Vec<f64>)>, pub overflow: i64,
}

struct RowRes { list: RowList, cl_win: Vec<i64>, cl_dmin: Vec<f64>, ident: Vec<usize>, counters: (i64, i64, i64), m: Vec<f64>, mfl: Vec<f64>, overflow: i64 }

/// m and the loci voting for one pair from its per-locus counts: the pair kernel's first median
#[inline]
fn median_m(w: &[u32], s: &[u32], rr: Option<&[f64]>, ok_r: Option<&[bool]>, v: &mut [f64], wv: &mut [f64]) -> (f64, usize) {
    let (one, hundred) = (1.0f32, 100.0f32);
    let mut cnt = 0usize;
    for l in 0..w.len() {
        if let Some(ok) = ok_r { if !ok[l] { continue; } }
        let r = match rr { Some(rr) => rr[l], None => 1.0 };
        let wc = w[l] as f32;
        if wc > 0.0 {
            let sc = s[l] as f32;
            let x = (hundred * (one - sc / wc)) as f64 / r;
            if x.is_finite() { v[cnt] = x; wv[cnt] = wc as f64; cnt += 1; }
        }
    }
    (if cnt >= MIN_LOCI_FOR_M && cnt > 0 { wmed_select(v, wv, cnt) } else { f64::NAN }, cnt)
}

/// the pair's resolution floor: the same median of half a mismatch per locus
#[inline]
fn median_floor(w: &[u32], rr: Option<&[f64]>, ok_r: Option<&[bool]>, f: &mut [f64], wf: &mut [f64]) -> f64 {
    let c50f = 50.0f32; let mut cntf = 0usize;
    for l in 0..w.len() {
        if let Some(ok) = ok_r { if !ok[l] { continue; } }
        let r = match rr { Some(rr) => rr[l], None => 1.0 };
        let wc = w[l] as f32;
        if wc > 0.0 {
            let xf = (c50f / wc) as f64 / r;
            if xf.is_finite() { f[cntf] = xf; wf[cntf] = wc as f64; cntf += 1; }
        }
    }
    if cntf >= MIN_LOCI_FOR_M && cntf > 0 { wmed_select(f, wf, cntf) } else { f64::NAN }
}

/// m, its floor and the loci voting for one pair from its per-locus counts
#[inline]
fn pair_medians(w: &[u32], s: &[u32], rr: Option<&[f64]>, ok_r: Option<&[bool]>, v: &mut [f64], wv: &mut [f64], f: &mut [f64], wf: &mut [f64]) -> (f64, f64, usize) {
    let (m, cnt) = median_m(w, s, rr, ok_r, v, wv);
    (m, median_floor(w, rr, ok_r, f, wf), cnt)
}

/// Without computing it: is the pair's m below thr, and how many loci vote? The weighted median is below thr exactly
/// when the loci whose value is below thr carry at least half the weight, so no selection is needed; and a locus's
/// value x / r is compared with thr through the products tlo = thr r (1 - 1e-12) and thi = thr r (1 + 1e-12), the
/// division kept for the values in between (the same answer as comparing the divided value, since r > 0)
#[inline]
fn m_below(w: &[u32], s: &[u32], rr: Option<&[f64]>, ok_r: Option<&[bool]>, thr: f64, tlo: &[f64], thi: &[f64]) -> (bool, usize) {
    let (one, hundred) = (1.0f32, 100.0f32);
    let (mut cnt, mut tot, mut below) = (0usize, 0u64, 0u64);
    for l in 0..w.len() {
        if let Some(ok) = ok_r { if !ok[l] { continue; } }
        if w[l] == 0 { continue; }
        cnt += 1; tot += w[l] as u64;
        let x = (hundred * (one - s[l] as f32 / w[l] as f32)) as f64;
        let lt = if x < tlo[l] { true } else if x > thi[l] { false } else { x / (match rr { Some(rr) => rr[l], None => 1.0 }) < thr };
        if lt { below += w[l] as u64; }
    }
    (2 * below >= tot, cnt)
}

/// per-thread scratch of a pair: its counts at every locus, the medians' buffers, the kernel's masks
struct Scratch { w: Vec<u32>, s: Vec<u32>, v: Vec<f64>, wv: Vec<f64>, f: Vec<f64>, wf: Vec<f64>, pa: Vec<u64>, eq: Vec<u64> }
impl Scratch {
    fn new(core: &Core) -> Scratch {
        let (nloc, wpr) = (core.nloc, core.bits.wpr);
        Scratch { w: vec![0; nloc], s: vec![0; nloc], v: vec![0.0; nloc], wv: vec![0.0; nloc], f: vec![0.0; nloc], wf: vec![0.0; nloc], pa: vec![0; wpr], eq: vec![0; wpr] }
    }
}

/// the closest-match band of one (taxon, locus) as the columns go by: the entries (rate d, vals, column) of the
/// distinct rates within 1e-9 of the lowest so far, ascending, each the (vals, column)-smallest of its rate. One more
/// than BAND_SLOTS is kept, so a band with more distinct rates than slots is known (the overflow the log counts);
/// the result is that of the Python's rule over all columns at once, whatever the order the columns come in
#[derive(Clone, Copy)]
struct Band { n: u8, e: [(f64, f64, i64); BAND_SLOTS + 1] }
const BAND_EMPTY: Band = Band { n: 0, e: [(f64::INFINITY, f64::INFINITY, -1); BAND_SLOTS + 1] };

/// a column whose rate d is within the band (d <= *limit): into the band; the limit follows the minimum
#[inline(never)]
fn band_insert(b: &mut Band, limit: &mut f64, d: f64, mj: f64, c: i64) {
    let v = d + 1e-9 * mj / 100.0; let n = b.n as usize; const CAP: usize = BAND_SLOTS + 1;
    if n == 0 || d < b.e[0].0 {                              // a new minimum: the entries still within 1e-9 of it stay
        let band = d + 1e-9; let old = b.e; let mut k = 1; b.e[0] = (d, v, c);
        for e in old.iter().take(n) { if e.0 <= band && k < CAP { b.e[k] = *e; k += 1; } }
        b.n = k as u8; *limit = band; return;
    }
    for i in 0..n {
        if b.e[i].0 == d { if v < b.e[i].1 || (v == b.e[i].1 && c < b.e[i].2) { b.e[i] = (d, v, c); } return; }
        if b.e[i].0 > d {                                    // a new distinct rate before entry i
            let last = if n < CAP { n } else { CAP - 1 };
            for j in (i + 1..=last).rev() { b.e[j] = b.e[j - 1]; }
            b.e[i] = (d, v, c); b.n = (last + 1) as u8; return;
        }
    }
    if n < CAP { b.e[n] = (d, v, c); b.n += 1; }
}

/// the rows of a chunk against every column. The columns are the outer loop: a column's planes are read from memory
/// once for the whole chunk, whose own rows stay in the cache (the sweep is bound by memory traffic otherwise).
/// Nothing is kept per pair but what the results need: the cmax smallest keys so far and, in sweep 2, the band of
/// the closest match at every locus.
/// lazy: a pair's medians are computed only where their values are used. A row's list needs m only for a pair whose
/// key can be among its cmax smallest: once the list has been cut back, thr is its largest kept key, and a pair
/// whose m is not below thr cannot enter (its key is at least m, and at an equal key the earlier column stays), which
/// m_below decides without the median. The closest match needs m only for a pair that enters a band. Everything else
/// (whether m is finite, the counters, the identical pairs) follows from the counts. In lazy mode the nfin of an M0
/// sweep counts the keys computed, not every eligible key (prepare settles the rows where that matters)
#[allow(clippy::too_many_arguments)]
fn sweep_chunk(core: &Core, rows: &[usize], rr: Option<&[f64]>, ok_r: Option<&[bool]>, key: Key, cmax: usize, fit: f64, want_closest: bool, want_m: bool, lazy: bool, sc: &mut Scratch) -> Vec<RowRes> {
    let (n, nloc) = (core.n, core.nloc); let nr = rows.len();
    struct St { m: Vec<f64>, mfl: Vec<f64>, entries: Vec<(f64, usize, f64, f64)>, nfin: i64, ident: Vec<usize>, counters: (i64, i64, i64), thr: f64, tlo: Vec<f64>, thi: Vec<f64> }
    let mut st: Vec<St> = (0..nr).map(|_| St { m: if want_m { vec![f64::NAN; n] } else { vec![] }, mfl: if want_m { vec![f64::NAN; n] } else { vec![] }, entries: Vec::new(), nfin: 0, ident: vec![], counters: (0, 0, 0),
                                               thr: f64::INFINITY, tlo: if lazy { vec![f64::INFINITY; nloc] } else { vec![] }, thi: if lazy { vec![f64::INFINITY; nloc] } else { vec![] } }).collect();
    let by_key = |a: &(f64, usize, f64, f64), b: &(f64, usize, f64, f64)| a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1));
    let cap = if cmax >= n { usize::MAX } else { 2 * cmax.max(16) };       // the list is cut back to its cmax smallest whenever it reaches this
    let mut bands = if want_closest { vec![BAND_EMPTY; nr * nloc] } else { vec![] };
    let mut limits = if want_closest { vec![f64::INFINITY; nr * nloc] } else { vec![] };
    let mut limits1 = if want_closest { vec![f64::INFINITY; nr * nloc] } else { vec![] };          // limit (1 + 1e-12): the test without a division
    for c in 0..n {
        for (ri, &g) in rows.iter().enumerate() {
            core.pair_counts(g, c, &mut sc.w, &mut sc.s, &mut sc.pa, &mut sc.eq);
            let (w, s) = (&sc.w, &sc.s); let t = &mut st[ri];
            let diag = c == g;
            let tot_same: f64 = s.iter().map(|x| *x as u64).sum::<u64>() as f64; let tot_shared: f64 = w.iter().map(|x| *x as u64).sum::<u64>() as f64;
            let ani = if tot_shared > 0.0 && !diag { 100.0 * tot_same / tot_shared } else { f64::NAN };
            let is_ident = tot_shared > 0.0 && tot_same == tot_shared && !diag;
            let base = (if ani.is_nan() { 100.0 } else { 100.0 - ani }) * 1e-9;
            // mm: the pair's m, or not yet computed (have_m false) in lazy mode; m_fin: whether it is finite
            let (mut mm, mut mf, mut have_m, m_fin, nl, candidate);
            if !lazy {
                let r = if diag { (f64::NAN, f64::NAN, 0) } else { pair_medians(w, s, rr, ok_r, &mut sc.v, &mut sc.wv, &mut sc.f, &mut sc.wf) };
                mm = r.0; mf = r.1; nl = r.2; have_m = true; m_fin = mm.is_finite(); candidate = true;
            } else {
                let (below, cnt) = if diag { (false, 0) } else { m_below(w, s, rr, ok_r, t.thr, &t.tlo, &t.thi) };
                mm = f64::NAN; mf = f64::NAN; have_m = false; nl = cnt; m_fin = cnt >= MIN_LOCI_FOR_M; candidate = below;
            }
            if want_m { t.m[c] = mm; t.mfl[c] = mf; }
            if is_ident && g < c { t.ident.push(c); }
            if g < c && core.keep[g] && core.keep[c] {
                t.counters.1 += 1; if nl < MIN_LOCI_FOR_M { t.counters.0 += 1; } if m_fin { t.counters.2 += 1; }
            }
            let key_fin = m_fin && !is_ident;                    // the pair has a finite mj (and may have a finite m0)
            if key_fin && candidate {
                if !have_m { mm = median_m(w, s, rr, ok_r, &mut sc.v, &mut sc.wv).0; mf = median_floor(w, rr, ok_r, &mut sc.f, &mut sc.wf); have_m = true; }
                let k = match key {
                    Key::M0 => if mm > mf { mm } else { f64::INFINITY },
                    Key::Mj => mm + (base + c as f64 * 1e-12),
                };
                if k.is_finite() {
                    if !lazy || key == Key::M0 { t.nfin += 1; }
                    if k < t.thr {
                        t.entries.push((k, c, mm, mf));
                        if t.entries.len() >= cap {
                            t.entries.select_nth_unstable_by(cmax - 1, by_key); t.entries.truncate(cmax);
                            if lazy {
                                t.thr = t.entries[cmax - 1].0;
                                for l in 0..nloc { let tr = t.thr * (match rr { Some(rr) => rr[l], None => 1.0 }); t.tlo[l] = tr * (1.0 - 1e-12); t.thi[l] = tr * (1.0 + 1e-12); }
                            }
                        }
                    }
                }
            }
            if lazy && key == Key::Mj && key_fin { t.nfin += 1; }
            if want_closest && key_fin {
                let (bd, lm, lm1) = (&mut bands[ri * nloc..(ri + 1) * nloc], &mut limits[ri * nloc..(ri + 1) * nloc], &mut limits1[ri * nloc..(ri + 1) * nloc]);
                for li in 0..nloc {
                    let o = w[li] as f64;
                    if o >= fit {
                        // d = a / om is above the limit for certain when a is above limit (1 + 1e-12) om: no division then
                        let a = o - s[li] as f64 + HALF; let om = o.max(1.0);
                        if a > lm1[li] * om { continue; }
                        let d = a / om;
                        if d <= lm[li] {
                            if !have_m { mm = median_m(w, s, rr, ok_r, &mut sc.v, &mut sc.wv).0; have_m = true; }
                            band_insert(&mut bd[li], &mut lm[li], d, mm + (base + c as f64 * 1e-12), c as i64);
                            lm1[li] = lm[li] * (1.0 + 1e-12);
                        }
                    }
                }
            }
        }
    }
    let mut out = Vec::with_capacity(nr);
    for (ri, mut t) in st.into_iter().enumerate() {
        t.entries.sort_by(by_key); t.entries.truncate(cmax);
        let list = RowList { cand: t.entries.iter().map(|e| e.1 as i64).collect(), key: t.entries.iter().map(|e| e.0).collect(),
                             m: t.entries.iter().map(|e| e.2).collect(), mfl: t.entries.iter().map(|e| e.3).collect(), nfin: t.nfin };
        let mut cl_win = vec![]; let mut cl_dmin = vec![]; let mut overflow = 0i64;
        if want_closest {
            cl_win = vec![-1i64; nloc]; cl_dmin = vec![f64::INFINITY; nloc];
            for li in 0..nloc {
                let b = &bands[ri * nloc + li]; let kept = (b.n as usize).min(BAND_SLOTS);
                if b.n as usize > BAND_SLOTS { overflow += 1; }
                if kept == 0 { continue; }
                // closest_winners: the (vals, column)-smallest of the kept entries
                let mut best = 0usize;
                for j in 1..kept { if b.e[j].1 < b.e[best].1 || (b.e[j].1 == b.e[best].1 && b.e[j].2 < b.e[best].2) { best = j; } }
                cl_win[li] = b.e[best].2; cl_dmin[li] = b.e[0].0;
            }
        }
        out.push(RowRes { list, cl_win, cl_dmin, ident: t.ident, counters: t.counters, m: t.m, mfl: t.mfl, overflow });
    }
    out
}

/// the sweep over the rows of the pair matrix (rows: the given taxa against everyone, for the second passes)
#[allow(clippy::too_many_arguments)]
pub fn sweep(core: &Core, rr: Option<&[f64]>, ok_r: Option<&[bool]>, key: Key, cmax: usize, fit: i64, want_closest: bool, want_full: bool, rows: Option<&[usize]>, stage: Option<&str>) -> SweepOut {
    let (n, nloc) = (core.n, core.nloc);
    let all: Vec<usize> = (0..n).collect(); let rows_v: &[usize] = rows.unwrap_or(&all);
    let prog = crate::log::prog();
    // rows per chunk: about half a megabyte of planes (they stay in the cache while the columns stream by), few enough
    // that every thread has several chunks
    let threads = rayon::current_num_threads().max(1);
    let chunk = (SWEEP_CHUNK_BYTES / core.bits.row_bytes().max(1)).clamp(4, 64).min(rows_v.len().div_ceil(4 * threads)).max(1);
    if let Some(st) = stage {
        log!("pair sweep {}: {} rows of the pair matrix against {} taxa in chunks of {}, {} thread{}, the pair counts {}", st, commas(rows_v.len() as i64), commas(n as i64), chunk, threads,
             if threads > 1 { "s" } else { "" }, if core.slab.is_some() { "kept for the run" } else { "from the residues" });
        prog.start(st, rows_v.len());
    }
    let cmax = cmax.min(n);
    // the medians only where they are used, unless every m is kept (the full matrix) or a locus's rate is so small that
    // a value divided by it could overflow (the shortcut reads the loci that vote from the counts alone)
    let lazy = !want_full && std::env::var_os("LOCUSMASKER_NO_LAZY").is_none()
        && (0..nloc).all(|l| ok_r.map_or(true, |ok| !ok[l]) || rr.map_or(true, |r| r[l].is_finite() && r[l] >= 1e-300));
    let res: Vec<Vec<RowRes>> = rows_v.par_chunks(chunk).map_init(|| Scratch::new(core),
        |sc, ch| { let r = sweep_chunk(core, ch, rr, ok_r, key, cmax, fit as f64, want_closest, want_full, lazy, sc); if stage.is_some() { prog.tick_k(ch.len(), None); } r }).collect();
    if stage.is_some() { prog.finish(); }
    let mut out = SweepOut { rows: Vec::with_capacity(rows_v.len()), cl_win: vec![], cl_dmin: vec![], ident: vec![], counters: (0, 0, 0), full: None, overflow: 0 };
    if want_closest { out.cl_win.reserve(rows_v.len() * nloc); out.cl_dmin.reserve(rows_v.len() * nloc); }
    let mut full = if want_full && rows.is_none() { Some((Vec::with_capacity(n * n), Vec::with_capacity(n * n))) } else { None };
    for (k, r) in res.into_iter().flatten().enumerate() {
        let g = rows_v[k];
        for c in r.ident { out.ident.push((g, c)); }
        out.counters.0 += r.counters.0; out.counters.1 += r.counters.1; out.counters.2 += r.counters.2; out.overflow += r.overflow;
        if want_closest { out.cl_win.extend(r.cl_win); out.cl_dmin.extend(r.cl_dmin); }
        if let Some(f) = full.as_mut() { f.0.extend(r.m); f.1.extend(r.mfl); }
        out.rows.push(r.list);
    }
    out.full = full; out
}

/// each taxon's class of identical genomes: the lowest index among the genomes it is linked to
pub fn classes_from_pairs(n: usize, pairs: &[(usize, usize)]) -> Vec<usize> {
    let mut cls: Vec<usize> = (0..n).collect();
    loop {
        let mut new = cls.clone();
        for &(g, c) in pairs { if cls[c] < new[g] { new[g] = cls[c]; } }
        for &(g, c) in pairs { if cls[g] < new[c] { new[c] = cls[g]; } }
        if new == cls { return cls; }
        cls = new;
    }
}

pub fn encode_alignment(n: usize, l: usize, rows: &[u8], dtype: &str, loci: &Loci) -> (Vec<i8>, Vec<usize>) {
    let alphabet = if dtype == "dna" { crate::input::DNA } else { crate::input::PROTEIN };
    let mut lut = [-1i8; 256];
    for (i, c) in alphabet.bytes().enumerate() { lut[c as usize] = i as i8; lut[c.to_ascii_lowercase() as usize] = i as i8; }
    if dtype == "dna" { lut[b'U' as usize] = 3; lut[b'u' as usize] = 3; }
    let mut letters = [0i64; 256];
    for &b in rows { letters[b as usize] += 1; }
    let (mut n_letters, mut n_in) = (0i64, 0i64);
    for c in (65..91).chain(97..123) { n_letters += letters[c]; if lut[c] >= 0 { n_in += letters[c]; } }
    if n_letters != 0 && (n_in as f64) < 0.95 * n_letters as f64 {
        log!("WARNING: {} of the {} letters ({:.1}%) are ambiguity codes and count as missing data", commas(n_letters - n_in), commas(n_letters), 100.0 * (n_letters - n_in) as f64 / n_letters as f64);
    }
    let mut off = Vec::with_capacity(loci.len() + 1); let mut acc = 0; for (_, st) in loci { off.push(acc); acc += st.len(); } off.push(acc);
    let lt = acc; let mut xl = vec![-1i8; n * lt];
    xl.par_chunks_mut(lt.max(1)).enumerate().for_each(|(g, row)| {
        if lt == 0 { return; }
        let src = &rows[g * l..(g + 1) * l]; let mut k = 0;
        for (_, st) in loci { for &site in st { row[k] = lut[src[site] as usize]; k += 1; } }
    });
    (xl, off)
}

pub fn locus_plateaus(core: &Core) -> Vec<f64> {
    (0..core.nloc).into_par_iter().map(|li| {
        let mut cnt = vec![0u64; core.a_];
        for g in 0..core.n { for &x in core.seq(g, li) { if x >= 0 { cnt[x as usize] += 1; } } }
        let cntf: Vec<f64> = cnt.iter().map(|c| *c as f64).collect();
        let tot = pymax(pairwise_sum(&cntf), 1.0);
        let f2: Vec<f64> = cntf.iter().map(|c| { let f = c / tot; f * f }).collect();
        pymax(100.0 * (1.0 - pairwise_sum(&f2)), 1.0)
    }).collect()
}

impl Prep {
    #[inline] pub fn n(&self) -> usize { self.core.n }
    /// m and its floor for a pair: from the full matrix when it is kept, else by the kernel on the pair's counts (NaN for a taxon with itself)
    pub fn pair_m(&self, g: usize, c: usize) -> (f64, f64) {
        if let (Some(m), Some(f)) = (&self.m_full, &self.mfl_full) { return (m[g * self.core.n + c], f[g * self.core.n + c]); }
        if g == c { return (f64::NAN, f64::NAN); }
        thread_local! { static SC: std::cell::RefCell<Option<Scratch>> = const { std::cell::RefCell::new(None) }; }
        SC.with(|cell| {
            let mut slot = cell.borrow_mut();
            if slot.as_ref().map_or(true, |sc| sc.w.len() != self.core.nloc || sc.pa.len() != self.core.bits.wpr) { *slot = Some(Scratch::new(&self.core)); }
            let sc = slot.as_mut().unwrap();
            self.core.pair_counts(g, c, &mut sc.w, &mut sc.s, &mut sc.pa, &mut sc.eq);
            let (m, mf, _) = pair_medians(&sc.w, &sc.s, Some(&self.rr), Some(&self.ok_r), &mut sc.v, &mut sc.wv, &mut sc.f, &mut sc.wf);
            (m, mf)
        })
    }
    /// m and its floor for a comparator pair, from the focal taxon's list (or its extended row, or the full matrix)
    pub fn pair_m_lists(&self, g: usize, c: usize) -> (f64, f64) {
        if let (Some(m), Some(f)) = (&self.m_full, &self.mfl_full) { return (m[g * self.core.n + c], f[g * self.core.n + c]); }
        let l = &self.lists[g];
        if let Some(p) = l.cand.iter().position(|x| *x == c as i64) { return (l.m[p], l.mfl[p]); }
        if let Some(e) = self.ext.get(&g) { if let Some(p) = e.cand.iter().position(|x| *x == c as i64) { return (e.m[p], e.mfl[p]); } }
        (f64::NAN, f64::NAN)
    }
    /// the candidates of a taxon in order of pair divergence: its extended row when it has one, else its list
    pub fn row_cands(&self, g: usize) -> &[i64] { match self.ext.get(&g) { Some(e) => &e.cand, None => &self.lists[g].cand } }
}

pub fn prepare(a: &Args, names: &[String], n: usize, l: usize, rows: &[u8], loci: &Loci, dtype: &str, k_eff: usize, t0: std::time::Instant) -> Prep {
    let a_ = if dtype == "dna" { 4 } else { 20 }; let nloc = loci.len();
    let (xl, off) = encode_alignment(n, l, rows, dtype, loci);
    let lt = *off.last().unwrap();
    let bits = Bits::build(n, &xl, lt, &off, a_);
    let mut core = Core { n, nloc, lt, xl, off, keep: vec![true; n], a_, slab: None, bits };
    let mut focal_sites = vec![0u32; n * nloc];
    focal_sites.par_chunks_mut(nloc).enumerate().for_each(|(g, fs)| { for li in 0..nloc { fs[li] = core.seq(g, li).iter().filter(|x| **x >= 0).count() as u32; } });
    let present: Vec<i64> = (0..n).map(|g| focal_sites[g * nloc..(g + 1) * nloc].iter().filter(|x| **x as i64 >= a.fit_overlap).count() as i64).collect();
    let thin: Vec<usize> = (0..n).filter(|&g| present[g] < MIN_LOCI_PER_TAXON).collect();
    for &g in &thin { core.keep[g] = false; }
    if !thin.is_empty() {
        log!("{} taxa with fewer than {} loci present are excluded from the test and left unmasked: {}{}", thin.len(), MIN_LOCI_PER_TAXON,
             thin.iter().take(20).map(|&i| names[i].as_str()).collect::<Vec<_>>().join(", "),
             if thin.len() > 20 { format!(" and {} more (the report lists them)", thin.len() - 20) } else { String::new() });
    }
    let plateaus = locus_plateaus(&core);
    let longest = loci.iter().map(|x| x.1.len()).max().unwrap_or(0);
    assert!(longest < (1 << 24), "a locus of 2^24 sites or more: the float32 products would not be exact");
    let c = a.candidates.unwrap_or(CANDIDATES).max(1);
    let full = n <= a.full_m_max_taxa;
    log!("pair sweep: rows of the pair matrix in {} thread{}, {} candidates per taxon, the full m matrix {}", a.threads, if a.threads > 1 { "s" } else { "" }, c, if full { "kept" } else { "not kept" });
    let kn = k_eff.min(n - 1);
    if n <= a.sweep_block_min {
        core.build_slab();
        let bytes = core.slab.as_ref().map(|s| 8 * s.w.len()).unwrap_or(0);
        log!("the pair counts of every locus ({}) are kept for the run", crate::fmem(bytes as u64));
    }
    let header = crate::restart::header(a, n, nloc, l, dtype, kn, c, &thin, full);
    let mut pair_names = vec!["ident", "rate", "ok_r", "cand", "key", "m", "mfl", "nfin", "cl_win", "cl_dmin", "counters"]; if full { pair_names.extend(["m_full", "mfl_full"]); }
    let pt2 = crate::restart::read(a, "pairs", &header, &pair_names);
    let pt1 = if pt2.is_none() { crate::restart::read(a, "rates", &header, &["ident", "rate", "ok_r"]) } else { None };
    let c_len = c.min(n);
    let mut restart = None; let mut n_starved = 0usize;
    let (ident, mut rate, ok_r): (Vec<(usize, usize)>, Vec<f64>, Vec<bool>);
    let from_point = pt2.as_ref().or(pt1.as_ref());
    if let Some(pt) = from_point {
        restart = Some(((if pt2.is_some() { "the pair sweeps" } else { "the locus rates" }).to_string(), pt.written.clone()));
        ident = pt.i64s("ident").chunks_exact(2).map(|x| (x[0] as usize, x[1] as usize)).collect(); rate = pt.f64s("rate"); ok_r = pt.bools("ok_r");
    } else {
        // ---- sweep 1: the raw divergences, for the rate-vote neighbours
        let s1 = sweep(&core, None, None, Key::M0, c, a.fit_overlap, false, false, None, Some("pairs0"));
        ident = s1.ident.clone();
        let cls0 = classes_from_pairs(n, &ident);
        let rep: Vec<bool> = (0..n).map(|i| cls0[i] == i).collect();
        // the locus rates, voted by each taxon's nearest resolved non-identical neighbours
        let neighbours = |list: &RowList, g: usize, out: &mut Vec<(usize, usize, f64)>| -> bool {
            let mut picked = 0;
            for (j, &cd) in list.cand.iter().enumerate() {
                if picked == kn { break; }
                if rep[cd as usize] { out.push((g, cd as usize, npmax(list.m[j], list.mfl[j]))); picked += 1; }
            }
            picked < kn
        };
        // a taxon whose list ran out before Kn representatives although it has more candidates than the list holds gets
        // its full row. The sweep does not count every candidate (it skips the pairs that cannot enter the list), so
        // the rows that are short with a full list are swept in full: more entries than the list means it ran out
        let mut votes: Vec<(usize, usize, f64)> = vec![]; let mut short_full: Vec<usize> = vec![]; let mut spans: Vec<(usize, usize)> = vec![];
        for g in 0..n {
            let before = votes.len();
            let short = neighbours(&s1.rows[g], g, &mut votes);
            if short && s1.rows[g].cand.len() >= c_len { short_full.push(g); spans.push((before, votes.len())); }
        }
        if !short_full.is_empty() {
            let e = sweep(&core, None, None, Key::M0, n, a.fit_overlap, false, false, Some(&short_full), None);
            let mut drop = vec![false; votes.len()]; let mut extra: Vec<(usize, usize, f64)> = vec![];
            for (k, &g) in short_full.iter().enumerate() {
                if e.rows[k].cand.len() > c_len { n_starved += 1; for d in drop[spans[k].0..spans[k].1].iter_mut() { *d = true; } neighbours(&e.rows[k], g, &mut extra); }
            }
            if n_starved > 0 {
                log!("rate neighbours: {} taxa ran out of candidates (copies of one genome fill the list): full rows computed", n_starved);
                votes = votes.into_iter().zip(drop).filter(|(_, d)| !d).map(|(v, _)| v).collect(); votes.extend(extra);
            }
        }
        drop(s1);
        let prog = crate::log::prog();
        prog.start("rates", nloc);
        let fit_f = a.fit_overlap as f32;
        rate = (0..nloc).into_par_iter().map(|li| {
            let p_ = plateaus[li];
            let corr = |p: f64| -> f64 { -p_ * (1.0 - npmin(p, 0.98 * p_) / p_).ln() };
            let mut vals: Vec<f64> = Vec::with_capacity(votes.len());
            for &(g, cc, m0) in &votes {
                let (o, s) = core.at(li, g, cc); let (o_, s_) = (o as f32, s as f32);
                let d_ = if o_ > 0.0 { 100.0f32 * (1.0f32 - s_ / o_.max(1.0)) } else { f32::NAN };
                if d_.is_finite() && o_ >= fit_f { vals.push(corr((d_ + 50.0f32 / o_) as f64) / corr(m0)); }
            }
            prog.tick();
            if vals.len() < MIN_FIT_PAIRS { f64::NAN } else { median(&vals) }
        }).collect();
        prog.finish();
        let lsites: Vec<f64> = loci.iter().map(|x| x.1.len() as f64).collect();
        ok_r = rate.iter().map(|r| r.is_finite() && *r > 0.0).collect();
        if ok_r.iter().filter(|x| **x).count() >= 3 {
            let (rv, lv): (Vec<f64>, Vec<f64>) = (0..nloc).filter(|&i| ok_r[i]).map(|i| (rate[i], lsites[i])).unzip();
            let scale = wmedian_sorted(&rv, &lv);
            let div = pymax(scale, 1e-6);
            for r in rate.iter_mut() { *r /= div; }
        }
    }
    let cls = classes_from_pairs(n, &ident); let has_dup = !ident.is_empty();
    if has_dup {
        log!("{} genomes are identical to another at every shared site: an identical genome is never a comparator, and no two comparators are identical to each other", (0..n).filter(|&i| cls[i] != i).count());
    }
    if restart.is_none() { log!("per-locus identities done, {:.0}s", t0.elapsed().as_secs_f64()); }
    let rr: Vec<f64> = (0..nloc).map(|i| if ok_r[i] { rate[i] } else { f64::NAN }).collect();
    let ident_flat: Vec<i64> = ident.iter().flat_map(|x| [x.0 as i64, x.1 as i64]).collect();
    use crate::restart::Arr;
    if restart.is_none() { crate::restart::write(a, "rates", &header, &[("ident", Arr::I64(&ident_flat, vec![ident.len(), 2])), ("rate", Arr::F64(&rate, vec![nloc])), ("ok_r", Arr::Bool(&ok_r))]); }
    let s2 = if let Some(pt) = &pt2 {
        let (cand, key, m, mfl, nfin) = (pt.i64s("cand"), pt.f64s("key"), pt.f64s("m"), pt.f64s("mfl"), pt.i64s("nfin"));
        let rows = (0..n).map(|g| { let k = (0..c_len).take_while(|&j| cand[g * c_len + j] >= 0).count(); let r = g * c_len..g * c_len + k;
            RowList { cand: cand[r.clone()].to_vec(), key: key[r.clone()].to_vec(), m: m[r.clone()].to_vec(), mfl: mfl[r].to_vec(), nfin: nfin[g] } }).collect();
        let cn = pt.i64s("counters");
        SweepOut { rows, cl_win: pt.i64s("cl_win"), cl_dmin: pt.f64s("cl_dmin"), ident: vec![], counters: (cn[0], cn[1], cn[2]), full: if full { Some((pt.f64s("m_full"), pt.f64s("mfl_full"))) } else { None }, overflow: 0 }
    } else {
        // ---- sweep 2: the rate-normalised divergences, the lists, the closest match at every locus
        let s2 = sweep(&core, Some(&rr), Some(&ok_r), Key::Mj, c, a.fit_overlap, true, full, None, Some("pairs"));
        let pad_i = |f: &dyn Fn(&RowList) -> &Vec<i64>| -> Vec<i64> { let mut v = vec![-1i64; n * c_len]; for g in 0..n { let x = f(&s2.rows[g]); v[g * c_len..g * c_len + x.len()].copy_from_slice(x); } v };
        let pad_f = |f: &dyn Fn(&RowList) -> &Vec<f64>, fill: f64| -> Vec<f64> { let mut v = vec![fill; n * c_len]; for g in 0..n { let x = f(&s2.rows[g]); v[g * c_len..g * c_len + x.len()].copy_from_slice(x); } v };
        let (cand, key, m, mfl) = (pad_i(&|r| &r.cand), pad_f(&|r| &r.key, f64::INFINITY), pad_f(&|r| &r.m, f64::NAN), pad_f(&|r| &r.mfl, f64::NAN));
        let nfin: Vec<i64> = s2.rows.iter().map(|r| r.nfin).collect(); let cn = [s2.counters.0, s2.counters.1, s2.counters.2];
        let mut arrays = vec![("ident", Arr::I64(&ident_flat, vec![ident.len(), 2])), ("rate", Arr::F64(&rate, vec![nloc])), ("ok_r", Arr::Bool(&ok_r)),
                              ("cand", Arr::I64(&cand, vec![n, c_len])), ("key", Arr::F64(&key, vec![n, c_len])), ("m", Arr::F64(&m, vec![n, c_len])), ("mfl", Arr::F64(&mfl, vec![n, c_len])),
                              ("nfin", Arr::I64(&nfin, vec![n])), ("cl_win", Arr::I64(&s2.cl_win, vec![n, nloc])), ("cl_dmin", Arr::F64(&s2.cl_dmin, vec![n, nloc])), ("counters", Arr::I64(&cn, vec![3]))];
        if let Some((mf, ff)) = &s2.full { arrays.push(("m_full", Arr::F64(mf, vec![n, n]))); arrays.push(("mfl_full", Arr::F64(ff, vec![n, n]))); }
        crate::restart::write(a, "pairs", &header, &arrays);
        s2
    };
    let (nl_low, npk, fin_pairs) = s2.counters; let fb = if npk != 0 { nl_low as f64 / npk as f64 } else { 0.0 };
    let n_ok = ok_r.iter().filter(|x| **x).count();
    log!("pair divergences done, {:.0}s{}", t0.elapsed().as_secs_f64(),
         if n_ok == 0 { "; note: no locus has a rate (no nearest-neighbour pairs resolved beyond their resolution floor), so no pair has a divergence".to_string() }
         else if fb > 0.05 { format!("; note: {:.0}% of taxon pairs overlap in <{} loci with a rate and are never compared", 100.0 * fb, MIN_LOCI_FOR_M) } else { String::new() });
    if s2.overflow != 0 { log!("note: {} cells had more distinct closest-match rates within 1e-9 of the minimum than {} slots", commas(s2.overflow), BAND_SLOTS); }
    let (m_full, mfl_full) = match s2.full { Some((m, f)) => (Some(m), Some(f)), None => (None, None) };
    Prep { core, thin, present, focal_sites, plateaus, rate, rr, ok_r, cls, has_dup, lists: s2.rows, c_len, cl_win: s2.cl_win, cl_dmin: s2.cl_dmin,
           m_full, mfl_full, fin_pairs, ext: HashMap::new(), n_starved, restart }
}

// ------------------------------------------------------------------ the broad-sample pools
pub struct BsPools { pub k: usize, pub pn: usize, pub pool: Vec<i32>, pub mstar: Vec<f64> }      // pool, mstar: n x pn

pub fn bs_keys_taxa(i: usize) -> u64 { sm64((i as u64 + 1) ^ BS_SEED) }
pub fn bs_keys_cell(g: usize, li: usize) -> u64 {
    let mix = (li as u64 + 1).wrapping_mul(0x9E3779B97F4A7C15);
    sm64(sm64((g as u64).wrapping_add(1)) ^ mix ^ BS_SEED)
}

pub fn bs_prepare(p: &Prep, k_eff: usize) -> Option<BsPools> {
    let n = p.n();
    let k = (BS_PER_COMP * k_eff).min(n - 1); let pn = (BS_POOL * k_eff).min(n - 1);
    if (k as i64) < BS_MIN { return None; }
    let prog = crate::log::prog(); let on_demand = p.m_full.is_none();
    if on_demand { prog.start("broad-sample pool m", n); }
    let rows: Vec<(Vec<i32>, Vec<f64>)> = (0..n).into_par_iter().map(|g| {
        let mut pool: Vec<i64> = if pn >= n - 1 { (0..(n - 1) as i64).collect() } else { hash_draw(bs_keys_taxa(g), n - 1, pn) };
        for v in pool.iter_mut() { if *v >= g as i64 { *v += 1; } }
        let ms: Vec<f64> = pool.iter().map(|&c| { let (m_, f_) = p.pair_m(g, c as usize); if f_.is_finite() { m_.max(f_) } else { m_ } }).collect();
        if on_demand { prog.tick(); }
        (pool.iter().map(|x| *x as i32).collect(), ms)
    }).collect();
    if on_demand { prog.finish(); }
    let mut pool = Vec::with_capacity(n * pn); let mut mstar = Vec::with_capacity(n * pn);
    for (a, b) in rows { pool.extend(a); mstar.extend(b); }
    Some(BsPools { k, pn, pool, mstar })
}
