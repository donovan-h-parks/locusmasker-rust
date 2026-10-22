#!/usr/bin/env python3
"""Per-locus outlier detection and masking in a multi-locus alignment.

    python locus_outliers.py --aln alignment.nex --out DIR
    python locus_outliers.py --aln alignment.fasta --partitions partitions.txt --out DIR
    python locus_outliers.py --aln alignment.phy   --partitions partitions.txt --out DIR

Input: a NEXUS alignment whose SETS block defines the loci, one charset per locus (read with
Biopython's Bio.Nexus), or a FASTA / PHYLIP / NEXUS alignment without charsets plus a RAxML-style
partition file ("DNA, locus1 = 1-100, 250-300" and "WAG, locus2 = 101-249\\3" forms; the model
prefix is optional). A NEXUS file with charsets and a partition file together is an error: give
one or the other. Sites are 1-based and inclusive. Data type is taken from the NEXUS file when
it states one, otherwise inferred from the residues; DNA is analysed on A, C, G, T only (U read
as T), protein on the 20 standard amino acids; everything else is missing. At least 3 loci and
4 taxa are required.

Method (--test z, the default). Every pair of taxa gets a measured divergence m: the median,
across loci, of the pair's rate-normalised divergence — one vote per overlapping locus, so one
bad locus cannot move its own baseline no matter how long it is. Each taxon is compared, locus by locus, with its
--n-comp closest comparators by m (odd, >= 3, auto-shrunk on tiny datasets). The expected
divergence of a pair in a locus is m times a ratio curve: a running median of
(locus divergence / m) against log10(m), flat outside its support (constant relative rate),
capped at the locus's random-identity plateau. A pair fails when all three hold: its binomial z — observed
mismatches against expected, with variance o*q*(1-q) on o shared sites — exceeds the locus
median z by --k robust SDs of z (floored at 1); its divergence is at least --min-ratio times
expected; and the absolute excess is at least --min-excess mismatches. A cell (taxon, locus)
is masked when a majority of its comparators fail; cells with fewer comparators than requested
that still fail a majority are reported as ambiguous and masked too. Masking blanks the cell's
whole stretch of the locus, ambiguity codes included.

--test band is the legacy rule kept for comparison: residual = locus identity minus a trend of
locus identity over alignment-wide identity; a pair fails when the residual is more than
max(--min-band, --k robust SDs) below the locus median. It also keeps the legacy skip rule
(a locus needs 2 x --bin-min testable cells).

Outputs (in --out; <stem> is the alignment file's name without extension)
  <stem>.masked.<ext>  the alignment in its input format with every outlier cell replaced by gaps
  masked.csv           every masked cell with the evidence
  cells.tsv            every tested cell (flagged or not)
  pairs.tsv            every comparison
  per_locus.tsv        one row per locus, including skipped loci (status column)
  trends.tsv           the locus expectation on a grid, with the effect-size cutoff (z mode)
  ani.npy, m.npy       taxon x taxon alignment-wide identity (%) and pair divergence m (points)
  taxa.txt             the row order of the matrices
Options
  --test z|band  --n-comp 5 (odd, >= 3)  --fit-overlap 20  --k 5  --min-ratio 2 --min-excess 4
  --min-band 25 (band mode)  --bin 0.1  --bin-min 40  --seed 1
"""
import argparse, csv, os, re, sys, time
import numpy as np

PROTEIN = "ACDEFGHIKLMNPQRSTVWY"
DNA = "ACGT"


# ------------------------------------------------------------------ input
def read_fasta(txt):
    names, seqs, cur = [], [], []
    for line in txt.splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur))
            names.append(line[1:].split()[0]); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return names, seqs


def read_phylip(txt):
    lines = [l for l in txt.splitlines() if l.strip()]
    n, L = map(int, lines[0].split()); names, seqs = [], []
    for line in lines[1:]:
        p = line.split(None, 1)
        if len(p) == 2 and len(names) < n: names.append(p[0]); seqs.append(p[1].replace(" ", ""))
        elif names: seqs[-1] += line.replace(" ", "")     # interleaved continuation
    if any(len(s) != L for s in seqs):                      # interleaved with names only in the first block
        blocks = [l for l in lines[1:]]; names, seqs = [], []
        for i, line in enumerate(blocks):
            p = line.split(None, 1)
            if i < n: names.append(p[0]); seqs.append(p[1].replace(" ", "") if len(p) == 2 else "")
            else: seqs[i % n] += line.replace(" ", "")
    assert len(names) == n and all(len(s) == L for s in seqs), "cannot parse PHYLIP"
    return names, seqs


def parse_raxml_partitions(path, L):
    """name -> 0-based site index array; ranges 'a-b', 'a-b\\s' (step) and comma lists"""
    loci = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line or "=" not in line: continue
        lhs, rhs = line.split("=", 1)
        name = lhs.split(",")[-1].strip()
        sites = []
        for part in rhs.replace(";", "").split(","):
            part = part.strip()
            if not part: continue
            m = re.match(r"^(\d+)\s*-\s*(\d+)(?:\s*\\\s*(\d+))?$", part)
            if m:
                a, b, s = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
                sites.extend(range(a - 1, b, s))
            elif re.match(r"^\d+$", part):
                sites.append(int(part) - 1)
            else:
                sys.exit(f"cannot parse partition range {part!r} in {line!r}")
        arr = np.array(sorted(set(sites)), dtype=int)
        if arr.max() >= L: sys.exit(f"partition {name} runs past the alignment ({arr.max() + 1} > {L})")
        loci.append((name, arr))
    return loci


def infer_dtype(seqs):
    """DNA when at least 90% of the non-missing characters are A, C, G, T or U; else protein"""
    res = "".join(s.upper() for s in seqs[: min(len(seqs), 100)]); res = re.sub(r"[-?XN.~*]", "", res)
    if not res: return "protein"
    return "dna" if sum(res.count(c) for c in "ACGTU") / len(res) >= 0.9 else "protein"


def read_input(path, partitions):
    txt = open(path).read(); head = txt.lstrip()[:7].upper()
    nexus = None; dtype = None
    if head.startswith("#NEXUS"):
        from Bio.Nexus import Nexus
        nexus = Nexus.Nexus(path)
        names = list(nexus.taxlabels); seqs = [str(nexus.matrix[t]) for t in names]; fmt = "nexus"
        if nexus.charsets and partitions:
            sys.exit("the NEXUS file has charsets and a partition file was given: use one or the other, not both")
        if nexus.charsets:
            loci = [(k, np.array(sorted(v), dtype=int)) for k, v in nexus.charsets.items()]
        elif partitions:
            loci = parse_raxml_partitions(partitions, len(seqs[0]))
        else:
            sys.exit("the NEXUS file has no charsets: add a SETS block or give --partitions")
        dt = (nexus.datatype or "").lower()
        if dt in ("dna", "rna", "nucleotide"): dtype = "dna"
        elif dt == "protein": dtype = "protein"
    else:
        if not partitions: sys.exit("a FASTA or PHYLIP alignment needs --partitions")
        names, seqs = read_fasta(txt) if head.startswith(">") else read_phylip(txt)
        fmt = "fasta" if head.startswith(">") else "phylip"
        loci = parse_raxml_partitions(partitions, len(seqs[0]))
    L = len(seqs[0]); assert all(len(s) == L for s in seqs), "sequences differ in length"
    if not loci: sys.exit("no loci defined")
    if dtype is None: dtype = infer_dtype(seqs)
    return names, seqs, loci, fmt, dtype, nexus


# ------------------------------------------------------------------ output
def write_masked(path, out_dir, fmt, names, seqs, nexus):
    stem, ext = os.path.splitext(os.path.basename(path)); out = os.path.join(out_dir, f"{stem}.masked{ext or '.' + fmt}")
    if fmt == "nexus":
        from Bio.Seq import Seq
        for t, s in zip(names, seqs): nexus.matrix[t] = Seq(s)
        nexus.write_nexus_data(filename=out, interleave=False, append_sets=True)
    elif fmt == "fasta":
        with open(out, "w") as fh:
            for t, s in zip(names, seqs): fh.write(f">{t}\n{s}\n")
    else:
        with open(out, "w") as fh:
            fh.write(f"{len(names)} {len(seqs[0])}\n")
            for t, s in zip(names, seqs): fh.write(f"{t} {s}\n")
    return out


# ------------------------------------------------------------------ shared machinery
def trend_fit(xf, yf, width, nmin):
    """running-median centres of yf along xf: windows of `width`, widened x1.5 until they hold
    nmin points; returns (centres, medians); a single centre when windows cannot form"""
    order = np.argsort(xf); xs, ys = xf[order], yf[order]
    lo, hi = xs[0], xs[-1]
    cx, cy = [], []
    for c in np.arange(lo + width / 2, hi, width):
        h = width / 2
        while True:
            sel = (xs >= c - h) & (xs <= c + h)
            if sel.sum() >= nmin or h > (hi - lo): break
            h *= 1.5
        if sel.sum() >= max(5, nmin // 4): cx.append(c); cy.append(np.median(ys[sel]))
    if not cx: cx, cy = [float(np.median(xs))], [float(np.median(ys))]
    return np.array(cx), np.array(cy)


def pair_median(vals, W, fit_overlap, fallback):
    """weighted median of vals across loci (axis 0): weights are shared sites, but no locus may
    carry more than ~a third of a pair's total weight (iterative cap), so a single bad locus can
    never reach the median however long it is, while short noisy fragments still count little.
    Pairs with fewer than 3 contributing loci get `fallback`."""
    fin = np.isfinite(vals)
    w = np.where(fin, W, 0.0)
    t = 0.34 * w.sum(0)
    for _ in range(8):
        t = 0.34 * np.minimum(w, t[None, :]).sum(0)
    w = np.minimum(w, t[None, :])
    v = wmedian_stack(vals, w)
    nl = (w > 0).sum(0)
    return np.where(nl >= 3, v, fallback), nl


def wmedian_stack(vals, wts):
    """weighted median of vals over axis 0, nan-aware; (L, P) -> (P,)"""
    ok = np.isfinite(vals) & (wts > 0)
    v = np.where(ok, vals, np.inf); w = np.where(ok, wts, 0.0)
    order = np.argsort(v, axis=0)
    v_s = np.take_along_axis(v, order, 0); w_s = np.take_along_axis(w, order, 0)
    cw = np.cumsum(w_s, axis=0, dtype=np.float64); tot = cw[-1]
    idx = np.minimum((cw < tot[None, :] / 2).sum(0), v.shape[0] - 1)
    out = np.take_along_axis(v_s, idx[None, :], 0)[0]
    return np.where(tot > 0, out, np.nan)


def top_comparators(score, kk, nres):
    """rank comparators by score (higher = closer); returns top indices, validity, counts"""
    top = np.argpartition(-score, kk - 1, axis=1)[:, :kk]
    top = np.take_along_axis(top, np.argsort(-np.take_along_axis(score, top, 1), axis=1), 1)
    valid = np.isfinite(np.take_along_axis(score, top, 1)) & (nres > 0)[:, None]
    return top, valid, valid.sum(1)


def flatten_pairs(idx, top, valid, kk):
    gi = np.repeat(idx, kk); ci = top[idx].ravel(); vi = valid[idx].ravel(); rk = np.tile(np.arange(kk), len(idx))
    return gi[vi], ci[vi], rk[vi]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True); ap.add_argument("--partitions"); ap.add_argument("--out", required=True)
    ap.add_argument("--test", choices=("z", "band"), default="z")
    ap.add_argument("--n-comp", type=int, default=5); ap.add_argument("--fit-overlap", type=int, default=20)
    ap.add_argument("--k", type=float, default=5.0)
    ap.add_argument("--min-ratio", type=float, default=2.0); ap.add_argument("--min-excess", type=float, default=4.0)
    ap.add_argument("--min-band", type=float, default=25.0)
    ap.add_argument("--bin", type=float, default=0.1); ap.add_argument("--bin-min", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); K = a.n_comp
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    names, seqs, loci, fmt, dtype, nexus = read_input(a.aln, a.partitions)
    if len(loci) < 3:
        sys.exit(f"only {len(loci)} locus defined: each locus is judged against the rest of the "
                 f"alignment, which needs at least 3 loci to mean anything")
    n, L = len(names), len(seqs[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    if n < 4: sys.exit(f"only {n} taxa: the comparator test needs at least 4")
    K_eff = min(K, n - 1); K_eff -= 1 - K_eff % 2                    # largest odd <= min(K, n-1)
    M_eff = (K_eff + 1) // 2                                          # fails that make a majority
    if K_eff < K: print(f"n-comp shrunk to {K_eff} ({n} taxa)", file=sys.stderr)
    print(f"{n} taxa, {L} sites, {len(loci)} loci, {dtype}, {fmt}, test={a.test}", file=sys.stderr)
    open(f"{a.out}/taxa.txt", "w").write("\n".join(names) + "\n")

    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = lut[np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L)]

    # ---- pass 1: per-locus shared/same counts; alignment-wide identity falls out of the sums
    nloc = len(loci)
    S = np.empty((nloc, n, n), np.float32); W = np.empty((nloc, n, n), np.float32)
    plateaus = np.empty(nloc)
    for li, (name, sites) in enumerate(loci):
        Xm = X[:, sites]; Pm = (Xm >= 0).astype(np.float32)
        sh = Pm @ Pm.T; sm = np.zeros((n, n), np.float32)
        for i in range(len(alphabet)):
            Ai = (Xm == i).astype(np.float32); sm += Ai @ Ai.T
        S[li], W[li] = sm, sh
        cnt = np.bincount(Xm[Xm >= 0], minlength=len(alphabet)).astype(float)
        f = cnt / max(cnt.sum(), 1.0)
        plateaus[li] = max(100.0 * (1.0 - float((f ** 2).sum())), 1.0)
    tot_same = S.sum(0); tot_shared = W.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(tot_shared > 0, 100.0 * tot_same / tot_shared, np.nan)
    np.fill_diagonal(ani, np.nan); np.save(f"{a.out}/ani.npy", ani)
    print(f"per-locus identities done, {time.time() - t0:.0f}s", file=sys.stderr)

    # ---- pair divergence m: weighted median across loci of rate-normalised divergence
    with np.errstate(invalid="ignore", divide="ignore"):
        Dv = np.where(W > 0, 100.0 * (1.0 - S / np.maximum(W, 1)), np.nan)
    Wm = W.copy()
    Wm[:, np.arange(n), np.arange(n)] = 0.0                            # self pairs carry no vote
    pooled = 100.0 - ani
    m0, _ = pair_median(Dv.reshape(nloc, -1), Wm.reshape(nloc, -1), a.fit_overlap,
                        pooled.reshape(-1))
    m0 = m0.reshape(n, n)
    iu = np.triu_indices(n, 1)
    m0c = np.maximum(m0[iu], 0.05)
    rate = np.full(nloc, np.nan)
    for li in range(nloc):
        d, w = Dv[li][iu], W[li][iu]
        val = np.isfinite(d) & (w >= a.fit_overlap)
        if val.sum() >= 8: rate[li] = float(np.median(d[val] / m0c[val]))
    lsites = np.array([len(s) for _, s in loci], dtype=float)
    ok_r = np.isfinite(rate) & (rate > 0)
    if ok_r.sum() >= 3:
        scale = wmedian_stack(rate[ok_r, None], lsites[ok_r, None])[0]
        rate = rate / max(scale, 1e-6)
    rr = np.where(ok_r, rate, np.nan)
    m, nl = pair_median((Dv / rr[:, None, None]).reshape(nloc, -1),
                        np.where(ok_r[:, None, None], Wm, 0.0).reshape(nloc, -1),
                        a.fit_overlap, pooled.reshape(-1))
    m = m.reshape(n, n); nl = nl.reshape(n, n)
    np.fill_diagonal(m, np.nan); np.save(f"{a.out}/m.npy", m)
    fb = float((nl[np.triu_indices(n, 1)] < 3).mean())
    print(f"pair divergences done, {time.time() - t0:.0f}s"
          + (f"; WARNING: {fb:.0%} of taxon pairs overlap in <3 loci and use pooled identity "
             f"instead of the robust pair median — expectations for those pairs can be dragged "
             f"by a bad locus" if fb > 0.05 else ""), file=sys.stderr)

    rng = np.random.default_rng(a.seed)
    jit = rng.random((n, n)) * 1e-9
    masked = [bytearray(s.encode()) for s in seqs]; nmasked = 0
    mcsv = open(f"{a.out}/masked.csv", "w", newline=""); mw = csv.writer(mcsv)

    if a.test == "z":
        run_z(a, names, loci, X, S, W, m, plateaus, rate, jit, masked, mw, K_eff, M_eff, t0)
    else:
        run_band(a, names, loci, X, S, W, ani, jit, masked, mw, K_eff, M_eff, t0)

    mcsv.close()
    nmasked = sum(1 for _ in open(f"{a.out}/masked.csv")) - 1
    out = write_masked(a.aln, a.out, fmt, names, [mb.decode() for mb in masked], nexus)
    print(f"masked {nmasked} cells -> {out}\ndone {time.time() - t0:.0f}s", file=sys.stderr)


# ------------------------------------------------------------------ the z test
def run_z(a, names, loci, X, S, W, m, plateaus, rate, jit, masked, mw, K_eff, M_eff, t0):
    n = len(names)
    mj = np.where(np.isfinite(m), m, np.inf) + jit
    pairs = open(f"{a.out}/pairs.tsv", "w")
    pairs.write("locus\ttaxon\trank\tcomparator\tpair_div\toverlap\tidentity\texpected\tz\tratio\texcess\tfail\n")
    cells = open(f"{a.out}/cells.tsv", "w")
    cells.write("locus\ttaxon\tn_comp\tcomparator\tpair_div\tfocal_sites\toverlap\tidentity\texpected"
                "\tz\tmedian_z\tmedian_ratio\tn_fail\tflag\tambiguous\n")
    per = open(f"{a.out}/per_locus.tsv", "w")
    per.write("locus\tsites\ttaxa\ttested\tpairs\trate\tplateau\tmed_z\ts_z\tcut\tflagged\t"
              + "\t".join(f"fail{j}" for j in range(M_eff, K_eff + 1)) + "\tambiguous\tstatus\n")
    trends = open(f"{a.out}/trends.tsv", "w"); trends.write("locus\tpair_div\ttrend\tcutoff\n")
    mw.writerow(["locus", "taxon", "sites_masked", "n_comp", "n_fail", "median_z", "median_ratio",
                 "ambiguous", "comparator", "pair_div", "overlap", "identity", "expected", "z"])
    grid_m = np.logspace(np.log10(0.05), np.log10(80.0), 200)
    for li, (name, sites) in enumerate(loci):
        ov, sm = W[li], S[li]
        nres = ov.diagonal().astype(int)
        need = np.minimum(a.fit_overlap, nres)[:, None]
        score = np.where((ov >= need) & np.isfinite(mj), -mj, -np.inf)
        np.fill_diagonal(score, -np.inf)
        skip = None
        top, valid, ncomp = top_comparators(score, K_eff, nres)
        idx = np.flatnonzero(ncomp > 0)
        if len(idx) == 0:
            skip = "skipped-no-cells"
        else:
            gi, ci, rk = flatten_pairs(idx, top, valid, K_eff)
            o = ov[gi, ci]; same = sm[gi, ci]; mis = o - same
            ident = 100.0 * same / o; dob = 100.0 - ident
            mp = m[gi, ci]
            fit = o >= a.fit_overlap
            if fit.sum() < 8: skip = f"skipped-few-fit-pairs({int(fit.sum())})"
        if skip:
            per.write(f"{name}\t{len(sites)}\t{int((ov.diagonal() > 0).sum())}\t0\t0\t"
                      f"{rate[li] if np.isfinite(rate[li]) else float('nan'):.3f}\t{plateaus[li]:.1f}"
                      f"\tnan\tnan\tnan\t0\t" + "\t".join("0" for _ in range(M_eff, K_eff + 1))
                      + f"\t0\t{skip}\n")
            print(f"{name}: {skip}", file=sys.stderr); continue
        nmin_eff = int(np.clip(fit.sum() // 10, 8, a.bin_min))
        mpc = np.maximum(mp, 0.05)
        rho = (dob + 50.0 / o) / (mpc + 50.0 / o)          # locus-vs-pair divergence ratio, stabilised
        cx, cy = trend_fit(np.log10(mpc[fit]), rho[fit], a.bin, nmin_eff)
        cy = np.maximum(cy, 1e-3)
        rho_hat = np.interp(np.log10(mpc), cx, cy)          # flat outside support = constant ratio
        dhat = np.minimum(np.maximum(rho_hat * mpc, 1e-3), plateaus[li])
        with np.errstate(invalid="ignore", divide="ignore"):
            q = np.clip(dhat / 100.0, 0.5 / np.maximum(o, 1), 0.99)
            z = (mis - o * q) / np.sqrt(o * q * (1 - q))
            ratio = (dob + 50.0 / o) / (dhat + 50.0 / o)
        excess = mis - o * q
        med_z = float(np.median(z[fit])); s_z = max(1.4826 * float(np.median(np.abs(z[fit] - med_z))), 1.0)
        cut = med_z + a.k * s_z
        fail = (z > cut) & (ratio >= a.min_ratio) & (excess >= a.min_excess)
        for j in range(len(gi)):
            pairs.write(f"{name}\t{names[gi[j]]}\t{rk[j] + 1}\t{names[ci[j]]}\t{mp[j]:.3f}\t{int(o[j])}"
                        f"\t{ident[j]:.2f}\t{100 - dhat[j]:.2f}\t{z[j]:.2f}\t{ratio[j]:.2f}\t{excess[j]:.1f}\t{int(fail[j])}\n")
        starts = np.r_[0, np.cumsum(ncomp[idx])[:-1]]
        fcount = {j: 0 for j in range(M_eff, K_eff + 1)}; namb = 0; nres_i = nres
        for g, s0, kc in zip(idx, starts, ncomp[idx]):
            sl = slice(s0, s0 + kc)
            mz = float(np.median(z[sl])); mr = float(np.median(ratio[sl])); nf = int(fail[sl].sum())
            flag = int(kc >= K_eff and nf >= M_eff)
            amb = int(kc < K_eff and nf >= kc // 2 + 1)
            if flag: fcount[min(nf, K_eff)] += 1
            namb += amb
            cells.write(f"{name}\t{names[g]}\t{kc}\t{names[ci[s0]]}\t{mp[s0]:.3f}\t{nres_i[g]}\t{int(o[s0])}"
                        f"\t{ident[s0]:.2f}\t{100 - dhat[s0]:.2f}\t{z[s0]:.2f}\t{mz:.2f}\t{mr:.2f}\t{nf}\t{flag}\t{amb}\n")
            if flag or amb:
                row = masked[g]; cnt = 0
                for site in sites:
                    if row[site] != ord("-"): cnt += 1
                    row[site] = ord("-")
                mw.writerow([name, names[g], cnt, kc, nf, f"{mz:.2f}", f"{mr:.2f}", amb, names[ci[s0]],
                             f"{mp[s0]:.3f}", int(o[s0]), f"{ident[s0]:.2f}", f"{100 - dhat[s0]:.2f}", f"{z[s0]:.2f}"])
        gd = np.minimum(np.maximum(np.interp(np.log10(grid_m), cx, cy) * grid_m, 1e-3), plateaus[li])
        for gm, gdv in zip(grid_m, gd):
            trends.write(f"{name}\t{gm:.3f}\t{100 - gdv:.2f}\t{100 - min(a.min_ratio * gdv, plateaus[li]):.2f}\n")
        nflag = sum(fcount.values()); fstr = "/".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
        per.write(f"{name}\t{len(sites)}\t{int((nres_i > 0).sum())}\t{len(idx)}\t{len(gi)}\t"
                  f"{rate[li] if np.isfinite(rate[li]) else float('nan'):.3f}\t{plateaus[li]:.1f}"
                  f"\t{med_z:.2f}\t{s_z:.2f}\t{cut:.2f}\t{nflag}\t"
                  + "\t".join(str(fcount[j]) for j in range(M_eff, K_eff + 1)) + f"\t{namb}\ttested\n")
        print(f"{name}: {len(idx)} cells, {len(gi)} pairs, cut z {cut:.1f}, flagged {nflag} "
              f"(fails {M_eff}..{K_eff} of {K_eff}: {fstr}), ambiguous {namb}", file=sys.stderr)
    for fh in (pairs, cells, per, trends): fh.close()


# ------------------------------------------------------------------ the legacy band test
def run_band(a, names, loci, X, S, W, ani, jit, masked, mw, K_eff, M_eff, t0):
    n = len(names)
    ani_j = np.where(np.isnan(ani), -np.inf, ani + jit * 1e3)
    lg = lambda v: np.log10(np.maximum(100.0 - v, 0.05))
    pairs = open(f"{a.out}/pairs.tsv", "w"); pairs.write("locus\ttaxon\trank\tcomparator\tani\toverlap\tidentity\texpected\tresidual\tfail\n")
    cells = open(f"{a.out}/cells.tsv", "w")
    cells.write("locus\ttaxon\tn_comp\tcomparator\tani\tfocal_sites\toverlap\tidentity\texpected\tresidual\tmedian_residual\tn_fail\tflag\tambiguous\n")
    per = open(f"{a.out}/per_locus.tsv", "w")
    per.write("locus\tsites\ttaxa\ttested\tpairs\trate\tmed\tspread\tband\tflagged\t"
              + "\t".join(f"fail{j}" for j in range(M_eff, K_eff + 1)) + "\tambiguous\n")
    trends = open(f"{a.out}/trends.tsv", "w"); trends.write("locus\tani\ttrend\n")
    mw.writerow(["locus", "taxon", "sites_masked", "n_comp", "n_fail", "median_residual", "ambiguous",
                 "comparator", "ani", "overlap", "identity", "expected", "residual"])
    grid = 100.0 - np.logspace(np.log10(60), np.log10(0.05), 200)
    for li, (name, sites) in enumerate(loci):
        ov, sm = W[li], S[li]
        nres = ov.diagonal().astype(int)
        need = np.minimum(a.fit_overlap, nres)[:, None]
        score = np.where(ov >= need, ani_j, -np.inf); np.fill_diagonal(score, -np.inf)
        top, valid, ncomp = top_comparators(score, K_eff, nres)
        idx = np.flatnonzero(ncomp > 0)
        if len(idx) < 2 * a.bin_min:
            print(f"{name}: only {len(idx)} testable cells, skipped", file=sys.stderr); continue
        gi, ci, rk = flatten_pairs(idx, top, valid, K_eff)
        o = ov[gi, ci]; same = sm[gi, ci]
        ident = 100.0 * same / o; an = ani[gi, ci]; fit = o >= a.fit_overlap
        if fit.sum() < a.bin_min:
            print(f"{name}: only {int(fit.sum())} pairs with >= {a.fit_overlap} shared sites, skipped", file=sys.stderr); continue
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = float(np.median((100.0 - ident[fit]) / np.maximum(100.0 - an[fit], 0.1)))
        cx, cy = trend_fit(lg(an[fit]), ident[fit], a.bin, a.bin_min)
        expected = np.interp(lg(an), cx, cy); resid = ident - expected
        med = np.median(resid[fit]); spread = 1.4826 * np.median(np.abs(resid[fit] - med))
        band = max(a.min_band, a.k * spread); fail = resid < med - band
        for j in range(len(gi)):
            pairs.write(f"{name}\t{names[gi[j]]}\t{rk[j] + 1}\t{names[ci[j]]}\t{an[j]:.2f}\t{int(o[j])}\t{ident[j]:.1f}\t{expected[j]:.1f}\t{resid[j]:.1f}\t{int(fail[j])}\n")
        starts = np.r_[0, np.cumsum(ncomp[idx])[:-1]]
        fcount = {j: 0 for j in range(M_eff, K_eff + 1)}; namb = 0
        for g, s0, kc in zip(idx, starts, ncomp[idx]):
            r = resid[s0:s0 + kc]; mres = float(np.median(r)); nf = int(fail[s0:s0 + kc].sum())
            below = mres < med - band
            flag = int(below and kc >= K_eff); amb = int(below and kc < K_eff)
            if flag: fcount[min(max(nf, M_eff), K_eff)] += 1
            namb += amb
            cells.write(f"{name}\t{names[g]}\t{kc}\t{names[ci[s0]]}\t{an[s0]:.2f}\t{nres[g]}\t{int(o[s0])}\t{ident[s0]:.1f}\t{expected[s0]:.1f}\t{resid[s0]:.1f}"
                        f"\t{mres:.1f}\t{nf}\t{flag}\t{amb}\n")
            if flag or amb:
                row = masked[g]; cnt = 0
                for site in sites:
                    if row[site] != ord("-"): cnt += 1
                    row[site] = ord("-")
                mw.writerow([name, names[g], cnt, kc, nf, f"{mres:.1f}", amb, names[ci[s0]], f"{an[s0]:.2f}", int(o[s0]), f"{ident[s0]:.1f}", f"{expected[s0]:.1f}", f"{resid[s0]:.1f}"])
        for gx, gy in zip(grid, np.interp(lg(grid), cx, cy)): trends.write(f"{name}\t{gx:.2f}\t{gy:.2f}\n")
        nflag = sum(fcount.values()); fstr = "/".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
        per.write(f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t{len(idx)}\t{len(gi)}\t{rate:.3f}\t{med:.2f}\t{spread:.2f}\t{band:.1f}\t{nflag}\t"
                  + "\t".join(str(fcount[j]) for j in range(M_eff, K_eff + 1)) + f"\t{namb}\n")
        print(f"{name}: {len(idx)} cells, {len(gi)} pairs, band {band:.0f}, flagged {nflag} (fails {M_eff}..{K_eff} of {K_eff}: {fstr}), ambiguous {namb}", file=sys.stderr)
    for fh in (pairs, cells, per, trends): fh.close()


if __name__ == "__main__":
    main()
