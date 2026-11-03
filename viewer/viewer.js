/* LocusMasker report viewer. One classic script: reads the bundle (a ZIP of .npy members and JSON) from the page
   and draws everything. Nothing here computes a result: every number comes from the bundle the tool wrote. */
(() => {
"use strict";

// ------------------------------------------------------------------ small helpers
const $ = (sel, root) => (root || document).querySelector(sel);
function el(tag, attrs, ...children) {
  const e = document.createElement(tag);
  if (attrs) for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v; else if (k === "style") e.style.cssText = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v); else if (v !== null && v !== undefined) e.setAttribute(k, v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined) e.append(c.nodeType ? c : document.createTextNode(String(c)));
  return e;
}
const fmtInt = d3.format(","), fmt1 = d3.format(".1f"), fmt2 = d3.format(".2f"), fmt3 = d3.format(".3f");
const pct = (a, b) => b ? `${fmt1(100 * a / b)}%` : "–";
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const cssVar = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
function hexToRgb(h) { h = h.replace("#", ""); return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)]; }
function mixHex(a, b, t) { const A = hexToRgb(a), B = hexToRgb(b); return "#" + [0, 1, 2].map(i => Math.round(A[i] * (1 - t) + B[i] * t).toString(16).padStart(2, "0")).join(""); }
function luminance(hex) { const [r, g, b] = hexToRgb(hex).map(v => { v /= 255; return v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); }); return 0.2126 * r + 0.7152 * g + 0.0722 * b; }

// ------------------------------------------------------------------ base64, ZIP, npy
function b64ToBytes(b64) {
  b64 = b64.replace(/\s+/g, "");
  const out = new Uint8Array(Math.floor(b64.length * 3 / 4) - (b64.endsWith("==") ? 2 : b64.endsWith("=") ? 1 : 0));
  let o = 0; const CH = 1 << 20;                            // atob in chunks that are multiples of 4 characters
  for (let i = 0; i < b64.length; i += CH) { const bin = atob(b64.slice(i, i + CH)); for (let j = 0; j < bin.length; j++) out[o++] = bin.charCodeAt(j); }
  return out;
}
class Zip {
  constructor(bytes) {
    this.bytes = bytes; const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength); this.entries = new Map();
    let eocd = -1;
    for (let i = bytes.length - 22; i >= Math.max(0, bytes.length - 22 - 65536); i--) if (dv.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
    if (eocd < 0) throw new Error("not a ZIP file (no end-of-central-directory record)");
    let count = dv.getUint16(eocd + 10, true), cdOff = dv.getUint32(eocd + 16, true);
    if (count === 0xffff || cdOff === 0xffffffff) {           // ZIP64: the locator sits just before the EOCD
      const loc = eocd - 20;
      if (dv.getUint32(loc, true) !== 0x07064b50) throw new Error("ZIP64 locator missing");
      const e64 = Number(dv.getBigUint64(loc + 8, true));
      if (dv.getUint32(e64, true) !== 0x06064b50) throw new Error("ZIP64 end record missing");
      count = Number(dv.getBigUint64(e64 + 32, true)); cdOff = Number(dv.getBigUint64(e64 + 48, true));
    }
    let p = cdOff; const td = new TextDecoder();
    for (let k = 0; k < count; k++) {
      if (dv.getUint32(p, true) !== 0x02014b50) throw new Error("bad central directory");
      const method = dv.getUint16(p + 10, true); let csize = dv.getUint32(p + 20, true), usize = dv.getUint32(p + 24, true);
      const nlen = dv.getUint16(p + 28, true), xlen = dv.getUint16(p + 30, true), clen = dv.getUint16(p + 32, true); let off = dv.getUint32(p + 42, true);
      const name = td.decode(bytes.subarray(p + 46, p + 46 + nlen));
      if (csize === 0xffffffff || usize === 0xffffffff || off === 0xffffffff) {   // ZIP64 extra field
        let q = p + 46 + nlen; const qe = q + xlen;
        while (q < qe) { const id = dv.getUint16(q, true), sz = dv.getUint16(q + 2, true); if (id === 1) { let r = q + 4; if (usize === 0xffffffff) { usize = Number(dv.getBigUint64(r, true)); r += 8; } if (csize === 0xffffffff) { csize = Number(dv.getBigUint64(r, true)); r += 8; } if (off === 0xffffffff) off = Number(dv.getBigUint64(r, true)); } q += 4 + sz; }
      }
      this.entries.set(name, { method, csize, usize, off }); p += 46 + nlen + xlen + clen;
    }
  }
  has(name) { return this.entries.has(name); }
  async read(name) {
    const e = this.entries.get(name); if (!e) throw new Error("missing member " + name);
    const dv = new DataView(this.bytes.buffer, this.bytes.byteOffset, this.bytes.byteLength);
    const nlen = dv.getUint16(e.off + 26, true), xlen = dv.getUint16(e.off + 28, true); const start = e.off + 30 + nlen + xlen;
    const raw = this.bytes.subarray(start, start + e.csize);
    if (e.method === 0) return raw.slice();
    if (e.method !== 8) throw new Error("unsupported compression " + e.method);
    const ds = new DecompressionStream("deflate-raw");
    const buf = await new Response(new Blob([raw]).stream().pipeThrough(ds)).arrayBuffer();
    return new Uint8Array(buf);
  }
}
function parseNpy(u8) {
  if (u8[0] !== 0x93 || String.fromCharCode(...u8.subarray(1, 6)) !== "NUMPY") throw new Error("not an .npy member");
  const major = u8[6]; const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
  const hlen = major === 1 ? dv.getUint16(8, true) : dv.getUint32(8, true); const hstart = major === 1 ? 10 : 12;
  const header = new TextDecoder().decode(u8.subarray(hstart, hstart + hlen));
  const descr = /'descr':\s*'([^']+)'/.exec(header)[1]; const fortran = /'fortran_order':\s*(True|False)/.exec(header)[1] === "True";
  const shape = (/'shape':\s*\(([^)]*)\)/.exec(header)[1].split(",").map(s => s.trim()).filter(s => s.length).map(Number));
  if (fortran) throw new Error("fortran order not supported");
  const off = hstart + hlen; const n = shape.reduce((a, b) => a * b, 1);
  const T = { "<f4": Float32Array, "<f8": Float64Array, "<i4": Int32Array, "<i8": BigInt64Array, "<u4": Uint32Array, "<u2": Uint16Array, "<i2": Int16Array, "|u1": Uint8Array, "|i1": Int8Array, "|b1": Uint8Array, "<u8": BigUint64Array }[descr];
  if (!T) throw new Error("unsupported dtype " + descr);
  let data;
  if ((u8.byteOffset + off) % T.BYTES_PER_ELEMENT === 0) data = new T(u8.buffer, u8.byteOffset + off, n);
  else { const copy = u8.slice(off, off + n * T.BYTES_PER_ELEMENT); data = new T(copy.buffer, 0, n); }
  return { descr, shape, data };
}

// ------------------------------------------------------------------ the bundle
async function loadBundle(bytes) {
  const zip = new Zip(bytes); const td = new TextDecoder();
  const json = async name => JSON.parse(td.decode(await zip.read(name)));
  const npy = async name => parseNpy(await zip.read(name));
  const man = await json("manifest.json");
  if (man.format_version > 9) throw new Error(lmFmt(LM_STR.fatal.newer_format, { found: man.format_version, max: 9, version: man.version }));
  const taxa = await json("taxa.json"), loci = await json("loci.json");
  const D = { zip, man, format: man.format_version || 1, taxa, loci, n: taxa.names.length, nloc: loci.names.length, K: man.options.K_eff, M: man.options.M_eff, cols: {}, cache: { codes: new Map(), seq: new Map() } };
  const byteCols = ["n_comp", "n_fail", "n_similar", "side", "focal_sites"];
  if (zip.has("cells/side2.npy")) byteCols.push("side2");
  const numCols = ["med_identity", "med_expected", "excess_m", "closest_overlap", "closest_identity", "closest_expected", "best_match", "best_match_identity", "partner"];
  for (const k of byteCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  D.dense = man.carried.numeric_columns === "dense";
  for (const k of numCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  for (const k of ["lam", "lam_p", "row_factor", "div_log10p", "sim_log10p", "dir_log10p", "depth_factor", "locus_closest", "locus_closest_log10p"]) if (zip.has(`cells/${k}.npy`)) D.cols[k] = (await npy(`cells/${k}.npy`)).data;   // not in every format; D.num gives null without them
  const bsp = zip.has("cells/bs_locus.npy") || zip.has("cells/bs_log10p.npy") ? "bs" : "rc";   // older format-8 bundles wrote the broad-sample members as cells/rc_*
  for (const k of ["log10p", "rho", "stat"]) if (zip.has(`cells/${bsp}_${k}.npy`)) D.cols[`bs_${k}`] = (await npy(`cells/${bsp}_${k}.npy`)).data;   // dense runs: the broad-sample stage-1 p, rho and variance ratio of every tested cell
  D.s2 = i => D.cols.side2 ? D.cols.side2[i] : 0;                 // format 8: the broad-sample bits (1 tested, 2 past the bar, 4 close, 8 divergent in older bundles, 16 excluded as comparator, 32 close not re-reported, 64 the divergent-cell label)
  D.rcRows = null;
  if (zip.has(`cells/${bsp}_locus.npy`)) {                         // the sequences the verdict pass inspected, with their tail p's and partner
    const rl = (await npy(`cells/${bsp}_locus.npy`)).data, rt = (await npy(`cells/${bsp}_taxon.npy`)).data, rv = {};
    for (const k of ["var_log10p", "lo_log10p", "hi_log10p", "partner", "call"]) rv[k] = zip.has(`cells/${bsp}_${k}.npy`) ? (await npy(`cells/${bsp}_${k}.npy`)).data : null;   // call from format 9; hi also in old format-8 two-arm bundles
    D.rcRows = new Map();
    for (let j = 0; j < rl.length; j++) D.rcRows.set(rt[j] * D.nloc + rl[j], { var_p: rv.var_log10p ? rv.var_log10p[j] : NaN, lo: rv.lo_log10p ? rv.lo_log10p[j] : NaN, hi: rv.hi_log10p ? rv.hi_log10p[j] : NaN, pt: rv.partner ? rv.partner[j] : -1, call: rv.call ? rv.call[j] : 0 });
  }
  D.rcAt = i => D.rcRows ? D.rcRows.get(i) : undefined;
  D.comparators = (await npy("cells/comparators.npy")).data;
  // the closest comparator is the rank-1 comparator, so it is read from the comparators: up to 0.11.0 the tool wrote the
  // taxa's five closest relatives (taxa x 5, as in taxa/closest.npy) over cells/closest.npy
  { const cl = new Int32Array(D.comparators.length / D.K); for (let r = 0; r < cl.length; r++) cl[r] = D.comparators[r * D.K]; D.cols.closest = cl; }
  D.closest = (await npy("taxa/closest.npy")).data;
  if (!D.dense) {
    const sl = (await npy("cells/sparse_locus.npy")).data, st = (await npy("cells/sparse_taxon.npy")).data;
    D.sparse = new Map(); for (let i = 0; i < sl.length; i++) D.sparse.set(st[i] * D.nloc + sl[i], i);
  }
  D.num = (k, g, li) => {                                     // a numeric column of cell (g, li), or null when not carried
    const arr = D.cols[k]; if (!arr) return null;
    if (D.dense) { const v = arr[g * D.nloc + li]; return Number.isNaN(v) ? null : v; }
    const r = D.sparse.get(g * D.nloc + li); if (r === undefined) return null; const v = arr[r]; return Number.isNaN(v) ? null : v;
  };
  if (zip.has("cells/comp_div_log10p.npy")) D.tests = { div: (await npy("cells/comp_div_log10p.npy")).data, sim: (await npy("cells/comp_sim_log10p.npy")).data };   // from format 7: both tests of every comparator, by rank
  D.compTests = (g, li) => {                                  // comparator -> {div, sim}, log10 p of its two tests against the cell; null when not carried
    if (!D.tests) return null;
    let base; if (D.dense) base = (g * D.nloc + li) * D.K; else { const q = D.sparse.get(g * D.nloc + li); if (q === undefined) return null; base = q * D.K; }
    const out = new Map(); for (let k = 0; k < D.K; k++) { const c = D.comparators[base + k]; if (c >= 0) out.set(c, { div: D.tests.div[base + k], sim: D.tests.sim[base + k] }); }
    return out;
  };
  D.comps = (g, li) => {                                      // the cell's comparators in rank order
    let base; if (D.dense) base = (g * D.nloc + li) * D.K; else { const r = D.sparse.get(g * D.nloc + li); if (r === undefined) return []; base = r * D.K; }
    const out = []; for (let r = 0; r < D.K; r++) { const c = D.comparators[base + r]; if (c >= 0) out.push(c); } return out;
  };
  D.codes = async li => {
    if (D.cache.codes.has(li)) return D.cache.codes.get(li);
    const p = (async () => { if (!zip.has(`codes/${li}.npy`)) return null; const b = await npy(`codes/${li}.npy`), c = await npy(`codes/${li}_cells.npy`); return { cells: c.data, blk: b.data, len: b.shape[1], rows: b.shape[0] }; })();
    D.cache.codes.set(li, p); return p;
  };
  D.seq = async li => {
    if (D.cache.seq.has(li)) return D.cache.seq.get(li);
    const p = (async () => { if (!zip.has(`seq/${li}.npy`)) return null; const b = await npy(`seq/${li}.npy`), t = await npy(`seq/${li}_taxa.npy`); return { taxa: t.data, blk: b.data, len: b.shape[1], rows: b.shape[0] }; })();
    D.cache.seq.set(li, p); return p;
  };
  D.similar = async li => {                                   // the locus's too-similar tests (from format 4; a locus with none has no file), or null before
    if (!D.cache.similar) D.cache.similar = new Map();
    if (D.cache.similar.has(li)) return D.cache.similar.get(li);
    const p = (async () => zip.has(`similar/${li}.json`) ? json(`similar/${li}.json`) : D.format >= 4 ? [] : null)(); D.cache.similar.set(li, p); return p;
  };
  D.trees = async li => {                                     // {cells, expected, observed, moved, tips} or null
    if (!D.cache.trees) D.cache.trees = new Map();
    if (D.cache.trees.has(li)) return D.cache.trees.get(li);
    const p = (async () => zip.has(`trees/${li}.json`) ? json(`trees/${li}.json`) : null)(); D.cache.trees.set(li, p); return p;
  };
  deriveStates(D);
  return D;
}

// ------------------------------------------------------------------ states, palette, severity
// state code per cell: 0 absent, 1 untested, 2 + f pass with f fails (f < M), 20 + f masked divergent (f clamped to M..K),
// 40 + s masked similar with s pairs (1..3; from format 4 a copy, the identified recipient of a too-similar pair, with s
// significant too-similar partners of any role); a cell masked by both tests draws divergent and the hover says both
const ST = { ABSENT: 0, UNTESTED: 1, GREY: 2, RED: 20, BLUE: 40, BOTH: 60 };   // BOTH + pairs (1..3): masked by both tests
function deriveStates(D) {
  const { n, nloc, K, M } = D; const N = n * nloc; const st = new Uint8Array(N), sev = new Uint8Array(N);
  const fs = D.cols.focal_sites, nc = D.cols.n_comp, nf = D.cols.n_fail, ns = D.cols.n_similar, sd = D.cols.side, sd2 = D.cols.side2; const both = new Uint8Array(N);
  for (let i = 0; i < N; i++) {
    if (fs[i] === 0) { st[i] = ST.ABSENT; sev[i] = 0; continue; }
    both[i] = (sd[i] & 3) === 3 ? 1 : 0;
    if (nc[i] === 0) { st[i] = ST.UNTESTED; sev[i] = 1; continue; }
    const s = sd[i];
    const s2v = sd2 ? sd2[i] : 0;
    if ((s & 3) === 3 || ((s & 1) && (s2v & 64))) { const k = clamp(ns[i], 1, 3); st[i] = ST.BOTH + k; sev[i] = 60 + k; }   // both kinds, or the divergent cell the broad sample also calls too close
    else if (s & 1 || s & 16) { const f = clamp(nf[i], M, K); st[i] = ST.RED + f; sev[i] = 20 + f; }
    else if (s & 2) { const k = clamp(ns[i], 1, 3); st[i] = ST.BLUE + k; sev[i] = 40 + k; }
    else { st[i] = ST.GREY; sev[i] = 2; }
  }
  D.state = st; D.sev = sev; D.both = both;
  // one colour per outcome, no levels (levels would mean different things as K changes; the counts live in the labels)
  const col = { pass: cssVar("--pass"), passCell: cssVar("--pass-cell"), red: cssVar("--red"), blue: cssVar("--blue"), both: cssVar("--yellow"), cell: cssVar("--cell") };
  col.blueTint = mixHex("#ffffff", col.blue, 0.45);         // the alignment: where a too-similar sequence agrees with the focal at a red site
  col.redTint = mixHex("#ffffff", col.red, 0.45);           // the alignment: where a comparator differs from the focal
  D.col = col;
  D.colorOf = code => code === ST.ABSENT ? null : code === ST.UNTESTED ? "hatch" : code < ST.RED ? col.passCell : code < ST.BLUE ? col.red : code < ST.BOTH ? col.blue : col.both;
  D.pointColorOf = code => code >= ST.GREY && code < ST.RED ? col.pass : D.colorOf(code);   // a not-masked point is a shade darker than a not-masked cell
  D.famOf = code => code >= ST.BOTH ? ST.BOTH : code >= ST.BLUE ? ST.BLUE : code >= ST.RED ? ST.RED : code;
  D.famLabel = code => { const f = D.famOf(code); return f === ST.BOTH ? "masked by both tests" : f === ST.BLUE ? (D.format >= 4 ? "masked, copied from another taxon" : "masked, too similar") : f === ST.RED ? "masked, too divergent" : D.labelOf(code); };
  // the alignment's focal-cell grades by k, the comparators differing from the focal residue, on end points fixed whatever K
  // (only the number of steps changes): the plain cell at k = 0, grey up to 24% at M - 1, red from 38% at a bare majority
  // (M) to 55% at all K, each mixed over white. Chosen with the dataviz validator: the darkest grey and the palest red are
  // ΔE 8.4 apart for colour-blind readers, and black letters read at 7.3:1 on the darkest red.
  const gi = cssVar("--grade-grey"), a0 = (255 - hexToRgb(col.cell)[0]) / (255 - hexToRgb(gi)[0]);
  D.grades = d3.range(K + 1).map(k => k <= 0 ? col.cell : k < M ? mixHex("#ffffff", gi, a0 + (0.24 - a0) * k / Math.max(M - 1, 1)) : mixHex("#ffffff", col.red, K > M ? 0.38 + 0.17 * (k - M) / (K - M) : 0.38));
  const simLabel = k => D.format >= 4 ? `masked, copied from another taxon · ${k}${k === 3 ? "+" : ""} too-similar partner${k > 1 ? "s" : ""}` : `masked, too similar · ${k}${k === 3 ? "+" : ""} pair${k > 1 ? "s" : ""}`;
  D.labelOf = code => code === ST.ABSENT ? "absent" : code === ST.UNTESTED ? "present, not tested" : code < ST.RED ? "not masked" : code < ST.BLUE ? `masked, too divergent · fails ${code - ST.RED} of ${K}` : code < ST.BOTH ? simLabel(code - ST.BLUE) : `masked by both tests: too divergent from its relatives and copied from another taxon · ${code - ST.BOTH}${code - ST.BOTH === 3 ? "+" : ""} too-similar pair${code - ST.BOTH > 1 ? "s" : ""}`;
  // one cell's label: the colour class, with the exact partner count where the class says "3+"
  D.cellLabel = (code, i) => {
    const ns = D.cols.n_similar[i], s2v = D.s2(i), sd_ = D.cols.side[i];
    if (s2v & 64 && (sd_ & 3) === 1) return "masked, too divergent · also significantly too close to part of the broad sample";
    if ((sd_ & 3) === 2 && (s2v & 4) && !ns) return "masked, too close to the broad sample";
    if ((sd_ & 3) === 1 && (s2v & 8) && nf[i] < M) return "masked, too divergent from the broad sample";
    if ((sd_ & 3) === 3 && (s2v & 8) && nf[i] < M) return "masked: copied from another taxon · also too divergent from the broad sample";
    return D.format >= 4 && code - ST.BLUE === 3 && ns >= 3 ? `masked, copied from another taxon · ${ns >= 255 ? "255+" : ns} too-similar partners` : D.labelOf(code); };
  const surf = hexToRgb(cssVar("--surface")), hatch = hexToRgb(cssVar("--surface-3"));
  const lut = new Uint8ClampedArray(64 * 4);
  for (let c = 0; c < 64; c++) { let rgb = surf; const col = D.colorOf(c); if (col === "hatch") rgb = hatch; else if (col) rgb = hexToRgb(col); lut.set([rgb[0], rgb[1], rgb[2], 255], c * 4); }
  D.lut = lut;
  D.taxCount = { masked: new Int32Array(n), div: new Int32Array(n), sim: new Int32Array(n), tested: new Int32Array(n), present: new Int32Array(n), near: new Int32Array(n) };
  D.locCount = { masked: new Int32Array(nloc), div: new Int32Array(nloc), sim: new Int32Array(nloc), tested: new Int32Array(nloc), present: new Int32Array(nloc), near: new Int32Array(nloc) };
  for (let g = 0; g < n; g++) for (let li = 0; li < nloc; li++) {
    const i = g * nloc + li; if (fs[i] === 0) continue;
    D.taxCount.present[g]++; D.locCount.present[li]++;
    if (nc[i] === 0) continue;
    D.taxCount.tested[g]++; D.locCount.tested[li]++;
    const s = sd[i]; const m = (s & 3) > 0;
    if (m) { D.taxCount.masked[g]++; D.locCount.masked[li]++; }
    if (s & 1) { D.taxCount.div[g]++; D.locCount.div[li]++; }
    if (s & 2) { D.taxCount.sim[g]++; D.locCount.sim[li]++; }
    if (!m && nf[i] > 0) { D.taxCount.near[g]++; D.locCount.near[li]++; }
  }
  // the residual axis shared by every ribbon: the range of excess_m over the carried cells
  let lo = Infinity, hi = -Infinity; const ex = D.cols.excess_m;
  for (let i = 0; i < ex.length; i++) { const v = ex[i]; if (Number.isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } }
  D.exDomain = [Math.min(Number.isFinite(lo) ? lo * 1.3 : -30, -30), Math.max(Number.isFinite(hi) ? hi * 1.3 : 300, 300)];
}
function stateReason(D, g, li) {
  const i = g * D.nloc + li;
  if (D.cols.focal_sites[i] === 0) return "the taxon has no residue at this locus";
  if (D.cols.n_comp[i] === 0) {
    if (D.taxa.status[g] !== "tested") return "taxon " + D.taxa.status[g];
    if (D.loci.status[li] !== "tested") return "locus " + D.loci.status[li];
    return "no comparator shares enough sites with the taxon here";
  }
  return null;
}
// a cell's reported too-similar partner and the cell's role in that pair, in words. From format 4 every pair of taxa at a
// locus is tested: the identified recipient (the copy) is masked, its source (the donor) is not, and a pair with no
// identified direction masks nothing; the partner reported is the source when the cell is a recipient, else its strongest
// partner. Before format 4 only the comparators were tested. Null when the cell has no partner.
const bsMan = D => D.man.broad_sample || D.man.random_comparators;   // older format-8 manifests use the old key
// a masked cell's one significance: the most significant BAR-HELD p over the tests that masked it (the close-relative
// decisive p, the fired closest-match p, the broad-sample stage-1 variance p) - a rough level of belief, not a combined p
function maskP(D, g, li) {
  const i = g * D.nloc + li, sd = D.cols.side[i]; if (!(sd & 3)) return null;
  const s2v = D.s2 ? D.s2(i) : 0, rcr = D.rcAt ? D.rcAt(i) : null, c = [];
  const dp = D.num("div_log10p", g, li), sp = D.num("sim_log10p", g, li);
  if ((sd & 1) && dp !== null && !Number.isNaN(dp)) c.push([dp, "close-relative"]);
  if ((sd & 2) && sp !== null && !Number.isNaN(sp)) c.push([sp, "closest-match"]);
  if ((s2v & (4 | 8)) && rcr && !Number.isNaN(rcr.var_p)) c.push([rcr.var_p, "broad-sample"]);
  if (!c.length) return null;
  c.sort((a, b) => a[0] - b[0]); return { p: c[0][0], test: c[0][1] };
}
function maskIndex(D) {                                     // every masked cell once, with its one significance (maskP)
  if (D._masks) return D._masks;
  const out = [];
  for (let g = 0; g < D.n; g++) for (let li = 0; li < D.nloc; li++) {
    const i = g * D.nloc + li; if (!(D.cols.side[i] & 3)) continue;
    const mp = maskP(D, g, li) || { p: NaN, test: "" };
    out.push({ g, li, i, p: mp.p, test: mp.test, code: D.state ? D.state[i] : 0 });
  }
  return (D._masks = out);
}
const starsOf = (D, p) => { const bar = D.man.test && D.man.test.log10_bar != null ? D.man.test.log10_bar : null; return p == null || Number.isNaN(p) || bar === null ? "" : p < bar - Math.log10(50) ? "***" : p < bar - Math.log10(5) ? "**" : p < bar ? "*" : ""; };

function similarRole(D, g, li) {
  const i_ = g * D.nloc + li, s2v = D.s2 ? D.s2(i_) : 0, ns_ = D.cols.n_similar[i_];
  let pt = D.num("partner", g, li);
  if ((pt === null || pt < 0) && (s2v & 64)) { const rcr = D.rcAt && D.rcAt(i_); pt = rcr && rcr.pt >= 0 ? rcr.pt : null; }   // the labelled divergent cell: its donor-side comparator
  if (pt === null || pt < 0) return null;
  const sd = D.cols.side[g * D.nloc + li], P = D.taxa.names[pt];
  if ((s2v & (4 | 8 | 64)) && !ns_) {                               // masked or labelled by the broad-sample test: the partner is a broad comparator, not a pair member
    if (s2v & 64) return { pt, role: "none", tag: "close broad comparator", fired: true, rc: true, text: `${P}: a far broad-sample comparator this divergent sequence is significantly too CLOSE to — the donor side; the divergent mask stands, this names where it points`, note: "the most extreme too-close broad-sample comparator (the donor side)" };
    if ((s2v & 8) && (sd & 1)) return { pt, role: "none", tag: "far broad comparator", fired: true, rc: true, text: `${P}: significantly too DIVERGENT from this sequence for their genome-wide divergence — contamination or a chimera confined to sites its close relatives do not cover (the broad-sample test)`, note: "the most extreme too-divergent broad-sample comparator" };
    if ((s2v & 4) && (sd & 2)) return { pt, role: "copy", tag: "close broad comparator", fired: true, rc: true, text: `${P}: significantly too close to this sequence for their genome-wide divergence — a coalescence younger than the species split allows (the broad-sample test)`, note: "the most extreme too-close broad-sample comparator" };
  }
  if (D.format < 4) return { pt, role: sd & 8 ? "recipient" : "undetermined", tag: "partner", text: `${P}${sd & 8 ? " (identified as the recipient)" : " (direction not identified)"}`, note: "the partner of the too-similar pair" };
  if (D.format >= 5) {                                           // format 5: the partner of a too-similar pair, or the closest match that was tested and passed
    if (sd & 8) return { pt, role: "uncertain", tag: "partner", fired: true, text: `${P}: the two are significantly too similar, and neither fits its own relatives clearly worse than the other (within the donor threshold), so the direction is uncertain and both are masked${sd & 2 ? "" : "; similar masking was off"}`, note: "the other member of a too-similar pair whose direction is uncertain; both masked" };
    if (sd & (2 | 32)) return { pt, role: "copy", tag: "source", fired: true, text: `${P}: this sequence is closer to it than its own relatives are, and fits its relatives worse than ${P} does, so it is the copy${sd & 2 ? ", masked" : "; similar masking was off"}`, note: "the source the focal sequence was copied from" };
    if (sd & 128) return { pt, role: "source", tag: "copy", fired: true, text: `${P} is a copy of this sequence: ${P} is closer to it than its relatives are and fits its own relatives worse; this sequence is the source, not masked for it`, note: "a copy of the focal sequence" };
    return { pt, role: "none", tag: "closest match", fired: false, text: `${P}, its closest match at the locus, is not closer than its own relatives are`, note: "the closest match at the locus; not surprisingly close" };
  }
  if (sd & (2 | 32)) return { pt, role: "recipient", tag: "source", text: `${P}, its source: this sequence is the identified recipient (the copy)${sd & 2 ? ", masked" : "; similar masking was off"}`, note: "the source the focal sequence was copied from" };
  if (sd & 128) return { pt, role: "donor", tag: "copy", text: `${P}, its copy: ${P} is the identified recipient; this sequence is the source, not masked for it`, note: "the identified recipient, a copy of the focal sequence" };
  return { pt, role: "undetermined", tag: "partner", text: `${P}, direction not identified: neither sequence is masked for it`, note: "too similar to the focal sequence, direction not identified" };
}
// the tests behind a cell's alignment rows and trees: blue, every taxon significantly too similar to it (any test of the pair,
// either member's; before format 4 the significant partner alone); red, every comparator significantly too divergent from it
// (from format 7, which carries every comparator's tests: comp); the closest match at the locus with its too-similar test
// (format 7; before it, the closest match tested outside the set, format 5); stars, *, ** and *** below alpha / N, alpha / 5N
// and alpha / 50N (the bar, and a fifth and a fiftieth of it)
function cellTests(D, g, li, sim, sr) {
  const bar = D.man.test && D.man.test.log10_bar != null ? D.man.test.log10_bar : null;
  const blue = new Map(), red = new Map(), comp = D.compTests(g, li); let closest = null;
  if (sim) for (const e of sim) { if (e.fail && (e.taxon === g || e.partner === g)) { const o = e.taxon === g ? e.partner : e.taxon; if (!blue.has(o) || e.log10p < blue.get(o).p) blue.set(o, { p: e.log10p, test: e.test }); } }
  else if (sr && (sr.fired || sr.role === "recipient")) blue.set(sr.pt, { p: D.num("sim_log10p", g, li), test: null });
  if (comp && bar !== null) for (const [c, t] of comp) if (t.div < bar && !blue.has(c)) red.set(c, t.div);
  const s2v = D.s2 ? D.s2(g * D.nloc + li) : 0, rcr = D.rcAt && D.rcAt(g * D.nloc + li);
  if (rcr && rcr.pt >= 0) {                                             // the broad-sample partner: blue for a close mask or a donor-side
    if ((s2v & 4 && !(s2v & 32)) || s2v & 64) { if (!blue.has(rcr.pt)) blue.set(rcr.pt, { p: (s2v & 64) ? rcr.lo : rcr.var_p, test: "broad" }); }   // label, red for a divergent call;
    else if (s2v & 8) { if (!red.has(rcr.pt)) red.set(rcr.pt, rcr.var_p); }     // a mask's badge is the test's own p: stage 1, bar-held
  }
  const lc = D.num("locus_closest", g, li);
  if (lc !== null && lc >= 0) closest = { taxon: lc, p: D.num("locus_closest_log10p", g, li) };
  else if (sim) for (const e of sim) if ((e.test === "outside" || e.test === "closest") && e.taxon === g) { closest = { taxon: e.partner, p: e.log10p }; break; }   // "closest" in pre-format-9 bundles
  const stars = p => p === null || p === undefined || Number.isNaN(p) || bar === null ? "" : p < bar - Math.log10(50) ? "***" : p < bar - Math.log10(5) ? "**" : p < bar ? "*" : "";
  return { bar, blue, red, comp, closest, stars, pcol: D.format >= 7 && bar !== null };
}
// every copied sequence (from format 4): the recipient, the cell identified as the copy of a too-similar pair, its donor
// (its source: the member that fits its own relatives better) and the locus
function copiesOf(D) {
  if (D._copies) return D._copies; const out = [];
  if (D.format >= 4) for (let g = 0; g < D.n; g++) for (let li = 0; li < D.nloc; li++) { const sd = D.cols.side[g * D.nloc + li]; if (!(sd & (2 | 32)) || sd & 8) continue; const p = D.num("partner", g, li); if (p !== null && p >= 0) out.push({ rec: g, don: p, li }); }   // a pair of uncertain direction names no donor
  return (D._copies = out);
}
// the copies grouped by donor, recipient or pair (mode), most copies first, ties by name
function copyGroups(D, mode) {
  const nm = D.taxa.names, keyOf = c => mode === "donor" ? String(c.don) : mode === "recipient" ? String(c.rec) : `${c.don},${c.rec}`;
  const label = k => mode === "pair" ? k.split(",").map(x => nm[+x]).join(" → ") : nm[+k];
  return d3.groups(copiesOf(D), keyOf).map(([k, v]) => ({ k, label: label(k), list: v })).sort((a, b) => b.list.length - a.list.length || a.label.localeCompare(b.label));
}

// the p every test is held to, when the bundle gives it (format 4), and " (bar 4.5e-6)" built on it.
// The machine form everywhere is log10 p; lmP is the one gate where it becomes words
const bar10 = D => D.man.test && D.man.test.log10_bar != null ? lmP(D.man.test.log10_bar) : null;
const barNote = D => bar10(D) !== null ? ` (bar ${bar10(D)})` : "";
// a p-value from a float32 column (lam_p); the column underflows to 0 below about 1e-45
const log10pTxt = p => p > 0 ? `p = ${lmP(Math.log10(p))}` : "p < 1e-45";
// the direction test's p: the exact column when the bundle has it (dir_log10p), else from lam_p
const dirTxt = (D, g, li, lp) => { const d = D.num("dir_log10p", g, li); return d !== null ? `p = ${lmP(d)}` : log10pTxt(lp); };
// the run's tests and their bar: format 4 counts every test, Bonferroni over all of them; older bundles count the distinct
// comparisons and give the bar as a normal quantile k
function testSummary(T) {
  if (T.tests == null) return `${fmtInt(T.comparisons)} distinct comparisons, Bonferroni k ${T.bonferroni_k}`;
  const crT = T.close_relative_tests != null ? T.close_relative_tests : T.divergent_tests, cmrT = T.closest_match_relative_tests != null ? T.closest_match_relative_tests : T.comparator_tests, cmoT = T.closest_match_outside_tests != null ? T.closest_match_outside_tests : T.closest_tests;
  if (cmrT != null) return `${fmtInt(T.tests)} tests (${fmtInt(crT)} close-relative, ${fmtInt(cmrT)} closest-match against a relative, ${fmtInt(cmoT)} against an outside match${T.broad_sample_tests ? `, ${fmtInt(T.broad_sample_tests)} broad-sample` : ""}), one Bonferroni bar for all of them, p < ${lmP(T.log10_bar)}`;
  return `${fmtInt(T.tests)} tests (${fmtInt(T.divergent_tests)} too divergent, ${fmtInt(T.similar_tests)} too similar, ${fmtInt(T.direction_tests)} direction), ` +
    (T.log10_bar != null ? `one Bonferroni bar for all of them, p < ${lmP(T.log10_bar)}` : `Bonferroni k ${T.bonferroni_k}`);
}

// ------------------------------------------------------------------ the app
const App = { D: null, panels: {}, hover: null, sel: { rows: new Set(), cols: new Set(), cells: new Set() }, order: { taxa: null, loci: null, taxaMode: "input", lociMode: "input" }, hidePass: false, view: "summary" };

function boot() {
  const tag = document.getElementById("lm-bundle"); const b64 = tag ? tag.textContent.trim() : "";
  if (b64 && b64 !== "@@" + "BUNDLE@@") {                            // the marker the tool replaces, split so this file never carries it
    let bytes = null; try { bytes = b64ToBytes(b64); } catch (e) { /* a cut inside the base64 */ }
    if (!bytes || !zipEnds(bytes)) { showIncomplete(b64.length); return; }     // the page was read before its file was complete
    start(bytes).catch(showError); return;
  }
  window.__lmBooted = true;
  $("#app").replaceChildren(el("div", { class: "error" }, LM_STR.fatal.no_data));
}
function showError(e) { console.error(e); window.__lmBooted = true; $("#app").replaceChildren(el("div", { class: "error" }, LM_STR.fatal.unreadable + "\n\n" + (e && e.stack || e))); }
// a complete bundle ends with the ZIP's end-of-central-directory record; a page cut before its end does not
function zipEnds(u8) { const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength); for (let i = u8.length - 22; i >= Math.max(0, u8.length - 22 - 65536); i--) if (dv.getUint32(i, true) === 0x06054b50) return true; return false; }
function showIncomplete(n) {                                          // the file on disk is short or damaged: the loading screen holds and retries by itself
  if (window.__lmRetryShort) { window.__lmRetryShort(n * 3 / 4, (window.__LM_META && __LM_META.bundle_bytes) || n * 3 / 4 + 1); return; }
  location.reload();                                                   // no loader on the page (an odd embed): the retry is a plain reload
}
async function start(bytes) {
  // the loading screen stays as-is while the data decodes: same line, full bar
  await new Promise(r => setTimeout(r, 0));                            // let the note paint before the decode blocks the thread
  try { sessionStorage.removeItem("lmRetry"); } catch (e) { }
  const D = await loadBundle(bytes); App.D = D;
  window.__lmBooted = true;
  App.order.taxa = d3.range(D.n); App.order.loci = d3.range(D.nloc);
  buildShell(); readHash(); renderAll();
  window.addEventListener("resize", () => { Overview.redraw(); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); });
}
window.LM = { start, App, openCell: (g, li) => Sequence.open(g, li), openUnit: (k, i) => Detail.open(k, i), openCopies: (m, k) => Detail.openCopies(m, k), openCard: (k, a, b) => Card.open(k, a, b), setView: v => Overview.setView(v) };   // the tests' API

// ------------------------------------------------------------------ shell: topbar, the two panel rows, hover card
function buildShell() {
  const D = App.D, m = D.man; const app = $("#app"); app.replaceChildren();
  const c = m.counts;
  const about = el("details", { class: "meta" }, el("summary", null, `${m.tool} ${m.version} · ${m.date}`),
    el("div", null, `commit ${m.commit} · python ${m.python} · numpy ${m.numpy} · ${m.system}`),
    el("div", null, `input ${m.input.file} (${(m.input.bytes / 1e6).toFixed(1)} MB, sha256 ${m.input.sha256.slice(0, 12)}…)${m.input.partitions ? ", partitions " + m.input.partitions : ""}`),
    el("div", null, `options: datatype ${m.options.datatype}, n-comp ${m.options.n_comp} (K ${m.options.K_eff}, majority ${m.options.M_eff}), fit-overlap ${m.options.fit_overlap}, alpha ${m.options.alpha}${m.options.mask_similar_only ? ", --mask-similar-only" : ""}${m.options.mask_different_only ? ", --mask-different-only" : ""}`),
    el("div", null, `test: ${testSummary(m.test)}, dataset dispersion ${m.test.phi_dataset}; carried: ${m.carried.sequences} sequences, codes for ${m.carried.codes} cells, best matches for ${m.carried.best_match_cells} cells`),
    el("div", null, `command: ${m.command}`));
  const logo = document.getElementById("lm-logo");                   // the Lo Shu grid, when the page carries it
  app.append(el("div", { class: "topbar" }, logo ? logo.content.firstElementChild.cloneNode(true) : "",
    el("h1", null, m.dataset || m.input.file), el("span", { class: "sub" }, `${fmtInt(c.masked)} of ${fmtInt(c.cells)} sequences masked (${pct(c.masked, c.cells)}) · ${fmtInt(c.taxa)} taxa × ${fmtInt(c.loci)} loci`),
    el("span", { class: "spacer" }), about));
  const panels = el("div", { class: "panels", id: "panels" }); app.append(panels);
  const top = el("div", { class: "top" }), split = el("div", { class: "hsplit", style: "display:none" }), bottom = el("div", { class: "bottom", style: "display:none" });
  panels.append(top, split, bottom); App.panelsEl = panels; App.topEl = top; App.bottomEl = bottom; App.splitEl = split;
  let drag = null;
  split.addEventListener("pointerdown", ev => { split.setPointerCapture(ev.pointerId); drag = { y: ev.clientY, h: bottom.getBoundingClientRect().height }; });
  split.addEventListener("pointermove", ev => { if (!drag) return; const h = clamp(drag.h - (ev.clientY - drag.y), 120, panels.getBoundingClientRect().height - 160); bottom.style.flexBasis = h + "px"; });
  split.addEventListener("pointerup", () => { drag = null; Overview.redraw(); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); });
  const vs = el("div", { class: "vsplit", style: "display:none" }); App.vsplitEl = vs; let vdrag = null;
  vs.addEventListener("pointerdown", ev => { vs.setPointerCapture(ev.pointerId); vdrag = { x: ev.clientX, w: App.panels.overview.el.getBoundingClientRect().width }; });
  vs.addEventListener("pointermove", ev => { if (!vdrag) return; const w = clamp(vdrag.w + (ev.clientX - vdrag.x), 260, top.getBoundingClientRect().width - 260); App.panels.overview.el.style.flex = `0 0 ${w}px`; });
  vs.addEventListener("pointerup", () => { vdrag = null; relayout(); });
  App.hover = el("div", { class: "hover" }); document.body.append(App.hover);
  App.panels.overview = Overview.create(); top.append(App.panels.overview.el);
}
function renderAll() { Overview.render(); for (const k of ["detail", "sequence"]) if (App.panels[k]) App.panels[k].render(); }
function showHover(html, x, y) { const h = App.hover; h.innerHTML = ""; h.append(html); h.style.display = "block"; const r = h.getBoundingClientRect(); const px = x + 14 + r.width > innerWidth ? x - r.width - 10 : x + 14; const py = y + 14 + r.height > innerHeight ? Math.max(4, y - r.height - 10) : y + 14; h.style.left = px + "px"; h.style.top = py + "px"; }
function hideHover() { App.hover.style.display = "none"; }
function sw(color) { const s = el("span", { class: "sw" }); if (color === "hatch") s.style.background = "var(--hatch)"; else if (color) s.style.background = color; else { s.style.border = "1px solid var(--spine)"; } return s; }
function relayout() { Overview.redraw(); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); }

// a panel frame with breadcrumb and the collapse / maximise / close tools
function makePanel(key, crumbParts, onClose) {
  const body = el("div", { class: "panel-body" });
  const crumb = el("span", { class: "crumb" });
  const setCrumb = parts => { crumb.replaceChildren(); parts.forEach((p, i) => { if (i) crumb.append(el("span", { class: "sep" }, "›")); crumb.append(typeof p === "string" ? p : p); }); };
  setCrumb(crumbParts);
  const frame = el("div", { class: "panel", "data-key": key });
  const tools = el("div", { class: "tools" });
  const collapse = el("button", { title: "collapse", onclick: () => { frame.classList.toggle("collapsed"); collapse.textContent = frame.classList.contains("collapsed") ? "▸" : "▾"; relayout(); } }, "▾");
  const maxi = el("button", { title: "maximise", onclick: () => { const on = frame.classList.toggle("maximised"); App.panelsEl.classList.toggle("has-max", on); relayout(); } }, "⤢");
  tools.append(collapse, maxi);
  if (onClose) tools.append(el("button", { title: "close", onclick: onClose }, "✕"));
  frame.append(el("div", { class: "panel-head" }, crumb, tools), body);
  return { el: frame, body, setCrumb, tools };
}
function closePanel(key) {
  const p = App.panels[key]; if (!p) return; p.el.remove(); App.panels[key] = null;
  if (key === "sequence") { App.bottomEl.style.display = "none"; App.splitEl.style.display = "none"; }
  if (key === "detail") { App.vsplitEl.style.display = "none"; App.panels.overview.el.style.flex = ""; }
  App.panelsEl.classList.remove("has-max"); writeHash(); relayout();
}

// ------------------------------------------------------------------ the overview panel: one view at a time, sized to the panel
const Overview = (() => {
  let P, matrix, glCanvas, ctx, rowCanvas, colCanvas, strip, stripCtx, wrap, viewEl = {};
  const V = { x0: 0, y0: 0, cw: 4, ch: 4, user: false };      // first display row/col, the cell size in CSS px, and whether the user has zoomed
  let dragging = null;
  const VIEWS = [["summary", "summary"], ["matrix", "matrix"], ["loci", "loci ribbons"], ["taxa", "taxon ribbons"], ["table", "table"]];
  function create() {
    P = makePanel("overview", [App.D.man.dataset || App.D.man.input.file, "overview"]);
    P.el.classList.add("overview");
    const D = App.D;
    const viewSel = el("select", { onchange: e => setView(e.target.value) }, ...VIEWS.map(([v, t]) => el("option", { value: v }, t)));
    viewEl.viewSel = viewSel;
    const taxaSel = el("select", { onchange: e => setOrder("taxa", e.target.value) }, ...[["input", "input order"], ["masked", "most masked first"], ["name", "name"], ["grouped", "grouped by shared masks"]].map(([v, t]) => el("option", { value: v }, t)));
    const lociSel = el("select", { onchange: e => setOrder("loci", e.target.value) }, ...[["input", "input order"], ["masked", "most masked first"], ["name", "name"], ["length", "length"], ["grouped", "grouped by shared masks"]].map(([v, t]) => el("option", { value: v }, t)));
    viewEl.taxaSel = taxaSel; viewEl.lociSel = lociSel;
    const search = el("input", { type: "search", placeholder: "find a taxon or locus  /", size: 22, oninput: e => doSearch(e.target.value) });
    const viewCtl = el("span", { class: "viewctl" }); viewEl.viewCtl = viewCtl;             // the open view's own controls
    const controls = el("div", { class: "controls one-row" },
      el("label", null, "view ", viewSel), viewCtl,
      el("button", { title: "download the cell table as TSV", onclick: Table.download }, "⤓ cells.tsv"),
      el("button", { title: "export the matrix as PNG", onclick: exportPng }, "⤓ PNG"),
      el("button", { class: "cardbtn", onclick: () => Card.open("dataset") }, "⎙ " + LM_STR.card.button));
    viewEl.matrixCtl = [el("label", null, "taxa ", taxaSel), el("label", null, "loci ", lociSel), search];
    P.body.append(controls);
    // the summary view: stat tiles and the worst-ten lists
    viewEl.summary = el("div", { class: "summary-view", style: "display:none" });
    viewEl.tiles = el("div", { class: "tiles" }); viewEl.worst = el("div", { class: "worst" }); viewEl.summary.append(viewEl.tiles, viewEl.worst); P.body.append(viewEl.summary);
    // the matrix
    wrap = el("div", { class: "matrix-wrap" });
    rowCanvas = el("canvas", { class: "rowlabels" }); colCanvas = el("canvas", { class: "collabels" });
    glCanvas = el("canvas", { class: "matrix-gl" }); matrix = el("canvas", { class: "matrix", tabindex: 0 });
    wrap.append(colCanvas, rowCanvas, glCanvas, matrix); P.body.append(wrap); ctx = matrix.getContext("2d");
    const zq = (fx, fy) => { const { w, h } = canvasArea(); Overview._zoomAt(w / 2, h / 2, fx, fy); };
    const mctl = el("div", { class: "mx-ctl" },
      el("button", { title: LM_STR.matrix.zoom_in, onclick: () => zq(1.5, 1.5) }, "+"),
      el("button", { title: LM_STR.matrix.zoom_out, onclick: () => zq(1 / 1.5, 1 / 1.5) }, "−"),
      el("button", { title: LM_STR.matrix.zoom_fit, onclick: () => { V.user = false; fitAll(); redraw(); writeHash(); } }, "⤢"));
    wrap.append(mctl);
    let hinted = false; try { hinted = !!sessionStorage.lmMxHint; } catch (e) { }
    if (!hinted) { const hint = el("div", { class: "mx-hint" }, LM_STR.matrix.hint); wrap.append(hint);
      Overview._dropHint = () => { hint.remove(); Overview._dropHint = null; try { sessionStorage.lmMxHint = 1; } catch (e) { } }; }
    const legend = el("div", { class: "legend" });
    legend.append(el("span", { class: "fam" }, sw(D.col.passCell), LM_STR.legend.not_masked),
      el("span", { class: "fam" }, sw(D.col.red), LM_STR.legend.divergent),
      el("span", { class: "fam" }, sw(D.col.blue), LM_STR.legend.similar),
      el("span", { class: "fam" }, sw(D.col.both), LM_STR.legend.both),
      el("span", { class: "fam" }, sw("hatch"), LM_STR.legend.untested), el("span", { class: "fam" }, sw(null), LM_STR.legend.absent));
    const help = el("button", { class: "help", title: "pinch or ctrl+wheel: zoom (shift: loci only, alt: taxa only) · two-finger scroll or drag: move · shift+drag: zoom to a box · double-click: fit · click a label: open its panel · click a cell: open the alignment · arrows move, + and − zoom, / finds, Esc closes" }, "?");
    const stripWrap = el("div", { class: "strip-wrap" }); strip = el("canvas"); stripCtx = strip.getContext("2d");
    stripWrap.append(strip, legend, help);
    P.body.append(stripWrap);
    viewEl.ribbons = el("div", { class: "ribbons-view", style: "display:none" }); P.body.append(viewEl.ribbons);
    viewEl.table = el("div", { class: "table-wrap", style: "display:none; flex: 1 1 auto;" }); P.body.append(viewEl.table);
    viewEl.matrixParts = [wrap, stripWrap];
    P.body.style.display = "flex"; P.body.style.flexDirection = "column";
    bindMatrix(); bindStrip();
    new ResizeObserver(debounce(() => { if (App.D) refresh(); }, 120)).observe(P.body);
    return { el: P.el, redraw, render };
  }
  function setView(v) {
    if (!VIEWS.some(([k]) => k === v)) v = "summary";          // an unknown view (an old link's hash) lands on the summary
    App.view = v; viewEl.viewSel.value = v;
    viewEl.viewCtl.replaceChildren(); if (v === "matrix") viewEl.viewCtl.append(...viewEl.matrixCtl);
    viewEl.matrixParts.forEach(e => e.style.display = v === "matrix" ? "" : "none");
    viewEl.summary.style.display = v === "summary" ? "" : "none";
    viewEl.ribbons.style.display = (v === "loci" || v === "taxa") ? "" : "none"; viewEl.table.style.display = v === "table" ? "" : "none";
    refresh(); writeHash();
  }
  function refresh() {                                        // draw the open view at the panel's current size
    const v = App.view;
    if (v === "matrix") { if (!V.user) fitAll(); redraw(); } else if (v === "summary") { renderTiles(); renderWorst(); }
    else if (v === "loci" || v === "taxa") Ribbons.render(viewEl.ribbons, v === "loci" ? "locus" : "taxon", viewEl.viewCtl);
    else if (v === "table") Table.render(viewEl.table);
  }
  function render() { fitAll(); setView(App.view || "summary"); }
  function renderTiles() {
    const D = App.D, c = D.man.counts, S = LM_STR.overview; const t = viewEl.tiles; t.replaceChildren();
    const taxaWith = d3.sum(D.taxCount.masked, v => v > 0), lociWith = d3.sum(D.locCount.masked, v => v > 0);
    const bs = bsMan(D) || {}; const bsN = (bs.close || 0) + (bs.divergent || 0) + (bs.both_sides || 0);
    const hero = el("div", { class: "hero" },
      el("div", { class: "hero-n" }, fmtInt(c.masked)),
      el("div", { class: "hero-k" }, el("div", { class: "hero-word" }, "masked"),
        el("div", { class: "hero-sub" }, lmFmt(S.hero_sub, { tested: fmtInt(c.cells), pct: pct(c.masked, c.cells) }))),
      el("div", { class: "hero-stats" },
        el("div", null, lmFmt(S.hero_taxa, { with: fmtInt(taxaWith), all: fmtInt(c.taxa) })),
        el("div", null, lmFmt(S.hero_loci, { with: fmtInt(lociWith), all: fmtInt(c.loci) })),
        c.excluded_taxa ? el("div", null, lmFmt(S.hero_excluded, { n: fmtInt(c.excluded_taxa) })) : null,
        c.skipped_loci ? el("div", null, lmFmt(S.hero_skipped, { n: fmtInt(c.skipped_loci) })) : null));
    const card = (color, name, n, line, extra, test) => el("div", { class: "card", onclick: () => Table.open(test) },
      el("div", { class: "card-head" }, color === "split" ? el("span", { class: "sw split" }) : sw(color), el("span", { class: "card-name" }, name)),
      el("div", { class: "card-n" }, fmtInt(n)),
      el("div", { class: "card-line" }, line),
      extra ? el("div", { class: "card-extra" }, extra) : null);
    const cards = el("div", { class: "cards" },
      card(D.col.red, S.card_close_relative, Math.max(0, c.masked_divergent - (bs.divergent || 0)), S.card_close_relative_line, null, "close_relative"),
      card(D.col.blue, S.card_closest_match, Math.max(0, c.masked_similar - (bs.close || 0)), S.card_closest_match_line, null, "closest_match"),
      ...(bsMan(D) ? [card("split", S.card_broad_sample, bsN, S.card_broad_sample_line,
        bs.divergent_labels ? lmFmt(S.card_divergent_labels, { n: fmtInt(bs.divergent_labels) }) : null, "broad_sample")] : []));
    t.append(hero, cards);
  }
  function renderWorst() {
    const D = App.D; const w = viewEl.worst; w.replaceChildren();
    const col = (title, items, open, head) => {
      const c = el("div", { class: "col" }, el("h3", head || null, title));
      const max = d3.max(items, d => d.v) || 1;
      if (new Set(items.map(d => d.v)).size < 2) {                                 // ties carry no magnitude: a chip list, no bars
        const box = el("div", { class: "chiprow" });
        for (const it of items) box.append(el("button", { class: "chip", title: it.title, onclick: () => open(it), onmouseenter: () => it.kind && highlight(it.kind, it.i), onmouseleave: () => it.kind && highlight(null) },
          sw(it.div >= it.sim ? D.col.red : D.col.blue), ` ${it.name} · ${it.label}`));
        c.append(box); return c;
      }
      for (const it of items) {
        const bar = el("div", { class: "bar" });
        const div = it.div / (it.div + it.sim || 1);
        bar.append(el("i", { style: `width:${100 * it.v / max}%; background: linear-gradient(90deg, ${D.col.red} ${100 * div}%, ${D.col.blue} ${100 * div}%)` }));
        c.append(el("div", { class: "row", title: it.title, onclick: () => open(it), onmouseenter: () => it.kind && highlight(it.kind, it.i), onmouseleave: () => it.kind && highlight(null) },
          el("span", { class: "nm" }, it.name), bar, el("span", { class: "num" }, it.label)));
      }
      return c;
    };
    const taxa = d3.range(D.n).filter(g => D.taxCount.masked[g] > 0).map(g => ({ i: g, kind: "row", name: D.taxa.names[g], v: D.taxCount.masked[g] / Math.max(D.taxCount.tested[g], 1), div: D.taxCount.div[g], sim: D.taxCount.sim[g], label: pct(D.taxCount.masked[g], D.taxCount.tested[g]), title: `${D.taxCount.masked[g]} of ${D.taxCount.tested[g]} tested loci masked` }))
      .sort((a, b) => b.v - a.v || b.div + b.sim - a.div - a.sim).slice(0, 10);
    const loci = d3.range(D.nloc).filter(l => D.locCount.masked[l] > 0).map(l => ({ i: l, kind: "col", name: D.loci.names[l], v: D.locCount.masked[l], div: D.locCount.div[l], sim: D.locCount.sim[l], label: fmtInt(D.locCount.masked[l]), title: `${D.locCount.masked[l]} of ${D.locCount.tested[l]} tested taxa masked` }))
      .sort((a, b) => b.v - a.v).slice(0, 10);
    const ct = col(LM_STR.overview.worst_taxa, taxa, it => Detail.open("taxon", it.i)), cl = col(LM_STR.overview.worst_loci, loci, it => Detail.open("locus", it.i));
    w.append(el("div", { class: "worst-row" }, ct, cl));
    // the distribution plots and the significance histogram, sized to the room left in the panel
    const room = P.body.getBoundingClientRect().height - viewEl.tiles.getBoundingClientRect().height - w.getBoundingClientRect().height - 60;
    const hh = clamp(room, 190, 420), hw = (w.getBoundingClientRect().width || 800) / 2 - 20;
    const ht = el("div", { class: "hist-wrap" }), hl = el("div", { class: "hist-wrap" });
    w.append(el("div", { class: "worst-row" }, el("div", { class: "col" }, ht), el("div", { class: "col" }, hl)));
    Hist.render(ht, "taxon", hw, hh); Hist.render(hl, "locus", hw, hh);
  }
  // ---- ordering
  function setOrder(axis, mode) {
    const D = App.D; App.order[axis + "Mode"] = mode;
    const cnt = axis === "taxa" ? D.taxCount : D.locCount, names = axis === "taxa" ? D.taxa.names : D.loci.names, N = axis === "taxa" ? D.n : D.nloc;
    let ord = d3.range(N);
    if (mode === "masked") ord.sort((a, b) => (axis === "taxa" ? cnt.masked[b] / Math.max(cnt.tested[b], 1) - cnt.masked[a] / Math.max(cnt.tested[a], 1) : cnt.masked[b] - cnt.masked[a]) || a - b);
    else if (mode === "name") ord.sort((a, b) => names[a].localeCompare(names[b]));
    else if (mode === "length") ord.sort((a, b) => D.loci.length[b] - D.loci.length[a] || a - b);
    else if (mode === "grouped") ord = groupedOrder(axis);
    App.order[axis] = ord; fitAll(); redraw(); writeHash();
  }
  function groupedOrder(axis) {
    const D = App.D; const sd = D.cols.side;
    if (axis === "taxa") {
      const lociByCount = d3.range(D.nloc).sort((a, b) => D.locCount.masked[b] - D.locCount.masked[a] || a - b); const rank = new Int32Array(D.nloc); lociByCount.forEach((l, r) => rank[l] = r);
      const key = d3.range(D.n).map(g => { let first = D.nloc; for (let li = 0; li < D.nloc; li++) if (sd[g * D.nloc + li] & 3) first = Math.min(first, rank[li]); return first; });
      return d3.range(D.n).sort((a, b) => key[a] - key[b] || D.taxCount.masked[b] - D.taxCount.masked[a] || a - b);
    }
    const taxaByCount = d3.range(D.n).sort((a, b) => D.taxCount.masked[b] - D.taxCount.masked[a] || a - b); const rank = new Int32Array(D.n); taxaByCount.forEach((g, r) => rank[g] = r);
    const key = d3.range(D.nloc).map(li => { let first = D.n; for (let g = 0; g < D.n; g++) if (sd[g * D.nloc + li] & 3) first = Math.min(first, rank[g]); return first; });
    return d3.range(D.nloc).sort((a, b) => key[a] - key[b] || D.locCount.masked[b] - D.locCount.masked[a] || a - b);
  }
  // ---- geometry
  function size() { const r = wrap.getBoundingClientRect(); return { W: Math.max(50, r.width), H: Math.max(50, r.height) }; }
  let nameW = null;                                           // the longest taxon and locus names at 12 px, measured once
  function labelSpace() {
    if (!nameW) { const c = rowCanvas.getContext("2d"); c.font = `12px ${cssVar("--sans")}`; nameW = { taxa: d3.max(App.D.taxa.names, nm => c.measureText(nm).width) || 40, loci: d3.max(App.D.loci.names, nm => c.measureText(nm).width) || 40 }; }
    const rl = V.ch >= 10 ? Math.ceil(Math.min(260, nameW.taxa * Math.min(12, V.ch - 2) / 12 + 10)) : 0, cl = V.cw >= 10 ? Math.ceil(Math.min(160, nameW.loci * Math.min(12, V.cw - 2) / 12 + 10)) : 0;
    wrap.style.setProperty("--rowlab", rl + "px"); wrap.style.setProperty("--collab", cl + "px"); return { rl, cl };
  }
  function fitAll() { const D = App.D; const { W, H } = size(); let rl = 0, cl = 0; for (let k = 0; k < 3; k++) { V.cw = Math.max(0.1, (W - rl) / D.nloc); V.ch = Math.max(0.1, (H - cl) / D.n); ({ rl, cl } = labelSpace()); } V.x0 = 0; V.y0 = 0; }   // the label gutters depend on the cell size, so iterate to a fixed point
  function canvasArea() { const { W, H } = size(); const { rl, cl } = labelSpace(); return { w: Math.max(10, W - rl), h: Math.max(10, H - cl), rl, cl }; }
  function redraw() {
    if (!App.D || App.view !== "matrix" || P.el.classList.contains("collapsed")) return;
    const D = App.D; const { w, h, rl, cl } = canvasArea(); const dpr = devicePixelRatio || 1;
    matrix.width = Math.round(w * dpr); matrix.height = Math.round(h * dpr); matrix.style.width = w + "px"; matrix.style.height = h + "px";
    V.x0 = clamp(V.x0, 0, Math.max(0, D.nloc - w / V.cw)); V.y0 = clamp(V.y0, 0, Math.max(0, D.n - h / V.ch));
    glCanvas.width = matrix.width; glCanvas.height = matrix.height; glCanvas.style.width = w + "px"; glCanvas.style.height = h + "px";
    const glOn = GLMatrix.draw(glCanvas, D, V, dpr);                 // the GPU paints the grid when in its regime; the 2d path owns the readable zoom range
    glCanvas.style.visibility = glOn ? "" : "hidden";
    if (glOn) {
      ctx.clearRect(0, 0, matrix.width, matrix.height);
      drawSelection(dpr); drawLabels(rl, cl, dpr); drawStrip();
      return;
    }
    const img = ctx.createImageData(matrix.width, matrix.height);
    const hit = new Uint8Array(matrix.width * matrix.height); const px = img.data; const lut = D.lut; const sev = D.sev, st = D.state;
    const ot = App.order.taxa, ol = App.order.loci; const nloc = D.nloc;
    const pw = V.cw * dpr, ph = V.ch * dpr;
    const colStart = new Int32Array(matrix.width), colEnd = new Int32Array(matrix.width);
    for (let x = 0; x < matrix.width; x++) { const a = Math.floor(V.x0 + x / pw), b = Math.max(a, Math.ceil(V.x0 + (x + 1) / pw) - 1); colStart[x] = a; colEnd[x] = Math.min(b, nloc - 1); }
    for (let y = 0; y < matrix.height; y++) {
      const r0 = Math.floor(V.y0 + y / ph), r1 = Math.min(Math.max(r0, Math.ceil(V.y0 + (y + 1) / ph) - 1), D.n - 1);
      if (r0 >= D.n) break;
      const rowBase = y * matrix.width * 4;
      for (let x = 0; x < matrix.width; x++) {
        const c0 = colStart[x], c1 = colEnd[x]; if (c0 >= nloc) break;
        let best = -1, bestCode = 0;
        for (let r = r0; r <= r1; r++) { const g = ot[r] * nloc; for (let c = c0; c <= c1; c++) { const i = g + ol[c]; const s = sev[i]; if (s > best) { best = s; bestCode = st[i]; } } }
        const o = rowBase + x * 4; const l = bestCode * 4; px[o] = lut[l]; px[o + 1] = lut[l + 1]; px[o + 2] = lut[l + 2]; px[o + 3] = 255;
        if (best >= 20) hit[y * matrix.width + x] = 1;
      }
    }
    if (V.cw < 1.2) for (let q = 0; q < hit.length; q++)                 // thinner than a pixel: a mark claims an unmasked right
      if (hit[q] && (q + 1) % matrix.width && !hit[q + 1]) { const o = q * 4; px[o + 4] = px[o]; px[o + 5] = px[o + 1]; px[o + 6] = px[o + 2]; px[o + 7] = 255; }   // neighbour, so single masks stay visible at any scale
    ctx.putImageData(img, 0, 0);
    if (V.cw >= 6 && V.ch >= 6) {
      ctx.save(); ctx.scale(dpr, dpr); ctx.strokeStyle = cssVar("--spine"); ctx.lineWidth = 1;
      const c0 = Math.floor(V.x0), c1 = Math.min(nloc - 1, Math.ceil(V.x0 + w / V.cw)), r0 = Math.floor(V.y0), r1 = Math.min(D.n - 1, Math.ceil(V.y0 + h / V.ch));
      for (let r = r0; r <= r1; r++) for (let c = c0; c <= c1; c++) if (st[ot[r] * nloc + ol[c]] === ST.UNTESTED) { const x = (c - V.x0) * V.cw, y = (r - V.y0) * V.ch; ctx.beginPath(); ctx.moveTo(x, y + V.ch); ctx.lineTo(x + V.cw, y); ctx.stroke(); }
      ctx.restore();
    }
    drawSelection(dpr); drawLabels(rl, cl, dpr); drawStrip();
  }
  function drawSelection(dpr) {
    const D = App.D; const { w, h } = canvasArea(); const ot = App.order.taxa, ol = App.order.loci;
    ctx.save(); ctx.scale(dpr, dpr); ctx.strokeStyle = cssVar("--focus"); ctx.lineWidth = 1.5; ctx.fillStyle = "rgba(42,111,203,0.12)";
    const rowOf = new Map(), colOf = new Map();
    if (App.sel.rows.size || App.sel.cells.size) ot.forEach((g, r) => rowOf.set(g, r));
    if (App.sel.cols.size || App.sel.cells.size) ol.forEach((l, c) => colOf.set(l, c));
    for (const g of App.sel.rows) { const r = rowOf.get(g); const y = (r - V.y0) * V.ch; if (y > -V.ch && y < h) { ctx.fillRect(0, y, w, Math.max(V.ch, 1)); ctx.strokeRect(0.5, y + 0.5, w - 1, Math.max(V.ch, 2) - 1); } }
    for (const l of App.sel.cols) { const c = colOf.get(l); const x = (c - V.x0) * V.cw; if (x > -V.cw && x < w) { ctx.fillRect(x, 0, Math.max(V.cw, 1), h); ctx.strokeRect(x + 0.5, 0.5, Math.max(V.cw, 2) - 1, h - 1); } }
    for (const key of App.sel.cells) { const g = Math.floor(key / D.nloc), l = key % D.nloc; const r = rowOf.get(g), c = colOf.get(l); const x = (c - V.x0) * V.cw, y = (r - V.y0) * V.ch; ctx.strokeRect(x - 1, y - 1, Math.max(V.cw, 2) + 2, Math.max(V.ch, 2) + 2); }
    ctx.restore();
  }
  function drawLabels(rl, cl, dpr) {
    const D = App.D; const { w, h } = canvasArea(); const ot = App.order.taxa, ol = App.order.loci;
    rowCanvas.width = Math.round(rl * dpr); rowCanvas.height = Math.round(h * dpr); rowCanvas.style.width = rl + "px"; rowCanvas.style.height = h + "px";
    colCanvas.width = Math.round(w * dpr); colCanvas.height = Math.round(cl * dpr); colCanvas.style.width = w + "px"; colCanvas.style.height = cl + "px";
    if (rl) { const c = rowCanvas.getContext("2d"); c.scale(dpr, dpr); c.font = `${Math.min(12, V.ch - 2)}px ${cssVar("--sans")}`; c.textBaseline = "middle"; c.textAlign = "right";
      const r0 = Math.floor(V.y0), r1 = Math.min(D.n - 1, Math.ceil(V.y0 + h / V.ch));
      for (let r = r0; r <= r1; r++) { const nm = D.taxa.names[ot[r]]; c.fillStyle = D.taxa.status[ot[r]] === "tested" ? cssVar("--ink") : cssVar("--faint"); c.fillText(trunc(c, nm, rl - 8), rl - 4, (r - V.y0) * V.ch + V.ch / 2); } }
    if (cl) { const c = colCanvas.getContext("2d"); c.scale(dpr, dpr); c.font = `${Math.min(12, V.cw - 2)}px ${cssVar("--sans")}`; c.textBaseline = "middle"; c.textAlign = "left";
      const c0 = Math.floor(V.x0), c1 = Math.min(D.nloc - 1, Math.ceil(V.x0 + w / V.cw));
      for (let k = c0; k <= c1; k++) { c.save(); c.translate((k - V.x0) * V.cw + V.cw / 2, cl - 4); c.rotate(-Math.PI / 2); c.fillStyle = D.loci.status[ol[k]] === "tested" ? cssVar("--ink") : cssVar("--faint"); c.fillText(trunc(c, D.loci.names[ol[k]], cl - 8), 0, 0); c.restore(); } }
  }
  function trunc(c, s, maxw) { if (c.measureText(s).width <= maxw) return s; let lo = 0, hi = s.length; while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (c.measureText(s.slice(0, mid) + "…").width <= maxw) lo = mid; else hi = mid - 1; } return s.slice(0, lo) + "…"; }
  function drawStrip() {
    const D = App.D; const dpr = devicePixelRatio || 1; const sw_ = Math.min(420, Math.max(160, size().W * 0.35)), sh = 40;
    strip.width = Math.round(sw_ * dpr); strip.height = Math.round(sh * dpr); strip.style.width = sw_ + "px"; strip.style.height = sh + "px";
    const img = stripCtx.createImageData(strip.width, strip.height); const px = img.data; const ot = App.order.taxa, ol = App.order.loci, nloc = D.nloc, sev = D.sev, st = D.state, lut = D.lut;
    for (let y = 0; y < strip.height; y++) { const r0 = Math.floor(y * D.n / strip.height), r1 = Math.max(r0, Math.floor((y + 1) * D.n / strip.height) - 1);
      for (let x = 0; x < strip.width; x++) { const c0 = Math.floor(x * nloc / strip.width), c1 = Math.max(c0, Math.floor((x + 1) * nloc / strip.width) - 1); let best = -1, code = 0;
        for (let r = r0; r <= r1; r++) { const g = ot[r] * nloc; for (let c = c0; c <= c1; c++) { const i = g + ol[c]; if (sev[i] > best) { best = sev[i]; code = st[i]; } } }
        const o = (y * strip.width + x) * 4, l = code * 4; px[o] = lut[l]; px[o + 1] = lut[l + 1]; px[o + 2] = lut[l + 2]; px[o + 3] = 255; } }
    stripCtx.putImageData(img, 0, 0);
    const { w, h } = canvasArea(); stripCtx.save(); stripCtx.scale(dpr, dpr); stripCtx.strokeStyle = cssVar("--focus"); stripCtx.lineWidth = 1.5;
    stripCtx.strokeRect(V.x0 / nloc * sw_, V.y0 / D.n * sh, Math.min(w / V.cw, nloc) / nloc * sw_, Math.min(h / V.ch, D.n) / D.n * sh); stripCtx.restore();
  }
  function exportPng() { const a = document.createElement("a"); a.download = `${App.D.man.dataset || "report"}_overview.png`; a.href = matrix.toDataURL("image/png"); a.click(); }
  // ---- interaction
  function cellAt(ev) { const r = matrix.getBoundingClientRect(); const x = ev.clientX - r.left, y = ev.clientY - r.top; return { x, y, col: Math.floor(V.x0 + x / V.cw), row: Math.floor(V.y0 + y / V.ch), cols: [Math.floor(V.x0 + x / V.cw), Math.max(Math.floor(V.x0 + x / V.cw), Math.ceil(V.x0 + (x + 1) / V.cw) - 1)], rows: [Math.floor(V.y0 + y / V.ch), Math.max(Math.floor(V.y0 + y / V.ch), Math.ceil(V.y0 + (y + 1) / V.ch) - 1)] }; }
  function bindMatrix() {
    const zoomAt = (x, y, fx, fy) => { if (Overview._dropHint) Overview._dropHint();
      if (fx !== 1) { const ncw = clamp(V.cw * fx, 0.05, 40); V.x0 = V.x0 + x / V.cw - x / ncw; V.cw = ncw; } if (fy !== 1) { const nch = clamp(V.ch * fy, 0.05, 40); V.y0 = V.y0 + y / V.ch - y / nch; V.ch = nch; } V.user = true; labelSpace(); redraw(); writeHash(); };
    Overview._zoomAt = zoomAt;
    matrix.addEventListener("wheel", ev => {
      ev.preventDefault(); const { x, y } = cellAt(ev);
      if (ev.ctrlKey || ev.metaKey) { const f = Math.exp(-ev.deltaY * 0.01); zoomAt(x, y, ev.altKey ? 1 : f, ev.shiftKey ? 1 : f); return; }   // a pinch arrives as ctrl+wheel
      V.x0 += ev.deltaX / V.cw; V.y0 += ev.deltaY / V.ch; redraw(); writeHash();                                                         // two fingers scroll
    }, { passive: false });
    const ptrs = new Map(); let pinch = null;
    matrix.addEventListener("pointerdown", ev => { matrix.setPointerCapture(ev.pointerId); ptrs.set(ev.pointerId, [ev.clientX, ev.clientY]); const c = cellAt(ev);
      if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; pinch = { d: Math.hypot(a[0] - b[0], a[1] - b[1]), cw: V.cw, ch: V.ch, x0: V.x0, y0: V.y0, mid: cellAt({ clientX: (a[0] + b[0]) / 2, clientY: (a[1] + b[1]) / 2 }) }; dragging = null; }
      else dragging = { x: c.x, y: c.y, x0: V.x0, y0: V.y0, box: ev.shiftKey, moved: false }; });
    matrix.addEventListener("pointercancel", ev => { ptrs.delete(ev.pointerId); pinch = null; dragging = null; });
    matrix.addEventListener("pointermove", ev => {
      if (ptrs.has(ev.pointerId)) ptrs.set(ev.pointerId, [ev.clientX, ev.clientY]);
      if (pinch && ptrs.size === 2) { const [a, b] = [...ptrs.values()]; const d = Math.hypot(a[0] - b[0], a[1] - b[1]); if (pinch.d > 0) { const f = d / pinch.d; V.cw = pinch.cw; V.ch = pinch.ch; V.x0 = pinch.x0; V.y0 = pinch.y0; zoomAt(pinch.mid.x, pinch.mid.y, f, f); } return; }
      const c = cellAt(ev);
      if (dragging) { const dx = c.x - dragging.x, dy = c.y - dragging.y; if (Math.abs(dx) + Math.abs(dy) > 3) dragging.moved = true;
        if (dragging.box) { redraw(); const dpr = devicePixelRatio || 1; ctx.save(); ctx.scale(dpr, dpr); ctx.strokeStyle = cssVar("--focus"); ctx.setLineDash([4, 3]); ctx.strokeRect(dragging.x, dragging.y, dx, dy); ctx.restore(); }
        else if (dragging.moved) { V.x0 = dragging.x0 - dx / V.cw; V.y0 = dragging.y0 - dy / V.ch; redraw(); }
        return; }
      hoverAt(c, ev);
    });
    matrix.addEventListener("pointerup", ev => {
      ptrs.delete(ev.pointerId); if (ptrs.size < 2) pinch = null;
      const c = cellAt(ev); const d = dragging; dragging = null; if (!d) return;
      if (d.box && d.moved) { const { w, h } = canvasArea(); const x0 = Math.min(d.x, c.x), x1 = Math.max(d.x, c.x), y0 = Math.min(d.y, c.y), y1 = Math.max(d.y, c.y);
        const c0 = V.x0 + x0 / V.cw, c1 = V.x0 + x1 / V.cw, r0 = V.y0 + y0 / V.ch, r1 = V.y0 + y1 / V.ch;
        V.cw = clamp(w / Math.max(c1 - c0, 1), 0.05, 40); V.ch = clamp(h / Math.max(r1 - r0, 1), 0.05, 40); V.x0 = c0; V.y0 = r0; V.user = true; labelSpace(); redraw(); writeHash(); return; }
      if (!d.moved) clickCell(c);
      writeHash();
    });
    matrix.addEventListener("dblclick", () => { V.user = false; fitAll(); redraw(); writeHash(); });
    matrix.addEventListener("pointerleave", () => hideHover());
    rowCanvas.addEventListener("click", ev => { const r = rowCanvas.getBoundingClientRect(); const row = Math.floor(V.y0 + (ev.clientY - r.top) / V.ch); if (row >= 0 && row < App.D.n) Detail.open("taxon", App.order.taxa[row]); });
    colCanvas.addEventListener("click", ev => { const r = colCanvas.getBoundingClientRect(); const col = Math.floor(V.x0 + (ev.clientX - r.left) / V.cw); if (col >= 0 && col < App.D.nloc) Detail.open("locus", App.order.loci[col]); });
    rowCanvas.style.cursor = "pointer"; colCanvas.style.cursor = "pointer";
    rowCanvas.addEventListener("mousemove", ev => { const r = rowCanvas.getBoundingClientRect(); const row = Math.floor(V.y0 + (ev.clientY - r.top) / V.ch); if (row >= 0 && row < App.D.n) showHover(unitHover("taxon", App.order.taxa[row]), ev.clientX, ev.clientY); });
    colCanvas.addEventListener("mousemove", ev => { const r = colCanvas.getBoundingClientRect(); const col = Math.floor(V.x0 + (ev.clientX - r.left) / V.cw); if (col >= 0 && col < App.D.nloc) showHover(unitHover("locus", App.order.loci[col]), ev.clientX, ev.clientY); });
    rowCanvas.addEventListener("mouseleave", hideHover); colCanvas.addEventListener("mouseleave", hideHover);
    document.addEventListener("keydown", ev => {
      if (ev.target.tagName === "INPUT" || ev.target.tagName === "SELECT" || ev.target.tagName === "TEXTAREA") { if (ev.key === "Escape") ev.target.blur(); return; }
      if (ev.key === "/") { ev.preventDefault(); const s = $("input[type=search]", P.el); if (s) s.focus(); return; }
      if (ev.key === "Escape") { if (App.panels.sequence) closePanel("sequence"); else if (App.panels.detail) closePanel("detail"); return; }
      if (App.view !== "matrix") return;
      const step = 40;
      if (ev.key === "ArrowLeft") V.x0 -= step / V.cw; else if (ev.key === "ArrowRight") V.x0 += step / V.cw; else if (ev.key === "ArrowUp") V.y0 -= step / V.ch; else if (ev.key === "ArrowDown") V.y0 += step / V.ch;
      else if (ev.key === "+" || ev.key === "=") { V.cw *= 1.25; V.ch *= 1.25; } else if (ev.key === "-") { V.cw /= 1.25; V.ch /= 1.25; } else return;
      V.user = true; ev.preventDefault(); labelSpace(); redraw(); writeHash();
    });
  }
  function bindStrip() {
    const move = ev => { const D = App.D; const r = strip.getBoundingClientRect(); const { w, h } = canvasArea(); const fx = (ev.clientX - r.left) / r.width, fy = (ev.clientY - r.top) / r.height; V.x0 = fx * D.nloc - w / V.cw / 2; V.y0 = fy * D.n - h / V.ch / 2; redraw(); };
    let down = false; strip.addEventListener("pointerdown", ev => { down = true; strip.setPointerCapture(ev.pointerId); move(ev); }); strip.addEventListener("pointermove", ev => down && move(ev)); strip.addEventListener("pointerup", () => { down = false; writeHash(); });
  }
  function clickCell(c) {
    const D = App.D; if (c.row < 0 || c.row >= D.n || c.col < 0 || c.col >= D.nloc) return;
    if (V.cw < 4 || V.ch < 4) { const { w, h } = canvasArea(); const fx = c.x / w, fy = c.y / h; V.cw = Math.max(V.cw * 4, 6); V.ch = Math.max(V.ch * 4, 6); V.x0 = c.col - fx * w / V.cw; V.y0 = c.row - fy * h / V.ch; V.user = true; labelSpace(); redraw(); return; }
    const g = App.order.taxa[c.row], li = App.order.loci[c.col];
    if (D.cols.focal_sites[g * D.nloc + li] === 0) return;
    Sequence.open(g, li);
  }
  function hoverAt(c, ev) {
    const D = App.D; if (c.row < 0 || c.row >= D.n || c.col < 0 || c.col >= D.nloc) { hideHover(); return; }
    const cells = []; for (let r = c.rows[0]; r <= Math.min(c.rows[1], D.n - 1); r++) for (let k = c.cols[0]; k <= Math.min(c.cols[1], D.nloc - 1); k++) cells.push([App.order.taxa[r], App.order.loci[k]]);
    cells.sort((a, b) => D.sev[b[0] * D.nloc + b[1]] - D.sev[a[0] * D.nloc + a[1]]);
    const box = el("div"); const show = cells.slice(0, 4);
    if (cells.length > 1) box.append(el("div", { class: "more" }, `${fmtInt(cells.length)} cells under the pointer, worst first:`));
    for (const [g, li] of show) box.append(cellHover(g, li, cells.length > 1));
    if (cells.length > show.length) box.append(el("div", { class: "more" }, `… and ${fmtInt(cells.length - show.length)} more; zoom in to see them`));
    showHover(box, ev.clientX, ev.clientY);
  }
  function highlight(kind, i) { App.sel.rows.clear(); App.sel.cols.clear(); if (kind === "row") App.sel.rows.add(i); if (kind === "col") App.sel.cols.add(i); redraw(); }
  function doSearch(q) {
    const D = App.D; App.sel.rows.clear(); App.sel.cols.clear(); q = q.trim(); if (!q) { redraw(); return; }
    let first = null; const lq = q.toLowerCase();
    D.taxa.names.forEach((nm, g) => { if (nm.toLowerCase().includes(lq)) { App.sel.rows.add(g); if (!first) first = ["row", g]; } });
    D.loci.names.forEach((nm, l) => { if (nm.toLowerCase().includes(lq)) { App.sel.cols.add(l); if (!first) first = ["col", l]; } });
    if (!first && typeof Fuse !== "undefined" && q.length >= 3) {                 // nothing contains it: the fuzzy pass (typos, partial ids)
      if (!D._fuse) D._fuse = new Fuse(
        D.taxa.names.map((nm, g) => ({ nm, kind: "row", i: g })).concat(D.loci.names.map((nm, l) => ({ nm, kind: "col", i: l }))),
        { keys: ["nm"], threshold: 0.35, ignoreLocation: true });
      for (const r of D._fuse.search(q, { limit: 12 })) {
        (r.item.kind === "row" ? App.sel.rows : App.sel.cols).add(r.item.i);
        if (!first) first = [r.item.kind, r.item.i];
      }
    }
    if (first) { const { w, h } = canvasArea(); if (first[0] === "row") { const r = App.order.taxa.indexOf(first[1]); V.y0 = r - h / V.ch / 2; } else { const c = App.order.loci.indexOf(first[1]); V.x0 = c - w / V.cw / 2; } }
    redraw();
  }
  return { create, render, redraw, setOrder, setView, refresh, V, highlight, labelSpace, fitAll, viewEl: () => viewEl };
})();
function debounce(fn, ms) { let t = null; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }

// ------------------------------------------------------------------ the matrix's GL fast path: the whole grid as one
// R8 state texture (taxa x loci, the display order baked in), colours from a 64x1 LUT texture, pan/zoom as uniforms.
// Minification picks CPU-built mip levels that carry the MAX-SEVERITY state of each 2x2, so a single masked cell
// survives to any zoom level. The 2d canvas above stays for selection, hatching and the zoomed-in detail
const GLMatrix = (() => {
  let gl = null, prog = null, quad = null, stateTex = null, lutTex = null, u = {}, cache = {};
  const VS = "attribute vec2 p; void main() { gl_Position = vec4(p, 0.0, 1.0); }";
  const FS = `precision highp float;
    uniform sampler2D state, lut; uniform vec4 view; uniform vec2 dims; uniform vec2 cv; uniform vec3 bg;
    void main() {
      vec2 px = vec2(gl_FragCoord.x, cv.y - gl_FragCoord.y);
      vec2 cell = view.xy + px / view.zw;
      if (cell.x < 0.0 || cell.y < 0.0 || cell.x >= dims.x || cell.y >= dims.y) { gl_FragColor = vec4(bg, 1.0); return; }
      float st = texture2D(state, cell / dims).r;
      gl_FragColor = texture2D(lut, vec2((st * 255.0 + 0.5) / 64.0, 0.5));
    }`;
  function init(canvas) {
    gl = canvas.getContext("webgl2", { antialias: false, preserveDrawingBuffer: true });
    if (!gl) return false;
    const sh = (t, src) => { const h = gl.createShader(t); gl.shaderSource(h, src); gl.compileShader(h); if (!gl.getShaderParameter(h, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(h)); return h; };
    prog = gl.createProgram(); gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS)); gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    quad = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, quad); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    for (const k of ["state", "lut", "view", "dims", "cv", "bg"]) u[k] = gl.getUniformLocation(prog, k);
    stateTex = gl.createTexture(); lutTex = gl.createTexture();
    return true;
  }
  function upload(D, ot, ol) {
    const W0 = D.nloc, H0 = D.n; const maxT = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    const kx = Math.max(1, Math.ceil(W0 / maxT)), ky = Math.max(1, Math.ceil(H0 / maxT));   // an axis past the texture limit
    const W = Math.ceil(W0 / kx), H = Math.ceil(H0 / ky);                                    // is pre-binned by max severity;
    cache.kx = kx; cache.ky = ky;                                                            // the CPU path owns the zoomed-in range
    let st = new Uint8Array(W * H), sv = new Uint8Array(W * H);
    if (sv.fill) sv.fill(0);
    for (let r = 0; r < H; r++) for (let c = 0; c < W; c++) {
      let b = -1, bc = 0;
      for (let dy = 0; dy < ky; dy++) for (let dx = 0; dx < kx; dx++) {
        const rr = r * ky + dy, cc = c * kx + dx; if (rr >= H0 || cc >= W0) continue;
        const i = ot[rr] * W0 + ol[cc]; if (D.sev[i] > b) { b = D.sev[i]; bc = D.state[i]; } }
      st[r * W + c] = bc; sv[r * W + c] = b < 0 ? 0 : b;
    }
    gl.bindTexture(gl.TEXTURE_2D, stateTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    let w = W, h = H, lvl = 0;
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, w, h, 0, gl.RED, gl.UNSIGNED_BYTE, st);
    while (w > 1 || h > 1) {                                           // the max-severity pyramid: a lone mask survives every zoom level
      const w2 = Math.max(1, w >> 1), h2 = Math.max(1, h >> 1);
      const st2 = new Uint8Array(w2 * h2), sv2 = new Uint8Array(w2 * h2);
      for (let y = 0; y < h2; y++) for (let x = 0; x < w2; x++) {
        let b = -1, bc = 0;
        for (let dy = 0; dy < 2; dy++) for (let dx = 0; dx < 2; dx++) {
          const sy = Math.min(h - 1, y * 2 + dy), sx = Math.min(w - 1, x * 2 + dx); const q = sy * w + sx;
          if (sv[q] > b) { b = sv[q]; bc = st[q]; } }
        st2[y * w2 + x] = bc; sv2[y * w2 + x] = b < 0 ? 0 : b;
      }
      lvl++; gl.texImage2D(gl.TEXTURE_2D, lvl, gl.R8, w2, h2, 0, gl.RED, gl.UNSIGNED_BYTE, st2);
      st = st2; sv = sv2; w = w2; h = h2;
    }
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST_MIPMAP_NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.bindTexture(gl.TEXTURE_2D, lutTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 64, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(D.lut.buffer.slice(0)));
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE); gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return true;
  }
  function draw(canvas, D, V, dpr) {
    try {
      if (!gl && !init(canvas)) return false;
      if (gl.isContextLost && gl.isContextLost()) { gl = null; if (!init(canvas)) return false; }
      const ot = App.order.taxa, ol = App.order.loci;
      if (cache.D !== D || cache.ot !== ot || cache.ol !== ol) {
        if (!upload(D, ot, ol)) return false; cache.D = D; cache.ot = ot; cache.ol = ol;
      }
      const kx = cache.kx || 1, ky = cache.ky || 1;
      const fidX = kx === 1 || V.cw * dpr * kx <= 1.2, fidY = ky === 1 || V.ch * dpr * ky <= 1.2;   // a pre-binned axis is exact only below a pixel per bin
      const inRegime = !(V.cw >= 6 && V.ch >= 6) && fidX && fidY;      // fully readable cells go to the 2d path for gridlines and hatching
      if (!inRegime) return false;
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.useProgram(prog);
      gl.bindBuffer(gl.ARRAY_BUFFER, quad); const loc = gl.getAttribLocation(prog, "p"); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, stateTex); gl.uniform1i(u.state, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, lutTex); gl.uniform1i(u.lut, 1);
      gl.uniform4f(u.view, V.x0 / kx, V.y0 / ky, V.cw * dpr * kx, V.ch * dpr * ky);
      gl.uniform2f(u.dims, Math.ceil(D.nloc / kx), Math.ceil(D.n / ky));
      gl.uniform2f(u.cv, canvas.width, canvas.height);
      const bgc = hexToRgb(cssVar("--surface"));
      gl.uniform3f(u.bg, bgc[0] / 255, bgc[1] / 255, bgc[2] / 255);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      return true;
    } catch (e) { console.warn("GL matrix fell back to canvas:", e); gl = null; return false; }
  }
  return { draw };
})();

// hover content for a cell, shared by the matrix, the ribbons and the table
function cellHover(g, li, compact) {
  const D = App.D; const i = g * D.nloc + li; const code = D.state[i]; const col = D.colorOf(code);
  const t = el("table"); const row = (k, v) => { if (v === null || v === undefined || v === "") return; t.append(el("tr", null, el("td", { class: "k" }, k), el("td", null, v))); };
  const head = el("div", { class: "state" }, sw(col), `${D.taxa.names[g]} · ${D.loci.names[li]}`);
  const box = el("div", null, head);
  const reason = stateReason(D, g, li);
  const sd = D.cols.side[i], nf = D.cols.n_fail[i], ns = D.cols.n_similar[i], nc = D.cols.n_comp[i], v4 = D.format >= 4;
  let label = D.cellLabel(code, i); if ((sd & 3) === 3 && D.format < 5) label += v4 ? " · also copied from another taxon" : " · also masked as too similar"; if (sd & 4) label += " · fewer comparators than requested"; if (sd & 16) label += " (flagged; divergent masking was off)"; if (sd & 32) label += v4 ? " · identified as a copy (similar masking was off)" : " · a similar-side fail (similar masking was off)";
  row("state", reason ? `${D.cellLabel(code, i)}: ${reason}` : label);
  const mp_ = maskP(D, g, li);
  if (mp_) row("mask significance", `p = ${lmP(mp_.p)}${starsOf(D, mp_.p) ? " " + starsOf(D, mp_.p) : ""} (the ${mp_.test} test)`);
  if (nc > 0) {
    row("residues", fmtInt(D.cols.focal_sites[i]));
    // ---- the closest-match test (one direction: too similar; two arms by where the expectation comes from)
    if (D.format >= 5 && ns) row("too-similar pairs", `${ns}${ns >= 255 ? "+" : ""}${ns > 1 ? "; the strongest below" : ""}`);
    else if (v4 && ns) row("too-similar partners", `${ns}${ns >= 255 ? "+" : ""} significant, among every taxon at the locus${ns > 1 ? (sd & (2 | 32) ? "; its source below" : "; the strongest below") : ""}`);
    const sr = similarRole(D, g, li);
    if (sr && D.format >= 5) { const sp = D.num("sim_log10p", g, li); row(sr.fired ? (sr.role === "copy" ? "copied from" : sr.role === "uncertain" ? "too similar to" : sr.role === "none" ? (sr.tag === "far broad comparator" ? "too divergent from" : "too close to") : "copied by") : "closest match", sr.text); if (sp !== null && !sr.rc) row("closest-match test", `p = ${lmP(sp)}${barNote(D)}`);
      if (sr.fired && !sr.rc) { const xp = D.num("excess_m", sr.pt, li), x0 = D.num("excess_m", g, li), thr = D.man.options.donor_threshold; if (xp !== null && x0 !== null) row("which is the copy", `${sr.role === "uncertain" ? "uncertain" : "the member that fits its own relatives worse"}, by the close-relative excess: ${fmt2(Math.abs(x0))} m here against ${fmt2(Math.abs(xp))} m for ${D.taxa.names[sr.pt]}${thr ? ` (the larger must be more than ${thr} times the smaller to call it)` : ""}`); } }
    else if (sr) { const lam = D.num("lam", g, li), lp = D.num("lam_p", g, li), sp = D.num("sim_log10p", g, li); row("too similar to", sr.text); if (sp !== null) row("closest-match test", `p = ${lmP(sp)}${barNote(D)}`); if (lam !== null) row("direction LRT", v4 && lp !== null ? `|2 ln LR| = ${fmt1(lam)}, ${dirTxt(D, g, li, lp)} (approximate χ², 1 df${bar10(D) !== null ? `; bar ${bar10(D)}` : ""})` : `|2 ln LR| = ${fmt1(lam)}, p = ${lp < 1e-4 ? lp.toExponential(1) : fmt3(lp)} (χ², 1 df)`); const rf = D.num("row_factor", g, li); if (rf !== null && rf < 1) row("row factor", `${fmt3(rf)} (the pair was tested against ${fmt3(rf)} × its expectation: ${v4 ? "one of the two" : "the sequence"} sits below expectation against its other comparators too)`); const gf = D.num("depth_factor", g, li); if (gf !== null && gf < 1) row("depth factor", `${fmt3(gf)} (pairs this far apart sit at ${fmt3(gf)} × expectation at this locus; the pair was tested against row factor × depth factor × its expectation)`); }
    const bm = D.num("best_match", g, li);
    if (bm !== null && bm >= 0 && !compact) { const close = Array.from(D.closest.subarray(g * 5, g * 5 + 5)).includes(bm); row("best match", `${D.taxa.names[bm]} (${fmt1(D.num("best_match_identity", g, li))}% identity over shared sites${close ? ", one of the 5 closest relatives" : ""})`); }
    // ---- the broad-sample test (one direction at most: the stronger side is the call)
    const s2v = D.s2 ? D.s2(i) : 0;
    if (s2v & 1) {
      const rcr = D.rcAt && D.rcAt(i); const rcman = bsMan(D);
      let vp = rcr && !Number.isNaN(rcr.var_p) ? rcr.var_p : null;
      if (vp === null && D.cols.bs_log10p && !Number.isNaN(D.cols.bs_log10p[i])) vp = D.cols.bs_log10p[i];
      if (!compact || s2v & (2 | 4 | 8 | 64)) row("broad sample", `${rcman ? rcman.per_cell : ""} comparators drawn from the taxon's pool${vp !== null ? `; stage 1 p = ${lmP(vp)}` : ""}${s2v & 2 ? ` — past the run's bar${barNote(D)}` : ", nothing unusual"}`);
      if (D.cols.bs_rho && !Number.isNaN(D.cols.bs_rho[i]) && !compact) row("own rate ρ", `${fmt3(D.cols.bs_rho[i])} × the locus rate (the median over its comparators)`);
      if (s2v & 4 && rcr) row("broad: too close", `the far shortfall, p = ${lmP(rcr.lo)} at α${s2v & 32 ? " — another mask already covers the cell" : ""}`);
      if (s2v & 8 && rcr) row("broad: too divergent", `half the far sample leans high, p = ${lmP(rcr.hi)} at α${s2v & 32 ? " — another mask already covers the cell" : ""}`);
      if (s2v & 64 && rcr) row("broad: the donor side", `the divergent mask stands; without ρ the far shortfall names the donor side, p = ${lmP(rcr.lo)}${barNote(D)}`);
      if ((s2v & 2) && !(s2v & (4 | 8 | 64)) && rcr && rcr.pt >= 0 && rcr.call) row("broad: the lean", `${rcr.call > 0 ? "toward too close" : "toward too divergent"}, ${D.taxa.names[rcr.pt]} the most extreme — not significant at α, nothing masked`);
    }
    // ---- the close-relative test (one direction: too divergent)
    row("comparators", `${nc} tested, ${nf} failed as too divergent${ns && !v4 ? `, ${ns} as too similar` : ""}`);
    const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li), xm = D.num("excess_m", g, li);
    if (mi !== null) row("divergence, median over comparators", `${fmt2(100 - mi)}% observed, ${fmt2(100 - me)}% expected`);
    const cl = D.num("closest", g, li);
    if (cl !== null && cl >= 0 && !compact) row("closest comparator", `${D.taxa.names[cl]}: ${fmt2(100 - D.num("closest_identity", g, li))}% observed vs ${fmt2(100 - D.num("closest_expected", g, li))}% expected over ${fmtInt(D.num("closest_overlap", g, li))} shared sites`);
    if (xm !== null) row("excess beyond expectation", `${fmt2(xm)} m`);
    const dp = D.num("div_log10p", g, li);
    if (dp !== null && !compact) row("close-relative test", `p = ${lmP(dp)}, the comparison that decides the majority${barNote(D)}`);
    if (s2v & 16 && !compact) row("broad-sample role", "excluded as a comparator at this locus (too divergent, or a masked copy)");
    if (!compact) row("", el("span", { class: "more" }, "click: open the alignment"));
  }
  box.append(t); return box;
}
function unitHover(kind, i) {
  const D = App.D; const c = kind === "taxon" ? D.taxCount : D.locCount; const nm = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
  const t = el("table"); const row = (k, v) => t.append(el("tr", null, el("td", { class: "k" }, k), el("td", null, v)));
  row(kind, nm); row("status", kind === "taxon" ? D.taxa.status[i] : D.loci.status[i]);
  row("sequences", `${fmtInt(c.present[i])} present, ${fmtInt(c.tested[i])} tested`);
  row("masked", `${fmtInt(c.masked[i])} (${pct(c.masked[i], c.tested[i])}): ${c.div[i]} too divergent, ${c.sim[i]} ${D.format >= 4 ? "copied from another taxon" : "too similar"}; ${c.near[i]} near misses`);
  if (kind === "locus") { row("length", `${fmtInt(D.loci.length[i])} sites`); if (D.loci.rate[i] !== null) row("rate", `${fmt3(D.loci.rate[i])} × the typical locus`); row("plateau", `${fmt1(D.loci.plateau[i])}%`); if (D.loci.phi[i] !== null) row("dispersion φ", fmt3(D.loci.phi[i])); }
  else { const cl = Array.from(D.closest.subarray(i * 5, i * 5 + 5)).filter(x => x >= 0).map(x => D.taxa.names[x]); if (cl.length) row("closest by m", cl.join(", ")); }
  return el("div", null, t, el("div", { class: "more" }, "click: open the panel"));
}

// ------------------------------------------------------------------ residual ribbons: the cells of one unit on an asinh axis
// draws the residuals of the cells of a locus (by="locus") or a taxon into a canvas row; returns the drawn points for hover
function minCompSim(D, g, li) {                              // the strongest in-set pair evidence: min over the comparator slots
  if (!D.tests) return null;
  let base; if (D.dense) base = (g * D.nloc + li) * D.K; else { const q = D.sparse.get(g * D.nloc + li); if (q === undefined) return null; base = q * D.K; }
  let m = null; for (let k = 0; k < D.K; k++) { const v = D.tests.sim[base + k]; if (Number.isFinite(v) && (m === null || v < m)) m = v; }
  return m;
}
function cmArms(D, g, li) {                                  // the closest-match test's two evidences for one cell
  const i = g * D.nloc + li;
  const out = (D.cols.side[i] & 64) ? D.num("locus_closest_log10p", g, li) : null;   // the outside arm ran
  const pair = D.num("sim_log10p", g, li);                   // the strongest pair, either member's side, either arm
  const slots = minCompSim(D, g, li);
  const pairIsOut = pair !== null && out !== null && Math.abs(pair - out) < 1e-4;    // that pair WAS the outside test
  let rel = pairIsOut ? slots : pair !== null ? (slots !== null ? Math.min(pair, slots) : pair) : slots;
  return { rel, out };
}
function residualOf(D, g, li, mode) {                        // one residual per cell per TEST ARM, on its own axis
  if (mode === "cm_relative") { const v = cmArms(D, g, li).rel; return v === null ? null : -v; }
  if (mode === "cm_outside") { const v = cmArms(D, g, li).out; return v === null ? null : -v; }
  const i = g * D.nloc + li;
  if (mode === "broad_sample") {                             // observed / expected variance, log2
    const v = D.cols.bs_stat ? D.cols.bs_stat[i] : NaN;
    return Number.isFinite(v) && v > 0 ? Math.log2(v) : null;
  }
  return D.num("excess_m", g, li);                           // close_relative: the excess beyond expectation, in m
}
function maskKind(D, g, li) {                                // which test masked this cell (the ribbons' own-test colours)
  const i = g * D.nloc + li; const sd = D.cols.side[i]; if (!(sd & 3)) return null;
  const s2v = D.s2 ? D.s2(i) : 0;
  const bsC = (s2v & 4) !== 0 && (s2v & 32) === 0, bsD = (s2v & 8) !== 0 && (s2v & 32) === 0;
  const nc = D.cols.n_comp[i], nf = D.cols.n_fail[i];
  const cr = (sd & 1) !== 0 && (nc === D.K ? nf >= D.M : nf >= (nc >> 1) + 1);
  const cm = (sd & 2) !== 0 && D.cols.n_similar[i] > 0;
  let arm = null;
  if (cm) { const a = cmArms(D, g, li);
    arm = a.out !== null && (a.rel === null || a.out <= a.rel + 1e-4) ? "out" : "rel"; }
  return { cr, cm, arm, bsC, bsD };
}
function ribbonColour(D, g, li, mode) {                      // the mark's colour on this ribbon, null for a grey point
  const mk = maskKind(D, g, li); if (!mk) return null;
  if (mode === "close_relative") return mk.cr ? D.col.red : null;
  if (mode === "cm_relative") return mk.cm && mk.arm === "rel" ? D.col.blue : null;
  if (mode === "cm_outside") return mk.cm && mk.arm === "out" ? D.col.blue : null;
  if (mode === "broad_sample") return mk.bsD ? D.col.red : mk.bsC ? D.col.blue : null;
  return null;
}
function residualScale(D, mode, W, padL, padR) {
  if (mode === "cm_relative" || mode === "cm_outside") return d3.scaleLinear().domain([0, Math.max(6, -(D.man.test && D.man.test.log10_bar || -6) * 2.5)]).range([padL, W - padR]);
  if (mode === "broad_sample") { const d = D._bsDom || (D._bsDom = (() => { let m = 2; const v = D.cols.bs_stat; if (v) for (let q = 0; q < v.length; q++) { const x = v[q]; if (Number.isFinite(x) && x > 0) { const a = Math.abs(Math.log2(x)); if (a > m) m = a; } } return Math.min(m, 10); })());
    return d3.scaleLinear().domain([-d, d]).range([padL, W - padR]); }
  return d3.scaleSymlog().constant(1).domain(D.exDomain).range([padL, W - padR]);
}
function drawResidualRow(cv, kind, u, W, H, opts) {
  const D = App.D; const mode = opts.mode || "close_relative"; const dpr = devicePixelRatio || 1; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + "px"; cv.style.height = H + "px";
  const c = cv.getContext("2d"); c.scale(dpr, dpr); c.clearRect(0, 0, W, H);
  const xs = residualScale(D, mode, W, opts.padL || 6, opts.padR || 6);
  const pts = []; let identical = 0; const N = kind === "locus" ? D.n : D.nloc;
  for (let k = 0; k < N; k++) {
    const g = kind === "locus" ? k : u, li = kind === "locus" ? u : k; const i = g * D.nloc + li; if (D.cols.n_comp[i] === 0) continue;
    const x = residualOf(D, g, li, mode); if (x === null || !Number.isFinite(x)) continue;
    if (mode === "close_relative") { const mi = D.num("med_identity", g, li); if (mi !== null && 100 - mi <= 0) { identical++; continue; } }
    pts.push({ g, li, x, code: D.state[i], col: ribbonColour(D, g, li, mode), px: xs(x), py: 0 });
  }
  // the saturation band for a locus: where a saturated sequence would land, from the plateau, the rate and the expectation
  if (mode === "close_relative" && kind === "locus" && !opts.noBand && D.loci.rate[u]) {
    const P = D.loci.plateau[u], rate = D.loci.rate[u]; const corr = p => -P * Math.log(1 - Math.min(Math.max(p, 0), 0.98 * P) / P);
    const ce = pts.map(p => { const me = D.num("med_expected", p.g, p.li); return me === null ? NaN : (corr(0.98 * P) - corr(100 - me)) / rate; }).filter(Number.isFinite).sort(d3.ascending);
    if (ce.length) { const c10 = d3.quantileSorted(ce, 0.1); c.fillStyle = cssVar("--surface-3"); c.fillRect(xs(c10), 1, Math.max(0, W - xs(c10)), H - 2); }
  }
  c.strokeStyle = cssVar("--model"); c.lineWidth = 1; c.beginPath(); c.moveTo(xs(0) + 0.5, 1); c.lineTo(xs(0) + 0.5, H - 1); c.stroke();
  const jit = d3.randomUniform.source(d3.randomLcg(0.31 + (u % 97) / 100))(-(H / 2 - 4), H / 2 - 4); const cy = H / 2;
  pts.sort((a, b) => (a.col ? 1 : 0) - (b.col ? 1 : 0));                     // this test's masks paint last, on top
  const many = pts.length > 400, ring = cssVar("--surface");
  const mark = opts.mark || null;
  for (const p of pts) {
    p.py = cy + jit(); const col = p.col || cssVar("--pass");
    const dim = mark && !(p.g === mark[0] && p.li === mark[1]);
    c.globalAlpha = dim ? 0.25 : p.col ? 1 : many ? 0.35 : 0.7; c.fillStyle = col; c.beginPath(); c.arc(p.px, p.py, p.col ? 2.6 : many ? 1.8 : 2.4, 0, 2 * Math.PI); c.fill(); c.globalAlpha = 1;
    if (!dim && p.col) { c.strokeStyle = ring; c.lineWidth = 0.6; c.stroke(); }   // a fine white ring on this test's masks only
  }
  if (mark) { const p = pts.find(q => q.g === mark[0] && q.li === mark[1]);
    if (p) { c.beginPath(); c.arc(p.px, p.py, 5, 0, 2 * Math.PI); c.strokeStyle = cssVar("--focus"); c.lineWidth = 1.6; c.stroke(); } }
  return { pts, identical, xs };
}
function residualAxis(W, padL, padR, mode) {
  const D = App.D; mode = mode || "close_relative"; const xs = residualScale(D, mode, W, padL, padR);
  const svg = d3.create("svg").attr("width", W).attr("height", 40).attr("class", "resid-axis");
  const S = LM_STR.ribbons;
  const caption = t => svg.append("text").attr("x", (xs.range()[0] + xs.range()[1]) / 2).attr("y", 37).attr("text-anchor", "middle").attr("fill", cssVar("--muted")).attr("font-size", 10.5).text(t);
  if (mode === "cm_relative" || mode === "cm_outside") {
    svg.append("g").attr("class", "axis").attr("transform", "translate(0,2)").call(d3.axisBottom(xs).ticks(5).tickFormat(d => lmP(-d)));
    caption(S.axis_closest);
  } else if (mode === "broad_sample") {
    svg.append("g").attr("class", "axis").attr("transform", "translate(0,2)").call(d3.axisBottom(xs).ticks(7));
    caption(S.axis_broad);
  } else {
    const tv = [-1000, -100, -10, -1, 0, 1, 10, 100, 1000].filter(v => v >= D.exDomain[0] && v <= D.exDomain[1]);
    svg.append("g").attr("class", "axis").attr("transform", "translate(0,2)").call(d3.axisBottom(xs).tickValues(tv).tickFormat(d3.format(",")));
    caption(S.axis_close_relative);
  }
  return svg.node();
}
function bindResidualHover(cv, getPts) {
  cv.addEventListener("mousemove", ev => { const r = cv.getBoundingClientRect(); const mx = ev.clientX - r.left, my = ev.clientY - r.top; let best = null, bd = 64; for (const p of getPts()) { const d = (p.px - mx) ** 2 + (p.py - my) ** 2; if (d < bd) { bd = d; best = p; } }
    if (best) { showHover(cellHover(best.g, best.li), ev.clientX, ev.clientY); App.sel.cells.clear(); App.sel.cells.add(best.g * App.D.nloc + best.li); cv.style.cursor = "pointer"; } else { hideHover(); cv.style.cursor = ""; } });
  cv.addEventListener("mouseleave", () => { hideHover(); App.sel.cells.clear(); });
  cv.addEventListener("click", ev => { const r = cv.getBoundingClientRect(); const mx = ev.clientX - r.left, my = ev.clientY - r.top; let best = null, bd = 64; for (const p of getPts()) { const d = (p.px - mx) ** 2 + (p.py - my) ** 2; if (d < bd) { bd = d; best = p; } } if (best) Sequence.open(best.g, best.li); });
}

// all loci or all taxa as residual ribbons, one row each, virtualised
const Ribbons = (() => {
  let root, kind, sortMode = "masked", mode = "close_relative", units = [], list, inner, W = 600, axisBox;
  const RH = 22, NAME_W = 190 + 4;                                               // the name column and its grid gap
  let ctlHost = null;
  function render(container, k, host) { root = container; kind = k; ctlHost = host; root.replaceChildren(); build(); }
  function build() {
    const D = App.D; const cnt = kind === "locus" ? D.locCount : D.taxCount; const names = kind === "locus" ? D.loci.names : D.taxa.names; const N = kind === "locus" ? D.nloc : D.n;
    units = d3.range(N).filter(u => cnt.tested[u] > 0);
    if (sortMode === "masked") units.sort((a, b) => cnt.masked[b] - cnt.masked[a] || cnt.near[b] - cnt.near[a] || a - b); else if (sortMode === "name") units.sort((a, b) => names[a].localeCompare(names[b]));
    const sortSel = el("select", { onchange: e => { sortMode = e.target.value; build(); } }, ...[["masked", "most masked first"], ["name", "name"], ["input", "input order"]].map(([v, t]) => el("option", { value: v, selected: v === sortMode ? "" : null }, t)));
    const SR = LM_STR.ribbons;
    const testSel = el("select", { onchange: e => { mode = e.target.value; build(); } }, ...[["close_relative", SR.test_close_relative], ["cm_relative", SR.test_cm_relative], ["cm_outside", SR.test_cm_outside], ["broad_sample", SR.test_broad_sample]].map(([v, t]) => el("option", { value: v, selected: v === mode ? "" : null }, t)));
    if (ctlHost) { ctlHost.replaceChildren(el("label", null, SR.test_label + " ", testSel), el("label", null, "sort ", sortSel), el("span", { class: "k" }, `${fmtInt(units.length)} ${kind === "locus" ? "loci" : "taxa"}`)); }
    axisBox = el("div", { class: "rib-axis" }); list = el("div", { class: "rib-list" }); inner = el("div", { style: "position:relative" }); list.append(inner);
    root.replaceChildren(axisBox, list);
    W = Math.max(300, root.getBoundingClientRect().width - NAME_W - 24);
    axisBox.replaceChildren(el("div", { class: "rib-head", style: `width:${NAME_W}px; flex: 0 0 auto` }, el("span", null, kind === "locus" ? "locus" : "taxon")), residualAxis(W, 8, 8, mode));
    inner.style.height = units.length * RH + "px";
    list.addEventListener("scroll", paint); paint();
  }
  function paint() {
    const D = App.D; const cnt = kind === "locus" ? D.locCount : D.taxCount; const names = kind === "locus" ? D.loci.names : D.taxa.names;
    inner.replaceChildren(); const top = list.scrollTop; const i0 = Math.max(0, Math.floor(top / RH) - 3), i1 = Math.min(units.length, i0 + Math.ceil(list.clientHeight / RH) + 6);
    for (let i = i0; i < i1; i++) {
      const u = units[i]; const cv = el("canvas"); let drawn = null;
      const nm = el("span", { class: "nm", title: `${names[u]}: ${cnt.masked[u]} masked of ${cnt.tested[u]} tested`, onclick: () => Detail.open(kind, u) }, sw(cnt.masked[u] ? (cnt.div[u] >= cnt.sim[u] ? D.col.red : D.col.blue) : null), names[u]);
      const row = el("div", { class: "rib-row", style: `position:absolute; top:${i * RH}px; left:0; right:0; height:${RH}px` }, nm, cv);
      inner.append(row);
      drawn = drawResidualRow(cv, kind, u, W, RH - 2, { padL: 8, padR: 8, mode });
      if (drawn.identical) row.title = `${drawn.identical} sequences identical to their comparators are not drawn`;
      bindResidualHover(cv, () => drawn.pts);
    }
  }
  return { render };
})();

// ------------------------------------------------------------------ the taxon / locus panel: one plot at a time, sized to the panel
const Detail = (() => {
  let cur = null, view = "scatter", viewSel = null, body = null, note = null;      // note: the explanation of the plot on show, behind the ? mark
  function ensurePanel() {
    if (App.panels.detail) return;
    const P = makePanel("detail", [], () => closePanel("detail")); App.panels.detail = { el: P.el, P, render, redraw: render };
    App.topEl.append(App.vsplitEl, P.el); App.vsplitEl.style.display = "";
    new ResizeObserver(debounce(() => { if (App.panels.detail && cur) draw(); }, 120)).observe(P.body);
  }
  function open(kind, i) {
    cur = { kind, i }; ensurePanel();
    if (kind === "locus" && view === "matches") view = "scatter";
    render(); writeHash(); Overview.highlight(kind === "taxon" ? "row" : "col", i);
  }
  // the copies' site ribbons in blocks by donor, recipient or pair (mode); key: one block alone, null for all
  function openCopies(mode, key) { cur = { kind: "copies", mode, key: key === undefined ? null : key }; ensurePanel(); render(); writeHash(); }
  function cellsOf() {
    const D = App.D; const out = [];
    if (cur.kind === "taxon") { const g = cur.i; for (let li = 0; li < D.nloc; li++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    else { const li = cur.i; for (let g = 0; g < D.n; g++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    return out;
  }
  function area() { const r = body.getBoundingClientRect(); return { W: Math.max(300, r.width - 20), H: Math.max(220, r.height - 20) }; }
  function render() {
    if (!cur || !App.panels.detail) return;
    if (cur.kind === "copies") return renderCopies();
    const D = App.D; const { kind, i } = cur; const P = App.panels.detail.P; const name = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
    const cnt = kind === "taxon" ? D.taxCount : D.locCount;
    const views = [["scatter", "observed vs expected"], ...(kind === "taxon" ? [["matches", "best matches"]] : []), ["ribbons", "site ribbons"]];
    viewSel = el("select", { onchange: e => { view = e.target.value; draw(); } }, ...views.map(([v, t]) => el("option", { value: v, selected: v === view ? "" : null }, t)));
    const q = el("span", { class: "help-q", title: "what this plot shows", onmouseenter: ev => { if (note) showHover(note, ev.clientX, ev.clientY); }, onmouseleave: hideHover }, "?");
    P.setCrumb([el("span", null, kind, " ", el("b", null, name)), el("span", { class: "sub" }, `${fmtInt(cnt.masked[i])} of ${fmtInt(cnt.tested[i])} sequences masked (${pct(cnt.masked[i], cnt.tested[i])})`), el("label", { class: "crumb-ctl" }, "view ", viewSel, q, el("button", { class: "cardbtn", onclick: () => Card.open(kind, i) }, "⎙ " + LM_STR.card.button))]);
    P.body.replaceChildren(); body = P.body; body.classList.add("detail-body");
    draw();
  }
  function renderCopies() {
    const D = App.D; const P = App.panels.detail.P; const n = copiesOf(D).length;
    const modeSel = el("select", { onchange: e => openCopies(e.target.value, null) }, ...[["donor", "by donor"], ["recipient", "by recipient"], ["pair", "by donor-recipient pair"]].map(([v, t]) => el("option", { value: v, selected: v === cur.mode ? "" : null }, t)));
    const q = el("span", { class: "help-q", title: "what this shows", onmouseenter: ev => { if (note) showHover(note, ev.clientX, ev.clientY); }, onmouseleave: hideHover }, "?");
    const one = cur.key === null ? null : copyGroups(D, cur.mode).find(gr => gr.k === String(cur.key));
    P.setCrumb([el("span", null, "copies", one ? [" ", el("b", null, one.label)] : ""), el("span", { class: "sub" }, one ? `${one.list.length} of ${fmtInt(n)} copies` : `${fmtInt(n)} copies`), el("label", { class: "crumb-ctl" }, "blocks ", modeSel, q)]);
    P.body.replaceChildren(); body = P.body; body.classList.add("detail-body");
    draw();
  }
  // the site ribbons of the copies (the recipients' sequences) in blocks, most copies first
  function copiesRibbons() {
    const D = App.D; const box = el("div", { class: "plot ribbons-plot" }); const nm = D.taxa.names, ln = D.loci.names;
    const groups = copyGroups(D, cur.mode).filter(gr => cur.key === null || gr.k === String(cur.key));
    if (!groups.length) { box.append(el("div", { class: "note" }, "no copies here")); return box; }
    const maxLen = d3.max(groups, gr => d3.max(gr.list, c => D.loci.length[c.li]));
    const [leftHead, rightHead] = cur.mode === "donor" ? ["recipient", "locus"] : cur.mode === "recipient" ? ["locus", "donor"] : ["locus", ""];
    const wrap = el("div", { class: "ribbons-wrap" }); box.append(wrap);
    wrap.append(el("div", { class: "rib-head2" }, el("span", null, leftHead), el("span", null, ""), el("span", null, rightHead)));
    for (const gr of groups) {
      wrap.append(el("div", { class: "rule", title: cur.mode === "pair" ? "donor → recipient" : cur.mode }, `${gr.label} · ${gr.list.length} cop${gr.list.length > 1 ? "ies" : "y"}`));
      const rows = gr.list.slice().sort((a, b) => cur.mode === "donor" ? nm[a.rec].localeCompare(nm[b.rec]) || a.li - b.li : a.li - b.li);
      for (const c of rows) {
        const { rec: g, li } = c; const cv = el("canvas");
        const [left, right] = cur.mode === "donor" ? [nm[g], ln[li]] : cur.mode === "recipient" ? [ln[li], nm[c.don]] : [ln[li], ""];
        const row = el("div", { class: "rib-row2", title: `${nm[g]} at ${ln[li]}, copied from ${nm[c.don]} · click: the alignment`, onclick: () => Sequence.open(g, li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(g * D.nloc + li); Overview.redraw(); }, onmouseleave: () => { App.sel.cells.clear(); Overview.redraw(); } },
          el("span", { class: "nm" }, sw(D.colorOf(D.state[g * D.nloc + li])), left), cv, el("span", { class: "nm bm" }, right));
        wrap.append(row);
        D.codes(li).then(C => drawRibbon(cv, C, g, li, maxLen));
      }
    }
    note = el("div", { style: "max-width: 420px" }, el("b", null, "Copies."), " Each row is a copied sequence (the recipient, masked) at one locus, drawn along its sites: blue where it shares its donor's residue while two or more of its comparators differ, grey to red as more of its comparators differ from it. The donor is the member of the too-similar pair that fits its own relatives better. Blocks run from the most copies to the fewest. Click a row for the alignment.");
    return box;
  }
  function draw() {
    if (!cur || !App.panels.detail) return; body.replaceChildren(); note = null;
    if (cur.kind === "copies") { body.append(copiesRibbons()); return; }
    const cells = cellsOf();
    if (view === "scatter") body.append(scatterPlot(cells)); else if (view === "matches") body.append(matchesPlot(cells)); else body.append(ribbonsList(cells));
  }
  // observed vs expected divergence, log axes with a zero floor and the diagonal (the expectation); there are no effect floors
  function scatterPlot(cells) {
    const D = App.D; const box = el("div", { class: "plot" });
    const pts = cells.map(([g, li]) => { const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li); return mi === null ? null : { g, li, obs: 100 - mi, exp: 100 - me, code: D.state[g * D.nloc + li] }; }).filter(Boolean);
    if (!pts.length) { box.append(el("div", { class: "note" }, "the bundle carries no numbers for these cells")); return box; }
    const { W, H: H0 } = area(); const H = H0 - 36; const m = { l: 50, r: 14, t: 10, b: 42 };
    const posX = pts.filter(p => p.exp > 0).map(p => p.exp), posY = pts.filter(p => p.obs > 0).map(p => p.obs);
    const fx = posX.length ? d3.min(posX) / 2 : 0.01, fy = posY.length ? d3.min(posY) / 2 : 0.01;
    const X = p => p.exp > 0 ? p.exp : fx, Y = p => p.obs > 0 ? p.obs : fy;
    const xs = d3.scaleLog().domain([d3.min(pts, X) / 1.3, d3.max(pts, X) * 1.3]).range([m.l, W - m.r]); const ys = d3.scaleLog().domain([d3.min(pts, Y) / 1.3, d3.max(pts, Y) * 1.3]).range([H - m.b, m.t]);
    const wrap = el("div", { class: "scatter-wrap" }); box.append(wrap);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).ticks(5, "~g"));
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).ticks(5, "~g"));
    if (posX.length < pts.length) svg.append("text").attr("x", xs(fx)).attr("y", H - m.b + 14).attr("text-anchor", "middle").text("0");
    if (posY.length < pts.length) svg.append("text").attr("x", m.l - 6).attr("y", ys(fy) + 4).attr("text-anchor", "end").text("0");
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("expected divergence of the sequence, median over its comparators (%)");
    svg.append("text").attr("transform", `translate(13,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("observed divergence (%)");
    // the diagonal (observed = expected) and, map-style, contours at observed 100%, 200% and 500% above the expectation and
    // the expectation 100%, 200% and 500% above the observed: thin lines with the level written along them
    const seg = f => { const [x0, x1] = xs.domain(); const p0 = [Math.max(x0, ys.domain()[0] / f), Math.min(x1, ys.domain()[1] / f)]; return p0[0] < p0[1] ? [xs(p0[0]), ys(p0[0] * f), xs(p0[1]), ys(p0[1] * f)] : null; };
    const diag = seg(1); if (diag) svg.append("line").attr("x1", diag[0]).attr("y1", diag[1]).attr("x2", diag[2]).attr("y2", diag[3]).attr("stroke", cssVar("--model")).attr("stroke-width", 1.2);
    const ctG = svg.append("g").attr("class", "contours");
    for (const [f, lab] of [[2, "100%"], [3, "200%"], [6, "500%"], [1 / 2, "100%"], [1 / 3, "200%"], [1 / 6, "500%"]]) {
      const e = seg(f); if (!e) continue; const ang = Math.atan2(e[3] - e[1], e[2] - e[0]) * 180 / Math.PI;
      const at = f > 1 ? 0.78 : 0.22; const lx = e[0] + (e[2] - e[0]) * at, ly = e[1] + (e[3] - e[1]) * at;      // the label sits on the line, the line broken behind it
      ctG.append("line").attr("x1", e[0]).attr("y1", e[1]).attr("x2", e[2]).attr("y2", e[3]).attr("stroke", cssVar("--muted")).attr("stroke-width", 0.7).attr("stroke-opacity", 0.45).attr("stroke-dasharray", "9 6");
      const t = ctG.append("text").attr("transform", `translate(${lx},${ly}) rotate(${ang})`).attr("text-anchor", "middle").attr("dy", "0.35em").attr("font-size", 9).attr("fill", cssVar("--muted")).attr("fill-opacity", 0.8).attr("stroke", cssVar("--surface")).attr("stroke-width", 4).attr("paint-order", "stroke").text(lab);
      t.append("title").text(f > 1 ? `observed divergence ${lab} above the expectation (${f}× expected)` : `the expectation ${lab} above the observed divergence (observed = expected ÷ ${Math.round(1 / f)})`);
    }
    wrap.append(svg.node());
    const cv = el("canvas"); const dpr = devicePixelRatio || 1; cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px"; wrap.append(cv);
    const c = cv.getContext("2d"); c.scale(dpr, dpr);
    const ring = cssVar("--surface");
    const drawPt = (p, r) => { const x = xs(X(p)), y = ys(Y(p)); c.fillStyle = D.pointColorOf(p.code) || cssVar("--grey-1"); c.beginPath();
      c.globalAlpha = p.code >= ST.RED ? 1 : 0.8; c.arc(x, y, r, 0, 2 * Math.PI); c.fill(); c.globalAlpha = 1;
      if (p.code >= ST.RED) { c.strokeStyle = ring; c.lineWidth = 0.6; c.stroke(); } };   // one size and one colour per outcome; a fine white ring on the masked points only
    const order = pts.slice().sort((a, b) => D.sev[a.g * D.nloc + a.li] - D.sev[b.g * D.nloc + b.li]);
    const r = pts.length > 2000 ? 2.4 : pts.length > 500 ? 3.2 : 4.2;
    for (const p of order) drawPt(p, r);
    svg.style("pointer-events", "all");
    const nearest = ev => { const [mx, my] = d3.pointer(ev, svg.node()); let best = null, bd = 81; for (const p of pts) { const dx = xs(X(p)) - mx, dy = ys(Y(p)) - my; const d = dx * dx + dy * dy; if (d < bd) { bd = d; best = p; } } return best; };
    svg.on("mousemove", ev => { const p = nearest(ev); if (p) { showHover(cellHover(p.g, p.li), ev.clientX, ev.clientY); App.sel.cells.clear(); App.sel.cells.add(p.g * D.nloc + p.li); Overview.redraw(); svg.style("cursor", "pointer"); } else { hideHover(); svg.style("cursor", null); } });
    svg.on("mouseleave", () => { hideHover(); App.sel.cells.clear(); Overview.redraw(); });
    svg.on("click", ev => { const p = nearest(ev); if (p) Sequence.open(p.g, p.li); });
    note = el("div", { style: "max-width: 420px" }, el("b", null, "Observed vs expected divergence."), " Each point is one sequence: its observed divergence at the locus, the median over its comparators, against the divergence the genome predicts. Diagonal: observed = expected. The contours mark, above it, an observed divergence 100%, 200% and 500% above the expectation and, below it, an expectation 100%, 200% and 500% above the observed. Grey: not masked; red: masked, too divergent; blue: " + (D.format >= 4 ? "masked as a copy of another taxon (its source can be any taxon, not only a comparator, so a copy can sit anywhere here)" : "masked, too similar") + "; yellow: masked by both tests. Hover a point for its numbers; click it for the alignment.");
    return box;
  }
  // best matches of the masked sequences
  function matchesPlot(cells) {
    const D = App.D; const g0 = cur.i; const box = el("div", { class: "plot" });
    const close = new Set(Array.from(D.closest.subarray(g0 * 5, g0 * 5 + 5)).filter(x => x >= 0));
    const masked = [], unmasked = [];
    for (const [g, li] of cells) { const bm = D.num("best_match", g, li); if (bm === null || bm < 0) continue; ((D.cols.side[g * D.nloc + li] & 3) ? masked : unmasked).push({ g, li, bm, id: D.num("best_match_identity", g, li), code: D.state[g * D.nloc + li] }); }
    if (!masked.length) { box.append(el("div", { class: "note" }, unmasked.length ? "no masked sequence for this taxon" : "best matches were not carried for this taxon (the dataset is large: they are carried for cells with a fail only)")); return box; }
    const cm = d3.rollup(masked, v => v.length, d => d.bm); const nu = unmasked.length;
    let taxaSorted = Array.from(cm.keys()).sort((a, b) => cm.get(b) - cm.get(a) || D.taxa.names[a].localeCompare(D.taxa.names[b]));
    const { W, H: H0 } = area(); const H = H0 - 60; const MAXB = Math.max(6, Math.floor((W - 60) / 34));
    let other = null; if (taxaSorted.length > MAXB) { const rest = taxaSorted.slice(MAXB - 1); taxaSorted = taxaSorted.slice(0, MAXB - 1); other = { n: rest.length, cm: d3.sum(rest, u => cm.get(u)), set: new Set(rest) }; }
    const m = { l: 44, r: 10, t: 14, b: Math.min(130, H * 0.42) };
    const cats = taxaSorted.map(u => ({ key: u, name: D.taxa.names[u], cm: cm.get(u), rows: masked.filter(d => d.bm === u) })); if (other) cats.push({ key: -1, name: `other (${other.n} taxa)`, cm: other.cm, rows: masked.filter(d => other.set.has(d.bm)) });
    const xs = d3.scaleBand().domain(cats.map(c => c.name)).range([m.l, W - m.r]).padding(0.18); const top = Math.max(d3.max(cats, c => c.cm), 1);
    const ys = d3.scaleLinear().domain([0, top * 1.15]).range([H - m.b, m.t]);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).tickValues(d3.range(0, top + 1, Math.max(1, Math.ceil(top / 4)))).tickFormat(d3.format("d")));
    const ax = svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs));
    ax.selectAll("text").attr("transform", "rotate(-60)").attr("text-anchor", "end").attr("dx", "-0.6em").attr("dy", "0.2em").style("font-weight", d => cats.find(c => c.name === d) && close.has(cats.find(c => c.name === d).key) ? "700" : "400").style("fill", cssVar("--ink"));
    for (const c of cats) {
      let y0 = H - m.b; const byCode = d3.rollup(c.rows, v => v.length, d => D.famOf(d.code)); const codes = Array.from(byCode.keys()).sort(d3.ascending);   // one segment per outcome
      for (const code of codes) { const h = ys(0) - ys(byCode.get(code)); svg.append("rect").attr("x", xs(c.name)).attr("y", y0 - h).attr("width", xs.bandwidth()).attr("height", Math.max(h, 0)).attr("fill", D.colorOf(code)).attr("stroke", cssVar("--surface")).attr("stroke-width", 0.8).style("cursor", "pointer")
        .on("mousemove", ev => showHover(el("div", null, el("b", null, c.name), el("div", null, `${byCode.get(code)} masked sequence${byCode.get(code) > 1 ? "s" : ""} match it best (${D.famLabel(code)})`), el("div", { class: "more" }, `${c.cm} masked in all`)), ev.clientX, ev.clientY)).on("mouseleave", hideHover)
        .on("click", () => { App.sel.cells.clear(); for (const d of c.rows) App.sel.cells.add(d.g * D.nloc + d.li); Overview.redraw(); }); y0 -= h; }
      svg.append("text").attr("x", xs(c.name) + xs.bandwidth() / 2).attr("y", ys(c.cm) - 3).attr("text-anchor", "middle").style("fill", cssVar("--ink")).text(c.cm);
    }
    box.append(svg.node());
    const pm = 100 * d3.mean(masked, d => close.has(d.bm) ? 1 : 0), pu = nu ? 100 * d3.mean(unmasked, d => close.has(d.bm) ? 1 : 0) : NaN;
    note = el("div", { style: "max-width: 420px" }, el("b", null, "Best matches."), ` Masked sequences by the taxon each matches best at its locus, stacked by mask kind; bold names are the 5 closest relatives genome-wide. The best match is one of the 5 closest for ${fmt1(pm)}% of masked and ${Number.isFinite(pu) ? fmt1(pu) : "–"}% of unmasked sequences; median identity to it ${fmt1(d3.median(masked, d => d.id))}% masked vs ${nu ? fmt1(d3.median(unmasked, d => d.id)) : "–"}% unmasked.`, el("div", { class: "more", style: "margin-top: 4px" }, "Bars on a distant taxon point to contamination or a shared paralogue."));
    return box;
  }
  // site ribbons of the masked sequences: one pixel-row each, coloured by the site code
  function ribbonsList(cells) {
    const D = App.D; const box = el("div", { class: "plot ribbons-plot" });
    const masked = cells.filter(([g, li]) => D.cols.side[g * D.nloc + li] & 3);
    if (!masked.length) { box.append(el("div", { class: "note" }, "no masked sequence here")); return box; }
    const maxLen = cur.kind === "locus" ? D.loci.length[cur.i] : d3.max(masked, ([g, li]) => D.loci.length[li]);
    const wrap = el("div", { class: "ribbons-wrap" }); box.append(wrap);
    masked.sort((a, b) => D.sev[b[0] * D.nloc + b[1]] - D.sev[a[0] * D.nloc + a[1]]);
    wrap.append(el("div", { class: "rib-head2" }, el("span", null, cur.kind === "taxon" ? "locus" : "taxon"), el("span", null, ""), el("span", { title: "the taxon closest to the sequence at this locus (per-site divergence compared on the same sites); in the outside colour when it is not one of the sequence's comparators" }, "best match")));
    for (const [g, li] of masked) {
      const cv = el("canvas"); const bm = D.num("best_match", g, li); const inSet = bm !== null && bm >= 0 && D.comps(g, li).includes(bm);
      const bmCell = el("span", { class: "nm bm", style: bm !== null && bm >= 0 && !inSet ? `color:${TREE_COL.outside}` : null, title: bm !== null && bm >= 0 ? `${D.taxa.names[bm]}: ${fmt1(D.num("best_match_identity", g, li))}% identity over shared sites${inSet ? ", one of the comparators" : ", not one of the comparators"}` : "" }, bm !== null && bm >= 0 ? D.taxa.names[bm] : "");
      const row = el("div", { class: "rib-row2", onclick: () => Sequence.open(g, li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(g * D.nloc + li); Overview.redraw(); }, onmouseleave: () => { App.sel.cells.clear(); Overview.redraw(); } },
        el("span", { class: "nm", title: cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g] }, sw(D.colorOf(D.state[g * D.nloc + li])), cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g]), cv, bmCell);
      wrap.append(row);
      D.codes(li).then(C => drawRibbon(cv, C, g, li, maxLen));
    }
    note = el("div", { style: "max-width: 420px" }, el("b", null, "Site ribbons."), ` ${fmtInt(masked.length)} masked sequences, sites along x: red where a majority of the comparators differ, blue ${D.format >= 4 ? "where a copy shares its source's state while two or more of its comparators differ" : "where a pair shares a state against the others"}. ` + (cur.kind === "locus" ? "Rows share the alignment's columns: a vertical band is a column problem (misalignment), a horizontal block in one row is a taxon problem (a wrong exon, a frameshift)." : "Each row is a different locus drawn to its own length on a shared axis.") + " The right-hand column names the taxon closest to the sequence at that locus, in the outside colour when it is not one of its comparators. Hover a row for the site; click it for the alignment.");
    return box;
  }
  function drawRibbon(cv, C, g, li, maxLen) {
    const D = App.D; if (!C) { cv.title = "codes not carried for this cell"; return; }
    let r = -1; for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) { r = k; break; }
    if (r < 0) { cv.title = "codes not carried for this cell"; return; }
    const len = C.len; const W = Math.max(50, Math.floor(cv.getBoundingClientRect().width || 600)); cv.width = W; cv.height = 1;
    const img = cv.getContext("2d").createImageData(W, 1); const px = img.data; const base = r * len; const surf = hexToRgb(cssVar("--surface"));
    const span = len / maxLen;
    for (let x = 0; x < W; x++) {
      const s0 = Math.floor(x / W / span * len), s1 = Math.max(s0, Math.ceil((x + 1) / W / span * len) - 1);
      let rgb = surf;
      if (s0 < len) { let best = -1, code = 255; for (let s = s0; s <= Math.min(s1, len - 1); s++) { const c = C.blk[base + s]; const sv = codeSeverity(D, c); if (sv > best) { best = sv; code = c; } } rgb = codeRgb(D, code); }
      px[x * 4] = rgb[0]; px[x * 4 + 1] = rgb[1]; px[x * 4 + 2] = rgb[2]; px[x * 4 + 3] = 255;
    }
    cv.getContext("2d").putImageData(img, 0, 0);
    cv.style.width = "100%"; cv.style.height = "10px";
    cv.addEventListener("mousemove", ev => { const rr = cv.getBoundingClientRect(); const s = Math.floor((ev.clientX - rr.left) / rr.width / span * len); if (s >= 0 && s < len) showHover(el("div", null, el("b", null, `${D.taxa.names[g]} · ${D.loci.names[li]}`), el("div", null, `site ${s + 1} of ${fmtInt(len)} (alignment position ${fmtInt(alignmentPos(D, li, s) + 1)}): ${codeLabel(D, C.blk[base + s])}`)), ev.clientX, ev.clientY); });
    cv.addEventListener("mouseleave", hideHover);
  }
  return { open, openCopies, render, cur: () => cur, setView: v => { view = v; } };
})();

// site-code helpers shared by the ribbons and the alignment panel
function codeSeverity(D, c) { if (c === 255) return 0; if (c === 254) return 5; if (c >= 128) return 30 + (c - 128); if (c >= D.M) return 20 + c; return 2 + c; }
function codeHex(D, c) {                                    // the alignment's grades by k; blue where the partner agrees with the focal
  if (c === 255) return null; if (c === 254) return cssVar("--muted");
  if (c >= 128) return D.col.blueTint;
  return D.grades[Math.min(c, D.K)];
}
function codeRgb(D, c) { const h = codeHex(D, c); return hexToRgb(h || cssVar("--surface")); }
function codeLabel(D, c) { if (c === 255) return "no residue in the focal sequence"; if (c === 254) return D.format >= 4 ? "the copy differs from its source here" : "the pair differs here"; if (c >= 128) return D.format >= 5 ? `the partner agrees with the focal here while ${c - 128} comparators differ` : D.format >= 4 ? `the copy shares its source's state while ${c - 128} of its comparators differ` : `the pair shares this state while ${c - 128} other comparators differ`; return `${c} of the comparators present differ${c >= D.M ? " (a majority)" : ""}`; }
function alignmentPos(D, li, s) { let left = s; for (const [a, b, step] of D.loci.ranges[li]) { const n = Math.floor((b - a) / step) + 1; if (left < n) return a + left * step; left -= n; } return s; }

// ------------------------------------------------------------------ the two distance trees of one cell as a cophyloplot (drawn from the tool's Newick strings)
function parseNewick(s) {
  let i = 0; const node = () => { const n = { children: [], length: 0, name: null };
    if (s[i] === "(") { i++; n.children.push(node()); while (s[i] === ",") { i++; n.children.push(node()); } i++; }
    let name = ""; while (i < s.length && !",:);".includes(s[i])) name += s[i++]; n.name = name || null;
    if (s[i] === ":") { i++; let num = ""; while (i < s.length && !",);".includes(s[i])) num += s[i++]; n.length = parseFloat(num) || 0; }
    return n; };
  return node();
}
function tipsOf(n, out = []) { if (!n.children.length) out.push(n); else n.children.forEach(c => tipsOf(c, out)); return out; }
function depthOf(n, d = 0) { n.depth = d; return n.children.length ? Math.max(...n.children.map(c => depthOf(c, d + c.length))) : d; }
function orderLike(n, rank) {                                     // rotate a tree's nodes to follow the other tree's tip order
  if (!n.children.length) { n.rank = rank.get(n.name) ?? 0; return n.rank; }
  const r = n.children.map(c => orderLike(c, rank)); n.children = n.children.map((c, i) => [c, r[i]]).sort((a, b) => a[1] - b[1]).map(x => x[0]); n.rank = d3.mean(r); return n.rank;
}
function assignY(n, counter) { if (!n.children.length) { n.y = counter.y++; return n.y; } const ys = n.children.map(c => assignY(c, counter)); n.y = d3.mean(ys); return n.y; }
function splitsOf(root) {                                         // the bipartition of every edge, keyed by the smaller side's sorted tip names
  const all = tipsOf(root).map(t => t.name).sort(); const key = set => { const a = [...set].sort(); const b = all.filter(x => !set.has(x)); return (a.length < b.length || (a.length === b.length && a.join() < b.join()) ? a : b).join("|"); };
  const walk = n => { if (!n.children.length) { n.set = new Set([n.name]); return n.set; } n.set = new Set(); for (const c of n.children) for (const x of walk(c)) n.set.add(x); n.split = key(n.set); return n.set; };
  walk(root); return root;
}
const TREE_COL = { outside: "#c2185b" };                          // a best match from outside the comparator set, in the ribbons
function drawCophylo(host, D, entry, k, cell, T) {      // T: cellTests, for the tips' colours
  const { g, li } = cell; const exp = parseNewick(entry.expected[k]), obs = parseNewick(entry.observed[k]);
  const name = x => D.taxa.names[+x]; const partner = D.num("partner", g, li), bm = D.num("best_match", g, li), bmId = D.num("best_match_identity", g, li); const comps = new Set(D.comps(g, li));
  const v4 = D.format >= 4, sr = similarRole(D, g, li), nsim = D.cols.n_similar[g * D.nloc + li];   // from format 4 the trees carry up to five too-similar partners, not the best match
  const sigPair = !!(sr && (sr.fired || (D.format < 5 && v4 && sr.role !== "undetermined"))), F = String(g), Pn = partner !== null && partner >= 0 ? String(partner) : null;
  // ladderise the expected tree, then rotate the observed tree to follow it; the roots come from the tool (on the
  // outgroup when the bundle carries one, drawn first with a shortened, broken branch; midpoint otherwise)
  const ogPos = entry.outgroup ? entry.outgroup[k] : -1; const ogName = ogPos >= 0 ? String(entry.tips[k][ogPos]) : null;
  const isOg = n => !n.children.length && n.name === ogName; const ogFirst = n => { if (ogName) n.children.sort((a, b) => isOg(b) - isOg(a)); };
  const ladder = n => { n.children.forEach(ladder); n.children.sort((a, b) => tipsOf(a).length - tipsOf(b).length); }; ladder(exp); ogFirst(exp);
  const rank = new Map(); tipsOf(exp).forEach((t, i) => rank.set(t.name, i)); orderLike(obs, rank); ogFirst(obs);
  splitsOf(exp); splitsOf(obs); const expSplits = new Set(), obsSplits = new Set();
  const collect = (n, out) => { if (n.children.length && n.split && n.set.size > 1 && n.set.size < tipsOf(exp).length - 1) out.add(n.split); n.children.forEach(c => collect(c, out)); }; collect(exp, expSplits); collect(obs, obsSplits);
  const labelOf = t => name(t);
  // a tip's hover: inside or outside the focal group (the comparators) and, for the closest match at the locus, whether it is
  // significantly closer than expected. The closest match: the partner when it was tested as the closest match (from format 5,
  // a partner outside the comparators), else the best match
  const closestTip = Pn !== null && D.format >= 5 && !comps.has(partner) ? Pn : bm !== null && bm >= 0 ? String(bm) : Pn;
  const tipTitle = t => +t === g ? `${name(t)}: focal sequence` : `${name(t)}: ${comps.has(+t) ? "inside" : "outside"} focal group${t === closestTip ? `; closest match at locus; ${sigPair && t === Pn ? "" : "not "}significantly closer than expected` : ""}`;
  const clos = entry.closest ? entry.closest[k] : null;                                                       // [expected, observed]: the tip closest to the focal in each tree, by that tree's own distances
  const closName = j => clos && clos[j] >= 0 ? String(entry.tips[k][clos[j]]) : null;
  const meas = document.createElement("canvas").getContext("2d"); meas.font = `bold 11.5px ${cssVar("--sans")}`; const LAB = Math.ceil(d3.max(tipsOf(exp), t => meas.measureText(labelOf(t.name)).width) + 10);
  const nTip = tipsOf(exp).length; const RH = 20, GAP = 44, W = Math.min(820, Math.max(560, 2 * LAB + GAP + 240)), top = 10, H = top + nTip * RH + 40;
  depthOf(exp); depthOf(obs); const inDepth = root => d3.max(tipsOf(root).filter(t => !isOg(t)), t => t.depth) || 0;
  // one scale for both trees; each takes the width its own depth needs, so a shallow expected tree leaves no blank
  const dE = inDepth(exp), dO = inDepth(obs), maxD = Math.max(dE, dO, 1e-6), OGLEN = 22;
  const sx = d3.scaleLinear().domain([0, Math.max(dE + dO, 1e-6)]).range([0, W - 2 * LAB - GAP - 28 - 2 * OGLEN]);
  assignY(exp, { y: 0 }); assignY(obs, { y: 0 });
  const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "tangle");
  const yOf = n => top + n.y * RH + RH / 2;
  const focalCol = cssVar("--orange"), partnerCol = D.col.blue, inCol = cssVar("--ink"), outCol = cssVar("--faint");
  const stOf = t => T.blue.has(+t) ? "blue" : T.red.has(+t) ? "red" : null;      // the colour the alignment gives the taxon's name
  const draw = (root, x0, dir, otherSplits, bold) => {           // dir +1 draws left to right, -1 right to left; bold: the tip closest to the focal in this tree
    const gg = svg.append("g"); const walk = n => { for (const c of n.children) { const og = isOg(c); const xa = x0 + dir * sx(n.depth), xb = og ? x0 + dir * OGLEN : x0 + dir * sx(c.depth); c.xdraw = xb;
      const conflict = c.children.length && c.split && c.set.size > 1 && c.set.size < nTip - 1 && !otherSplits.has(c.split);
      const onPath = dir < 0 && sigPair && Pn !== null && c.set && (c.set.has(F) !== c.set.has(Pn));   // observed tree: the branches between the focal and its partner
      gg.append("path").attr("d", `M${xa},${yOf(n)}V${yOf(c)}H${xb}`).attr("fill", "none").attr("stroke", onPath ? partnerCol : conflict ? cssVar("--ink") : cssVar("--muted")).attr("stroke-width", onPath ? 3 : conflict ? 2.6 : 1.2).append("title").text(conflict ? "this split is not in the other tree" : og ? `the outgroup's branch, ${fmt1(c.length)}%, drawn shortened` : `branch length ${fmt1(c.length)}%`);
      if (og) { const xm = (xa + xb) / 2, y = yOf(c);                                                     // the break marks on the shortened branch
        gg.append("line").attr("x1", xm - 2.5).attr("x2", xm + 2.5).attr("y1", y + 4).attr("y2", y - 4).attr("stroke", cssVar("--surface")).attr("stroke-width", 4);
        gg.append("line").attr("x1", xm - 3.5).attr("x2", xm + 0.5).attr("y1", y + 4).attr("y2", y - 4).attr("stroke", cssVar("--muted")).attr("stroke-width", 1.2);
        gg.append("line").attr("x1", xm + 0.5).attr("x2", xm + 4.5).attr("y1", y + 4).attr("y2", y - 4).attr("stroke", cssVar("--muted")).attr("stroke-width", 1.2); }
      if (!c.children.length) { const isF = +c.name === g, isP = +c.name === partner, inSet = comps.has(+c.name), isBold = c.name === bold, st = isF ? null : stOf(c.name);
        const dot = isF ? focalCol : st === "blue" ? partnerCol : st === "red" ? D.col.red : inSet || isP ? inCol : outCol; const hollow = isP && st !== "blue";   // a closest match that is not too close: hollow
        if (og) gg.append("rect").attr("x", xb - 3).attr("y", yOf(c) - 3).attr("width", 6).attr("height", 6).attr("fill", cssVar("--surface")).attr("stroke", cssVar("--muted")).attr("stroke-width", 1.2);
        else gg.append("circle").attr("cx", xb).attr("cy", yOf(c)).attr("r", hollow ? 3.4 : isP ? 4 : 3.2).attr("fill", hollow ? cssVar("--surface") : dot).attr("stroke", dot).attr("stroke-width", hollow ? 1.6 : 1);
        const lab = gg.append("text").attr("x", xb + dir * 6).attr("y", yOf(c) + 4).attr("text-anchor", dir > 0 ? "start" : "end").attr("font-size", 11.5).attr("fill", og ? cssVar("--muted") : isF ? cssVar("--orange-ink") : st === "blue" ? partnerCol : st === "red" ? D.col.red : inSet || isP ? cssVar("--ink") : outCol).attr("font-weight", (isF || isBold || st) && !og ? 700 : 400).text(labelOf(c.name));
        lab.append("title").text(tipTitle(c.name)); }
      walk(c); } }; walk(root); return gg; };
  draw(exp, 8, 1, obsSplits, closName(0)); draw(obs, W - 8, -1, expSplits, closName(1));
  const obsTips = new Map(); tipsOf(obs).forEach(t => obsTips.set(t.name, t));
  for (const t of tipsOf(exp)) { const o = obsTips.get(t.name); if (!o) continue;      // the connector runs through the gap between the two label columns
    const x1 = t.xdraw + 6 + meas.measureText(labelOf(t.name)).width + 4, x2 = o.xdraw - 6 - meas.measureText(labelOf(t.name)).width - 4;
    svg.append("line").attr("x1", x1).attr("x2", x2).attr("y1", yOf(t)).attr("y2", yOf(o)).attr("stroke", cssVar("--floor")).attr("stroke-width", 1).attr("stroke-dasharray", "2 3"); }
  const bar = (() => { const lim = 0.28 * (W - 2 * LAB - GAP - 28 - 2 * OGLEN); let best = null;             // the largest 1, 2 or 5 x 10^k that fits in a quarter of the tree width
    for (let k = -2; k <= 3; k++) for (const f of [1, 2, 5]) { const v = f * Math.pow(10, k); if (sx(v) <= lim) best = v; } return best || 1; })();
  const by = H - 12; svg.append("line").attr("x1", 8).attr("x2", 8 + sx(bar)).attr("y1", by).attr("y2", by).attr("stroke", cssVar("--ink")).attr("stroke-width", 1.5);
  svg.append("text").attr("x", 12 + sx(bar)).attr("y", by + 4).attr("font-size", 10.5).attr("fill", cssVar("--muted")).text(`${bar}% ${D.format >= 3 ? "corrected divergence" : "divergence"}`);
  svg.append("text").attr("x", W - 36).attr("y", by + 4).attr("text-anchor", "end").attr("font-size", 10.5).attr("fill", cssVar("--muted")).text("expected  ·  observed at this locus");
  const q = svg.append("g").attr("transform", `translate(${W - 20},${by})`).style("cursor", "help");
  q.append("circle").attr("r", 8).attr("fill", "none").attr("stroke", cssVar("--muted")).attr("stroke-width", 1.2); q.append("text").attr("y", 4).attr("text-anchor", "middle").attr("font-size", 11).attr("font-weight", 700).attr("fill", cssVar("--muted")).text("?");
  const explain = el("div", { style: "max-width: 380px" },
    el("div", null, el("b", null, "Two neighbour-joining trees"), (D.format >= 5 ? ` over the focal sequence (orange), its comparators (black; red when significantly too divergent) and, when it has one, its too-similar partner or its closest match at the locus (blue when significantly too similar; a hollow dot when not${sr ? `; ${D.taxa.names[sr.pt]} is ${sr.note}` : ""}); the colours are the alignment's, other taxa are grey. Hover a tip for its role.` : v4 ? ` over the focal sequence (orange), its comparators (black; red when significantly too divergent) and up to five of the taxa it is significantly too similar to (blue), its reported partner first (${nsim ? fmtInt(nsim) : "none"} here${sr ? `; ${D.taxa.names[sr.pt]} is ${sr.note}` : ""}). Hover a tip for its role.` : " over the focal sequence (orange), its comparators (black dots), the partner of a too-similar pair (blue label) and, when it is not one of the comparators, the best match (crimson dot and label; the taxon closest to the focal by per-site divergence compared on the same sites" + (bmId !== null ? `, identity ${fmt1(bmId)}%` : "") + ").") + " In each tree the taxon closest to the focal by that tree's own distances is in bold: the closest relative genome-wide on the left, the closest at this locus on the right."),
    el("div", { style: "margin-top: 4px" }, el("b", null, "Left:"), (D.format >= 3 ? " the distances the test expects for every pair: each pair's genome-wide divergence, F81-corrected with this locus's plateau and scaled by its rate. " : " the distances the test expects for every pair, each pair's genome-wide divergence through this locus's expectation. "), el("b", null, "Right:"), (D.format >= 3 ? " the divergence observed at this locus over shared sites, corrected the same way (capped at 0.98 of the plateau). " + (ogName ? "Both rooted on the same outgroup (hollow square, muted): the taxon at this locus most distant from these taxa genome-wide; its branch is drawn shortened, with a break. One scale over the ingroup, in corrected divergence; negative branches drawn as zero." : "Both midpoint-rooted, one scale in corrected divergence; negative branches drawn as zero.") : " the mismatch fractions observed at this locus over shared sites. Both midpoint-rooted, one scale; negative branches drawn as zero.")),
    el("div", { style: "margin-top: 4px" }, "Thick branches are splits absent from the other tree; dotted lines join the same taxon on both sides. A black closest taxon on the right means the sequence has moved away from everything; a blue one means it has moved towards a taxon it is significantly too similar to."),
    el("div", { style: "margin-top: 4px", class: "more" }, v4 ? "A long focal branch with its neighbours unchanged is a damaged copy of the right sequence. A focal that moves across to a too-similar partner points to contamination or a paralogue; the direction test says which of the two is the copy. A whole clade pushed away from the rest on one long branch is a paralogue shared by the clade." : "A long focal branch with its neighbours unchanged is a damaged copy of the right sequence. A focal that moves across to its best match points to contamination or a paralogue. A whole clade pushed away from the rest on one long branch is a paralogue shared by the clade. A best match on its own long branch is merely the least different taxon, not a relative."));
  q.on("mouseenter", ev => showHover(explain, ev.clientX, ev.clientY)).on("mouseleave", hideHover);
  return svg.node();
}

// ------------------------------------------------------------------ the alignment panel (full width, bottom)
const Sequence = (() => {
  let cur = null, pps = null, x0 = 0, dots = true, drawFn = null, showTrees = false;
  const RULER_H = 18, ROW_H = 18, RULE_H = 7;                                       // the rule between the similarity rows and the divergence block
  const colW = { role: 72, name: 150 }; try { Object.assign(colW, JSON.parse(localStorage.getItem("lm.seqcols") || "{}")); } catch (e) { /* no storage */ }
  function open(g, li) {
    cur = { g, li }; pps = null; x0 = 0;
    if (!App.panels.sequence) {
      const P = makePanel("sequence", [], () => closePanel("sequence")); App.panels.sequence = { el: P.el, P, render, redraw: () => drawFn && drawFn() };
      App.bottomEl.append(P.el); App.bottomEl.style.display = ""; App.splitEl.style.display = ""; if (!App.bottomEl.style.flexBasis) App.bottomEl.style.flexBasis = "42%";
    }
    App.sel.cells.clear(); App.sel.cells.add(g * App.D.nloc + li); Overview.redraw(); render(); writeHash();
  }
  const isGap = b => b === 45 || b === 63 || b === 46 || b === 42 || b === 32;      // - ? . * space
  function inkFor(hex) { const L = luminance(hex); return (L + 0.05) / 0.05 >= 1.05 / (L + 0.05) ? "#0b0b0b" : "#ffffff"; }   // the letter colour with the higher contrast
  async function render() {
    if (!cur || !App.panels.sequence) return;
    const D = App.D; const { g, li } = cur; const P = App.panels.sequence.P; const i = g * D.nloc + li;
    P.setCrumb([D.man.dataset || D.man.input.file, el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("taxon", g); } }, D.taxa.names[g]), el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("locus", li); } }, D.loci.names[li])]);
    P.body.replaceChildren(el("div", { class: "loading" }, "loading the sequences…"));
    const [C, S, SIM] = await Promise.all([D.codes(li), D.seq(li), D.similar(li)]);
    if (!cur || cur.g !== g || cur.li !== li) return;
    const comps = D.comps(g, li); const partner = D.num("partner", g, li); const sd = D.cols.side[i]; const sr = similarRole(D, g, li);
    const T = cellTests(D, g, li, SIM, sr); let drawNow = null; const redrawNow = () => { if (drawNow) drawNow(); };
    // the rows, the same shape for every cell: the focal; its closest match at this locus, the candidate of the similarity
    // test; a rule; its relatives in genome-wide rank order, the divergence test. A taxon in both roles appears twice
    const Y = T.closest !== null ? T.closest.taxon : null, ord = ["closest", "2nd closest", "3rd closest", "4th closest", "5th closest"];
    const rcr0 = D.rcAt ? D.rcAt(i) : null; const bsRep = rcr0 && rcr0.pt >= 0 ? rcr0.pt : null;   // the sample's biggest-residual sequence
    const rows = [{ g, role: "focal", tag: "focal", tip: "the focal sequence" },
      { g: Y, role: "closest", tag: "closest", tip: Y === null ? (D.format >= 7 ? "no taxon shares enough sites with the focal at this locus" : "this report does not record the closest match of a sequence whose closest-match test did not run (it is carried from report format 7)") : `the focal's closest match at this locus${comps.includes(Y) ? ` (also relative ${comps.indexOf(Y) + 1})` : " (not one of its relatives)"}` },
      { g: bsRep, role: "broad", tag: "comparator", tip: bsRep === null ? "the report carries no broad-sample representative for this cell" : (sr && sr.rc ? sr.note : "the broad comparator with the most extreme standardised residual") }];
    comps.forEach((c, r) => rows.push({ g: c, role: "comp", tag: `relative ${r + 1}`, tip: `relative ${r + 1}: the ${ord[r] || `${r + 1}th closest`} taxon genome-wide` }));
    const divAt = rows.findIndex(r => r.role === "comp"), broadAt = rows.findIndex(r => r.role === "broad");
    const ruleAt = new Set([1, broadAt, divAt]);                                  // a rule between each test's block
    const rulesBefore = ri => { let n = 0; for (const q of ruleAt) if (q <= ri) n++; return n; };
    const seqRow = gg => { if (!S) return null; for (let k = 0; k < S.taxa.length; k++) if (S.taxa[k] === gg) return S.blk.subarray(k * S.len, (k + 1) * S.len); return null; };
    let codeRow = null; if (C) for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) codeRow = C.blk.subarray(k * C.len, (k + 1) * C.len);
    const len = D.loci.length[li]; for (const r of rows) r.bytes = r.g === null ? null : seqRow(r.g);
    const focal = rows[0].bytes;
    for (const r of rows.slice(1)) { let n = 0, d = 0; if (r.bytes && focal) for (let s = 0; s < len; s++) { const a = focal[s], b = r.bytes[s]; if (isGap(a) || isGap(b)) continue; n++; if (a !== b) d++; } r.shared = n; r.diff = n ? d / n : NaN; }
    // each row's colour is its own test's: the closest row blue when significantly too similar to the focal, a relative red when
    // significantly too divergent from it (the trees colour a taxon by whichever test fired: blue, red or black)
    const s2v0 = D.s2 ? D.s2(i) : 0; const bsC0 = (s2v0 & 4) && !(s2v0 & 32), bsD0 = (s2v0 & 8) && !(s2v0 & 32);
    for (const r of rows) r.st = r.role === "focal" ? "focal" : r.g === null ? null : r.role === "closest" ? (T.blue.has(r.g) ? "blue" : null) : r.role === "broad" ? (bsC0 || s2v0 & 64 ? "blue" : bsD0 ? "red" : null) : (T.red.has(r.g) ? "red" : null);
    // per site: n, the relatives present (the closest match is never counted), and d, those that differ from the focal. A red
    // site has n >= 3 and d >= half of n, rounded up; a blue site is a red site where the closest match agrees with the focal.
    // From the rows when they are carried, else roughly from the tool's site codes
    const countRows = rows.filter(r => r.role === "comp" && r.g !== Y && r.bytes), candRow = rows[1].bytes ? rows[1] : null;
    const nRow = new Uint8Array(len), dRow = new Int16Array(len).fill(-1), siteRed = new Uint8Array(len), siteBlue = new Uint8Array(len);
    if (focal) for (let s = 0; s < len; s++) {
      if (isGap(focal[s])) continue; let n = 0, d = 0;
      for (const r of countRows) { const b = r.bytes[s]; if (isGap(b)) continue; n++; if (b !== focal[s]) d++; }
      nRow[s] = n; dRow[s] = d; siteRed[s] = n >= 3 && d >= Math.ceil(n / 2) ? 1 : 0;
      if (siteRed[s] && candRow && candRow.bytes[s] === focal[s]) siteBlue[s] = 1;
    }
    else if (codeRow) for (let s = 0; s < len; s++) { const cc = codeRow[s]; if (cc === 255 || cc === 254) continue; const k = cc >= 128 ? cc - 128 : cc; nRow[s] = D.K; dRow[s] = k; siteRed[s] = k >= D.M ? 1 : 0; siteBlue[s] = cc >= 128 && k >= D.M ? 1 : 0; }
    // header: the verdict in the fewest words, with the stars; the numbers are in the hovers and the output files
    const head = el("div", { class: "seq-head" });
    const Pn = sr ? D.taxa.names[sr.pt] : partner !== null && partner >= 0 ? D.taxa.names[partner] : null;
    const masked = sd & 3, off = !masked && (sd & (16 | 32)), says = [];
    const divP = D.num("div_log10p", g, li), simP = D.num("sim_log10p", g, li), withStars = (txt, p) => T.stars(p) ? `${txt} ${T.stars(p)}` : txt;
    const dv = sd & (1 | 16), sm = sd & (2 | 32 | 128);
    const s2v = D.s2 ? D.s2(i) : 0, rcr = D.rcAt ? D.rcAt(i) : null;
    const divMain = D.cols.n_comp[i] === D.K ? D.cols.n_fail[i] >= D.M : D.cols.n_fail[i] >= (D.cols.n_comp[i] >> 1) + 1;
    if (dv && !((s2v & 8) && (divP === null || Number.isNaN(divP)))) says.push(withStars("significantly too divergent", divP));
    if (sr && sr.rc && Pn && (sd & 1) && (s2v & 8)) says.push(withStars(`significantly too divergent from ${Pn}, a broad-sample comparator`, rcr ? rcr.var_p : null));
    const donorTxt = D.format < 4 ? "" : sd & 8 ? " (donor uncertain: both masked)" : sd & (2 | 32) ? ` (${Pn} looks more likely to be the donor)` : sd & 128 ? ` (${D.taxa.names[g]} looks more likely to be the donor)` : "";
    if (sr && sr.rc && Pn && sm && (s2v & 4)) says.push(withStars(`significantly too close to ${Pn}, a broad-sample comparator`, rcr ? rcr.var_p : null));
    else if (Pn && (D.format < 4 || sm)) says.push(withStars(`significantly too similar to ${Pn}`, simP) + donorTxt);
    if ((s2v & 64) && Pn) says.push(withStars(`the broad sample also calls it too close to ${Pn}'s side`, rcr ? rcr.lo : null));
    const verdict = D.cols.n_comp[i] === 0 ? `Not tested: ${stateReason(D, g, li) || "no comparators"}` : `${masked ? "Masked" : off ? "Not masked (masking was off)" : "Not masked"}${says.length ? ": " + says.join("; ") : ""}${sd & 4 ? " (fewer comparators than requested)" : ""}`;
    // the tests in words, for the hovers: the focal's own two, or a row's against the focal
    const lp = (what, p) => p === null || p === undefined || Number.isNaN(p) ? `${what}: not tested` : `${what}: p ${lmP(p)}${T.stars(p) ? " " + T.stars(p) : ""}`, barTxt = T.bar !== null ? ` · bar ${lmP(T.bar)}` : "";
    const focalNote = [divP !== null ? lp("too divergent", divP) : null, simP !== null && Pn ? lp(`too similar to ${Pn}`, simP) : null].filter(Boolean).join(" · ");
    const testNote = r => {
      if (r.role === "focal") return focalNote ? focalNote + barTxt : "";
      if (r.g === null) return "";
      if (r.role === "closest") return lp("too similar to the focal", T.blue.has(r.g) ? T.blue.get(r.g).p : T.closest.p) + barTxt;
      if (r.role === "broad") return T.blue.has(r.g) ? lp("too close to the focal (the broad-sample test, stage 1)", T.blue.get(r.g).p) + barTxt : T.red.has(r.g) ? lp("too divergent from the focal (the broad-sample test, stage 1)", T.red.get(r.g)) + barTxt : "";
      const t = T.comp && T.comp.get(r.g); if (!t) return T.blue.has(r.g) ? lp("too similar to the focal", T.blue.get(r.g).p) + barTxt : "";
      return `${lp("too divergent from the focal", t.div)} · ${lp("too similar to the focal", T.blue.has(r.g) ? T.blue.get(r.g).p : t.sim)}${barTxt}`;
    };
    const parts = [el("span", { class: "state", title: focalNote ? focalNote + barTxt : null }, sw(D.colorOf(D.state[i])), verdict)];
    const dotsBox = el("input", { type: "checkbox", checked: dots ? "" : null, onchange: e => { dots = e.target.checked; redrawNow(); } });
    const treesBtn = el("button", { title: "the expected and the observed distances among these sequences as two neighbour-joining trees, face to face", onclick: () => toggleTrees() }, "trees");
    const ctl = el("span", { class: "seq-ctl" }, el("label", null, dotsBox, " dots where a sequence matches the focal residue"),
      el("button", { onclick: () => { pps = fitPps(); x0 = 0; redrawNow(); } }, "whole locus"), el("button", { onclick: () => { pps = 12; redrawNow(); } }, "letters"), treesBtn,
      el("button", { class: "cardbtn", onclick: () => Card.open("sequence", g, li) }, "⎙ " + LM_STR.card.button));
    head.append(el("div", { class: "seq-line" }, ...parts), ctl);
    // the body: test type · sequence · Taxon name · p<α · the alignment; the two name columns resize by their header's edge
    const body = el("div", { class: "seq-body" + (T.pcol ? " with-p" : "") }); body.style.setProperty("--wrole", colW.role + "px"); body.style.setProperty("--wname", colW.name + "px");
    const treesPane = el("div", { class: "trees-pane", style: "display:none" });
    const mid = el("div", { class: "seq-mid" }, body, treesPane);
    const cols = { gut: el("div", { class: "seq-col seq-gut" }), role: el("div", { class: "seq-col seq-role" }), names: el("div", { class: "seq-col seq-names" }), p: el("div", { class: "seq-col seq-p" }) };
    const cwrap = el("div", { class: "seq-canvas-wrap" }); const cv = el("canvas"); cwrap.append(cv);
    body.append(cols.gut, cols.role, cols.names, ...(T.pcol ? [cols.p] : []), cwrap);
    const handle = key => { const h = el("i", { class: "hdl", title: "drag to resize" }); let drag = null;
      h.addEventListener("pointerdown", ev => { ev.preventDefault(); h.setPointerCapture(ev.pointerId); drag = { x: ev.clientX, w0: colW[key] }; h.classList.add("on"); });
      h.addEventListener("pointermove", ev => { if (!drag) return; colW[key] = clamp(drag.w0 + ev.clientX - drag.x, key === "role" ? 44 : 60, key === "role" ? 160 : 400); body.style.setProperty("--w" + key, colW[key] + "px"); redrawNow(); });
      h.addEventListener("pointerup", () => { drag = null; h.classList.remove("on"); try { localStorage.setItem("lm.seqcols", JSON.stringify(colW)); } catch (e) { /* no storage */ } });
      return h; };
    const trk = (col, txt, extra) => cols[col].append(el("div", { class: "trk", style: `height:${RULER_H + 4}px`, ...(extra || {}) }, txt));
    trk("gut", "test type"); cols.role.append(el("div", { class: "trk", style: `height:${RULER_H + 4}px` }, "sequence", handle("role"))); cols.names.append(el("div", { class: "trk", style: `height:${RULER_H + 4}px` }, "Taxon name", handle("name")));
    if (T.pcol) trk("p", "p<α", { title: `significant: * below alpha / N (${lmP(T.bar)}), ** below alpha / 5N, *** below alpha / 50N` });
    const toggleTrees = async () => {                                  // the trees dock beside the alignment, never over it
      if (treesPane.style.display !== "none") { treesPane.style.display = "none"; treesPane.replaceChildren(); treesBtn.classList.remove("primary"); showTrees = false; redrawNow(); return; }
      showTrees = true; treesBtn.classList.add("primary");
      treesPane.style.display = ""; treesPane.replaceChildren(el("div", { class: "loading", style: "padding: 10px" }, "loading the trees…")); redrawNow();
      const entry = await D.trees(li); if (treesPane.style.display === "none") return; treesPane.replaceChildren();
      let k = -1; if (entry) for (let q = 0; q < entry.cells.length; q++) if (entry.cells[q] === g) { k = q; break; }
      if (k < 0) { treesPane.append(el("div", { class: "note", style: "padding: 8px" }, "trees were not carried for this cell")); return; }
      treesPane.append(el("button", { class: "close", title: "close", onclick: toggleTrees }, "✕"), drawCophylo(treesPane, D, entry, k, { g, li }, T));
    };
    if (showTrees) setTimeout(toggleTrees, 0);
    // one row per sequence in every column; the focal's p cell is its outcome as a colour (gold: both tests) with the stars of
    // the decisive divergence p, the closest-match p, or the more significant of the two; a row's stars are its own test's
    rows.forEach((r, ri) => {
      if (ruleAt.has(ri)) for (const c of Object.values(cols)) c.append(el("div", { class: "rule", style: `height:${RULE_H}px` }));
      const note = testNote(r), nm = r.g === null ? (r.role === "broad" ? LM_STR.seq.broad_none : D.format >= 7 ? "none" : "not recorded") : D.taxa.names[r.g];
      const tip = `${r.g === null ? "" : nm + ": "}${r.tip}${r.role !== "focal" && r.shared ? `; differs from the focal at ${fmt1(100 * r.diff)}% of ${fmtInt(r.shared)} shared sites` : ""}${note ? "\n" + note : ""}`;
      cols.gut.append(el("div", { class: "seqnm gut", style: `height:${ROW_H}px` },
        r.role === "closest" ? LM_STR.seq.gut_closest : r.role === "broad" ? LM_STR.seq.gut_broad : ri === divAt ? LM_STR.seq.gut_relatives : ""));
      cols.role.append(el("div", { class: "seqnm role", style: `height:${ROW_H}px`, title: tip }, r.tag));
      cols.names.append(el("div", { class: "seqnm", style: `height:${ROW_H}px`, title: tip }, el("span", { class: "nm" + (r.st ? " st-" + r.st : "") + (r.g === null ? " none" : "") }, nm)));
      if (!T.pcol) return;
      let bg = D.col.passCell, txt = "–", ink = cssVar("--muted");
      if (r.role === "focal") { const mp = maskP(D, g, li); const fp = mp ? mp.p : dv && sm ? Math.min(divP ?? Infinity, simP ?? Infinity) : dv ? divP : sm ? simP : null; bg = dv && sm ? D.col.both : dv ? D.col.red : sm ? D.col.blue : D.col.passCell; txt = dv || sm ? T.stars(fp) : "–"; if (dv || sm) ink = "#ffffff"; }
      else if (r.role === "closest") { if (r.g !== null && T.blue.has(r.g)) { bg = D.col.blue; txt = T.stars(T.blue.get(r.g).p); ink = "#ffffff"; } }
      else if (r.role === "broad") { let vp = rcr0 && !Number.isNaN(rcr0.var_p) ? rcr0.var_p : null;
        if (vp === null && D.cols.bs_log10p && !Number.isNaN(D.cols.bs_log10p[i])) vp = D.cols.bs_log10p[i];
        const bsC = (s2v & 4) && !(s2v & 32), bsD = (s2v & 8) && !(s2v & 32);
        if (bsC || bsD) { bg = bsC ? D.col.blue : D.col.red; ink = "#ffffff"; }
        txt = vp !== null && T.stars(vp) ? T.stars(vp) : vp !== null ? "–" : "–"; if ((bsC || bsD) && !T.stars(vp)) txt = "*"; }
      else if (r.st === "red") { bg = D.col.red; txt = T.stars(T.red.get(r.g)); ink = "#ffffff"; }
      cols.p.append(el("div", { class: "pcell", style: `height:${ROW_H}px`, title: note || "" }, el("span", { style: `background:${bg}; color:${ink}` }, txt)));
    });
    const swc = c => el("span", { class: "sw", style: `background:${c}` });
    const key = el("span", { class: "seq-key", title: "a red site: of the focal's relatives present there (at least 3; its closest match is never counted), half or more differ from the focal; hover a site for the count" },
      "focal ", swc(D.col.redTint), " half or more of its relatives differ ", swc(D.col.blueTint), " and its closest match agrees · closest ", swc(D.col.blueTint), " agrees with the focal at such a site · relatives ", swc(D.col.redTint), " differ from the focal");
    const foot = el("div", { class: "seq-foot" }, key, el("span", null, `${fmtInt(len)} sites · pinch or ctrl+wheel: zoom · drag: move · letters from 7 px per site`), el("span", { class: "k", style: "margin-left:auto" }, S ? "" : "sequences were not carried for this locus"));
    P.body.replaceChildren(head, mid, foot);
    if (!S && !codeRow) { cwrap.append(el("div", { class: "note", style: "padding: 10px" }, "neither sequences nor site codes were carried for this cell in this report")); return; }
    const dpr = devicePixelRatio || 1;
    const fitPps = () => Math.max(0.05, cwrap.getBoundingClientRect().width / len);
    const topOf = ri => RULER_H + 4 + ri * ROW_H + RULE_H * rulesBefore(ri);
    const draw = () => {
      const W = Math.max(100, cwrap.getBoundingClientRect().width), H = topOf(rows.length) + 2;
      if (pps === null) { const f = fitPps(); pps = f >= 7 ? Math.min(f, 12) : f; }
      pps = clamp(pps, fitPps(), 16); x0 = clamp(x0, 0, Math.max(0, len - W / pps));
      cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + "px"; cv.style.height = H + "px";
      const c = cv.getContext("2d"); c.setTransform(dpr, 0, 0, dpr, 0, 0); c.fillStyle = cssVar("--surface"); c.fillRect(0, 0, W, H);
      const s0 = Math.max(0, Math.floor(x0)), s1 = Math.min(len - 1, Math.ceil(x0 + W / pps)); const letters = pps >= 7;
      const sx = s => (s - x0) * pps;
      // ruler
      c.fillStyle = cssVar("--muted"); c.font = `10px ${cssVar("--sans")}`; c.textBaseline = "top"; c.textAlign = "center";
      const stepSites = Math.max(1, Math.pow(10, Math.ceil(Math.log10(60 / pps))));
      for (let s = Math.ceil(s0 / stepSites) * stepSites; s <= s1; s += stepSites) { const x = sx(s) + pps / 2; c.fillRect(x, RULER_H - 4, 1, 4); c.fillText(fmtInt(alignmentPos(D, li, s) + 1), x, 2); }
      c.strokeStyle = cssVar("--grid"); c.beginPath(); c.moveTo(0, RULER_H + 0.5); c.lineTo(W, RULER_H + 0.5); c.stroke();
      // sequence rows
      const cellCol = D.col.cell, faint = cssVar("--faint"), redT = D.col.redTint, blueT = D.col.blueTint;
      c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; c.textBaseline = "middle"; c.textAlign = "center";
      c.strokeStyle = cssVar("--spine"); for (const ri of ruleAt) { c.beginPath(); c.moveTo(0, topOf(ri) - RULE_H / 2 - 0.5); c.lineTo(W, topOf(ri) - RULE_H / 2 - 0.5); c.stroke(); }
      rows.forEach((r, ri) => {
        const y = topOf(ri); const bytes = r.bytes; const isFocal = ri === 0, isCand = ri === 1;
        if (!bytes) { c.fillStyle = faint; c.font = `11px ${cssVar("--sans")}`; c.textAlign = "left"; c.fillText(r.g === null ? (D.format >= 7 ? "no closest match at this locus" : "closest match not recorded in this report") : "row not carried in this report", 4, y + ROW_H / 2); c.textAlign = "center"; c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; return; }
        if (letters) {
          // the focal is red at a red site and blue at a blue one; the closest match blue where it agrees with the focal at a red
          // site; a relative red where it differs from the focal; every other residue sits on the plain cell
          for (let s = s0; s <= s1; s++) {
            const b = bytes[s]; if (isGap(b)) continue; const ch = String.fromCharCode(b).toUpperCase(); const x = sx(s);
            const fp = !isFocal && focal && !isGap(focal[s]), same = fp && focal[s] === b;
            const bg = isFocal ? (siteBlue[s] ? blueT : siteRed[s] ? redT : cellCol) : isCand ? (same && siteRed[s] ? blueT : cellCol) : fp && !same ? redT : cellCol;
            c.fillStyle = bg; c.fillRect(x, y + 1, pps, ROW_H - 2);
            if (same && dots) { c.fillStyle = faint; c.fillText("·", x + pps / 2, y + ROW_H / 2); }
            else { c.fillStyle = inkFor(bg); c.fillText(ch, x + pps / 2, y + ROW_H / 2); }
          }
        } else {
          // block zoom, per pixel: the focal blue when any of its sites is blue, else red when any is red; the closest match blue
          // where it agrees at a red site; a relative as red as the share of its sites that differ from the focal
          for (let x = 0; x < W; x++) {
            const a = Math.floor(x0 + x / pps), b = Math.min(len - 1, Math.max(a, Math.ceil(x0 + (x + 1) / pps) - 1)); if (a > b) break;
            let pres = 0, diff = 0, red = false, blue = false;
            for (let s = a; s <= b; s++) { const bb = bytes[s]; if (isGap(bb)) continue; pres++;
              if (isFocal) { if (siteBlue[s]) blue = true; else if (siteRed[s]) red = true; continue; }
              if (!focal || isGap(focal[s])) continue;
              if (focal[s] !== bb) diff++; else if (isCand && siteRed[s]) blue = true; }
            if (!pres) continue;
            c.fillStyle = blue ? blueT : red ? redT : cellCol; c.fillRect(x, y + 1, 1, ROW_H - 2);
            if (!isFocal && !isCand && diff) { c.globalAlpha = 0.3 + 0.7 * diff / pres; c.fillStyle = redT; c.fillRect(x, y + 1, 1, ROW_H - 2); c.globalAlpha = 1; }
          }
        }
      });
    };
    drawFn = draw; drawNow = draw; draw();
    // zoom with a pinch (ctrl+wheel) or two-pointer spread; scroll with the wheel or a drag
    cv.addEventListener("wheel", ev => { ev.preventDefault(); const r = cv.getBoundingClientRect(); const mx = ev.clientX - r.left;
      if (ev.ctrlKey || ev.metaKey) { const f = Math.exp(-ev.deltaY * 0.01); const np = clamp(pps * f, fitPps(), 16); x0 = x0 + mx / pps - mx / np; pps = np; }
      else { const d = Math.abs(ev.deltaX) > Math.abs(ev.deltaY) ? ev.deltaX : ev.deltaY; x0 += d / pps; }
      draw(); }, { passive: false });
    const ptrs = new Map(); let drag = null, pinch = null;
    cv.addEventListener("pointerdown", ev => { cv.setPointerCapture(ev.pointerId); ptrs.set(ev.pointerId, ev.clientX); if (ptrs.size === 2) { const [a, b] = [...ptrs.values()]; pinch = { d: Math.abs(a - b), pps, x0, mid: (a + b) / 2 - cv.getBoundingClientRect().left }; drag = null; } else drag = { x: ev.clientX, x0 }; });
    cv.addEventListener("pointerup", ev => { ptrs.delete(ev.pointerId); if (ptrs.size < 2) pinch = null; drag = null; });
    cv.addEventListener("pointercancel", ev => { ptrs.delete(ev.pointerId); pinch = null; drag = null; });
    cv.addEventListener("pointermove", ev => {
      if (ptrs.has(ev.pointerId)) ptrs.set(ev.pointerId, ev.clientX);
      if (pinch && ptrs.size === 2) { const [a, b] = [...ptrs.values()]; const d = Math.abs(a - b); if (pinch.d > 0) { const np = clamp(pinch.pps * d / pinch.d, fitPps(), 16); x0 = pinch.x0 + pinch.mid / pinch.pps - pinch.mid / np; pps = np; draw(); } return; }
      if (drag) { x0 = drag.x0 - (ev.clientX - drag.x) / pps; draw(); return; }
      const r = cv.getBoundingClientRect(); const s = Math.floor(x0 + (ev.clientX - r.left) / pps); if (s < 0 || s >= len) { hideHover(); return; }
      const t = el("table"); t.append(el("tr", null, el("td", { class: "k" }, "site"), el("td", null, `${s + 1} of ${fmtInt(len)} · alignment position ${fmtInt(alignmentPos(D, li, s) + 1)}`)));
      if (dRow[s] >= 0) t.append(el("tr", null, el("td", { class: "k" }, "difference"), el("td", null, `${dRow[s]} of ${nRow[s]} relatives differ from the focal${siteRed[s] ? " (half or more)" : nRow[s] < 3 ? " (fewer than 3 present)" : ""}${siteBlue[s] ? `; its closest match${candRow ? ` ${D.taxa.names[candRow.g]}` : ""} agrees with it` : ""}`)));
      for (const rr of rows) if (rr.bytes) t.append(el("tr", null, el("td", { class: "k" }, `${rr.tag} · ${D.taxa.names[rr.g]}`), el("td", { style: "font-family: var(--mono)" }, String.fromCharCode(rr.bytes[s]))));
      showHover(t, ev.clientX, ev.clientY);
    });
    cv.addEventListener("pointerleave", hideHover);
    App.panels.sequence.redraw = draw;
  }
  return { open, render, cur: () => cur };
})();

// ------------------------------------------------------------------ histograms: loci or taxa by masked count, one facet sized to the panel
const Hist = (() => {
  function render(root, by, width, height) { root.replaceChildren(); root.append(facet(by, root, width, height)); }
  function bins(kof) {
    const pos = kof.filter(k => k > 0); const kmax = d3.max(kof) || 0; let w = 1;
    if (pos.length >= 2) { const q = pos.slice().sort(d3.ascending); const iqr = d3.quantileSorted(q, 0.75) - d3.quantileSorted(q, 0.25); w = 2 * iqr / Math.cbrt(pos.length); if (w < 1) w = kmax / Math.max(1, 1 + Math.log2(pos.length)); w = Math.max(1, Math.ceil(w)); }
    const out = []; for (let lo = 1; lo <= kmax; lo += w) out.push([lo, Math.min(lo + w - 1, kmax)]); return out;   // units holding a mask; the zero bin is the subtitle's job
  }
  function facet(by, root, width, height) {
    const D = App.D; const N = by === "locus" ? D.nloc : D.n; const cnt = by === "locus" ? D.locCount : D.taxCount; const names = by === "locus" ? D.loci.names : D.taxa.names;
    const units = d3.range(N).filter(u => cnt.present[u] > 0 && cnt.masked[u] > 0); const kof = units.map(u => cnt.masked[u]);
    const present = d3.sum(d3.range(N), u => cnt.present[u] > 0 ? 1 : 0);
    const B = bins(kof); const heights = B.map(([lo, hi]) => units.filter((u, j) => kof[j] >= lo && kof[j] <= hi));
    const rotated = B.length > 14; const W = Math.max(300, width), H = Math.max(170, height), m = { l: 44, r: 12, t: 22, b: rotated ? 76 : 40 };
    const xs = d3.scaleBand().domain(B.map((b, i) => i)).range([m.l, W - m.r]).padding(0.15);
    const ymax = d3.max(heights, h => h.length) || 1; const ys = d3.scaleLinear().domain([0, ymax * 1.15]).range([H - m.b, m.t]);
    const box = el("div", { class: "plot" }, el("h3", null, by === "locus" ? LM_STR.overview.hist_loci : LM_STR.overview.hist_taxa, el("span", { class: "sub" }, `${fmtInt(units.length)} of ${fmtInt(present)}`)));
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).ticks(5).tickFormat(d3.format(",")));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).tickFormat(i => { const [lo, hi] = B[i]; return lo === hi ? String(lo) : `${lo}–${hi}`; })).selectAll("text").attr("transform", rotated ? "rotate(-60)" : null).attr("text-anchor", rotated ? "end" : "middle").attr("dx", rotated ? "-0.5em" : null).attr("dy", rotated ? "0.1em" : null);
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 4).attr("text-anchor", "middle").text(lmFmt(LM_STR.overview.hist_x, { unit: by }));
    svg.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text(by === "locus" ? "loci" : "taxa");
    B.forEach(([lo, hi], i) => {
      const us = heights[i]; if (!us.length) return;
      const comp = new Map();
      { for (const u of us) for (let k = 0; k < (by === "locus" ? D.n : D.nloc); k++) { const idx = by === "locus" ? k * D.nloc + u : u * D.nloc + k; const code = D.state[idx]; if (code >= ST.RED) { const f = D.famOf(code); comp.set(f, (comp.get(f) || 0) + 1); } } }   // one segment per outcome
      const tot = d3.sum(comp.values()) || 1; let y0 = ys(0); const codes = Array.from(comp.keys()).sort(d3.ascending);
      for (const code of codes) { const h = (ys(0) - ys(us.length)) * comp.get(code) / tot; svg.append("rect").attr("x", xs(i)).attr("y", y0 - h).attr("width", xs.bandwidth()).attr("height", h).attr("fill", D.colorOf(code)).attr("stroke", cssVar("--surface")).attr("stroke-width", 0.8).style("cursor", "pointer")
        .on("mousemove", ev => showHover(el("div", null, el("b", null, `${fmtInt(us.length)} ${by === "locus" ? "loci" : "taxa"} with ${lo === hi ? lo : `${lo} to ${hi}`} masked sequence${hi !== 1 ? "s" : ""}`), lo === 0 ? null : el("div", null, `${D.famLabel(code)}: ${fmtInt(comp.get(code))} sequences (${pct(comp.get(code), tot)})`), el("div", { class: "more" }, us.slice(0, 8).map(u => names[u]).join(", ") + (us.length > 8 ? ", …" : "") + " · click to highlight them in the matrix")), ev.clientX, ev.clientY))
        .on("mouseleave", hideHover).on("click", () => { if (by === "locus") { App.sel.cols = new Set(us); App.sel.rows.clear(); } else { App.sel.rows = new Set(us); App.sel.cols.clear(); } Overview.setView("matrix"); }); y0 -= h; }
      svg.append("text").attr("x", xs(i) + xs.bandwidth() / 2).attr("y", ys(us.length) - 3).attr("text-anchor", "middle").style("fill", cssVar("--ink")).text(fmtInt(us.length));
    });
    box.append(svg.node()); svg.append("title").text("zero is its own bar; the other bars are coloured by the kinds of mask they hold. The y axis is linear near zero and logarithmic above. Click a bar to highlight its units in the matrix.");
    return box;
  }
  return { render };
})();

// ------------------------------------------------------------------ the cards: one printable page per dataset,
// locus, taxon or sequence. One design language: identity strip, the verdict band, the three per-test ribbons on
// dataset-wide axes (the spine every card shares, so pages compare), the problems table, one true-axed shape plot,
// a traceable footer. The preview overlays the app; print CSS shows the card alone
const Card = (() => {
  let host = null;
  const S = () => LM_STR.card;
  function ensure() {
    if (host) return;
    host = el("div", { id: "lm-card", style: "display:none" });
    document.body.append(host);
  }
  function close() { host.style.display = "none"; host.replaceChildren(); document.body.classList.remove("card-open"); }
  function worstOf(masks) { let w = null; for (const q of masks) if (Number.isFinite(q.p) && (w === null || q.p < w)) w = q.p; return w; }
  function famChip(D, masks) {
    let best = 0; for (const q of masks) if (q.code > best) best = q.code;
    return best ? sw(D.colorOf(D.famOf(best) === ST.BOTH ? ST.BOTH : best)) : sw(null);
  }
  function verdictBand(D, masks, n, m, stateLine) {
    const line = stateLine !== null && stateLine !== undefined ? stateLine
      : masks.length === 0 ? S().verdict_clean : lmFmt(S().verdict_unit, { n: fmtInt(masks.length), m: fmtInt(m), pct: pct(masks.length, m) });
    return el("div", { class: "cd-verdict" + (masks.length ? "" : " clean") }, famChip(D, masks),
      el("span", { class: "cd-line" }, line));
  }
  function spineCounts(D, kind, u) {                          // per test arm: how many masks it made WITHIN this card's scope
    const key = kind + ":" + u; D._spineCounts = D._spineCounts || {};
    if (D._spineCounts[key]) return D._spineCounts[key];
    const c = { close_relative: 0, cm_relative: 0, cm_outside: 0, broad_sample: 0 };
    for (const q of maskIndex(D)) {
      if (kind === "locus" && q.li !== u) continue;
      if (kind === "taxon" && q.g !== u) continue;
      const mk = maskKind(D, q.g, q.li); if (!mk) continue;
      if (mk.cr) c.close_relative++;
      if (mk.cm) c[mk.arm === "out" ? "cm_outside" : "cm_relative"]++;
      if (mk.bsC || mk.bsD) c.broad_sample++; }
    return (D._spineCounts[key] = c);
  }
  function spine(D, kind, u, cell) {
    const box = el("div", { class: "cd-spine" }); const W = 700, RH = 32, PAD = 24;
    const counts = kind === "sequence" ? null : spineCounts(D, kind, u);
    for (const mode of ["close_relative", "cm_relative", "cm_outside", "broad_sample"]) {
      if (mode === "broad_sample" && !D.cols.bs_stat) continue;       // old bundles have no variance statistic to place
      const row = el("div", { class: "cd-rib" },
        el("span", { class: "cd-ribname" }, LM_STR.ribbons["test_" + mode], counts ? el("i", { class: "cd-ribn" }, fmtInt(counts[mode])) : null));
      const cv = el("canvas"); row.append(cv); box.append(row);
      if (kind === "dataset") drawResidualAll(cv, W, RH, mode);
      else drawResidualRow(cv, kind === "sequence" ? "locus" : kind, kind === "sequence" ? cell[1] : u, W, RH, { padL: PAD, padR: PAD, mode, noBand: true, mark: kind === "sequence" ? cell : null });
      box.append(el("div", { class: "cd-axis" }, residualAxis(W, PAD, PAD, mode)));
    }
    return box;
  }
  function drawResidualAll(cv, W, H, mode) {                   // the dataset card's strip: every cell once, the ribbon's OWN masks coloured
    const D = App.D; const dpr = devicePixelRatio || 1; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + "px"; cv.style.height = H + "px";
    const c = cv.getContext("2d"); c.scale(dpr, dpr);
    const xs = residualScale(D, mode, W, 20, 20); const cy = H / 2; const jit = d3.randomUniform.source(d3.randomLcg(0.41))(-(H / 2 - 4), H / 2 - 4);
    c.strokeStyle = cssVar("--model"); c.beginPath(); c.moveTo(xs(0) + 0.5, 1); c.lineTo(xs(0) + 0.5, H - 1); c.stroke();
    const N = D.n * D.nloc; const step = Math.max(1, Math.floor(N / 30000));
    c.fillStyle = cssVar("--pass"); c.globalAlpha = 0.25;
    for (let i = 0; i < N; i += step) { if (D.cols.n_comp[i] === 0 || (D.cols.side[i] & 3)) continue;
      const v = residualOf(D, Math.floor(i / D.nloc), i % D.nloc, mode); if (v === null) continue;
      c.beginPath(); c.arc(xs(v), cy + jit(), 1.5, 0, 2 * Math.PI); c.fill(); }
    c.globalAlpha = 1;
    for (const q of maskIndex(D)) { const v = residualOf(D, q.g, q.li, mode); if (v === null) continue;
      const col = ribbonColour(D, q.g, q.li, mode); if (!col) continue;        // not this test's mask: it stays in the grey cloud
      c.fillStyle = col;
      c.beginPath(); c.arc(xs(v), cy + jit(), 2.6, 0, 2 * Math.PI); c.fill(); c.strokeStyle = cssVar("--surface"); c.lineWidth = 0.6; c.stroke(); }
  }
  function maskLollipop(D, masks, kind, W, N) {               // who is involved, by significance: the card grammar's ranked form
    const rows = masks.filter(q => Number.isFinite(q.p)).sort((a, b) => a.p - b.p).slice(0, N);
    if (!rows.length) return null;
    const RH = 21, m = { l: 10, r: 46, t: 4, b: 36 }, nameW = 170; const H = m.t + rows.length * RH + m.b;
    const vmax = Math.max(6, -rows[0].p * 1.08);
    const xs = d3.scaleLinear().domain([0, vmax]).range([m.l + nameW, W - m.r]);
    const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "cd-rank");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).ticks(3).tickFormat(v => lmP(-v)));
    rows.forEach((q, i) => {
      const y = m.t + i * RH + RH / 2;
      const nm = kind === "locus" ? D.taxa.names[q.g] : D.loci.names[q.li];
      svg.append("text").attr("x", m.l + nameW - 8).attr("y", y + 3.5).attr("text-anchor", "end").attr("font-size", 11).attr("fill", cssVar("--ink")).text(nm);
      svg.append("line").attr("x1", xs(0)).attr("x2", xs(-q.p)).attr("y1", y).attr("y2", y).attr("stroke", cssVar("--spine")).attr("stroke-width", 1.4);
      const col = D.colorOf(q.code); svg.append("circle").attr("cx", xs(-q.p)).attr("cy", y).attr("r", 4).attr("fill", col === "hatch" || !col ? cssVar("--pass") : col);
      svg.append("text").attr("x", xs(-q.p) + 8).attr("y", y + 3.5).attr("font-size", 10.5).attr("fill", cssVar("--muted")).text(starsOf(D, q.p));
    });
    svg.append("text").attr("x", (m.l + nameW + W - m.r) / 2).attr("y", H - 3).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", cssVar("--faint")).text(S().mask_x);
    return el("div", { class: "cd-plotrow half" }, el("h4", null, lmFmt(kind === "locus" ? S().masked_taxa : S().masked_loci, { n: fmtInt(rows.length) })), svg.node());
  }
  function scopedFlowPlot(D, entries, groupBy, title, W, N) { // donors into / recipients out of this card's subject
    const nm = D.taxa.names;
    const groups = d3.groups(entries, e => String(groupBy(e))).map(([k, v]) => ({ label: nm[+k], n: v.length }))
      .sort((a, b) => b.n - a.n || a.label.localeCompare(b.label)).slice(0, N);
    if (!groups.length) return null;
    const total = entries.length || 1; const RH = 21, m = { l: 10, r: 46, t: 4, b: 36 }, nameW = 170; const H = m.t + groups.length * RH + m.b;
    const vmax = Math.max(5, 100 * groups[0].n / total * 1.12);
    const xs = d3.scaleLinear().domain([0, vmax]).range([m.l + nameW, W - m.r]);
    const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "cd-rank");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).ticks(4).tickFormat(v => v + "%"));
    groups.forEach((gr, i) => {
      const y = m.t + i * RH + RH / 2; const v = 100 * gr.n / total;
      svg.append("text").attr("x", m.l + nameW - 8).attr("y", y + 3.5).attr("text-anchor", "end").attr("font-size", 11).attr("fill", cssVar("--ink")).text(gr.label);
      svg.append("line").attr("x1", xs(0)).attr("x2", xs(v)).attr("y1", y).attr("y2", y).attr("stroke", cssVar("--spine")).attr("stroke-width", 1.4);
      svg.append("circle").attr("cx", xs(v)).attr("cy", y).attr("r", 4).attr("fill", D.col.blue);
      svg.append("text").attr("x", xs(v) + 8).attr("y", y + 3.5).attr("font-size", 10.5).attr("fill", cssVar("--muted")).text(fmtInt(gr.n));
    });
    svg.append("text").attr("x", (m.l + nameW + W - m.r) / 2).attr("y", H - 3).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", cssVar("--faint")).text(S().flow_x);
    return el("div", { class: "cd-plotrow half" }, el("h4", null, lmFmt(title, { n: fmtInt(groups.length) })), svg.node());
  }
  function rankStrip(D, kind, W, H) {                         // every unit's % masked, ranked: the shape of the damage, one row, a little title
    const cnt = kind === "taxon" ? D.taxCount : D.locCount; const names = kind === "taxon" ? D.taxa.names : D.loci.names;
    const units = [];
    for (let u = 0; u < (kind === "taxon" ? D.n : D.nloc); u++) if (cnt.tested[u] > 0) units.push([u, 100 * cnt.masked[u] / cnt.tested[u]]);
    units.sort((a, b) => b[1] - a[1]);
    if (units.length > 400) units.length = 400;                                     // thousands of near-zero units say nothing
    const m = { l: 44, r: 12, t: 6, b: 32 };
    const p99 = units.length ? units[Math.floor(units.length * 0.01)][1] : 0;      // one extreme unit must not flatten the curve:
    const ymax = Math.max(2, Math.min((units[0] ? units[0][1] : 0) * 1.08, Math.max(p99 * 2.5, 5)));   // the label carries the top
    const xs = d3.scaleLinear().domain([0, units.length]).range([m.l, W - m.r]), ys = d3.scaleLinear().domain([0, ymax]).range([H - m.b, m.t]);
    const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "cd-rank");
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).ticks(4).tickFormat(v => v + "%"));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).ticks(6).tickFormat(d3.format("~s")));
    const area = d3.area().x((d, i) => xs(i)).y0(ys(0)).y1(d => ys(Math.min(d[1], ymax)));
    const step = Math.max(1, Math.floor(units.length / (W - m.l - m.r)));
    const pts = []; for (let i = 0; i < units.length; i += step) pts.push(units[i]);
    svg.append("path").attr("d", area(pts)).attr("fill", cssVar("--pass")).attr("opacity", 0.75);
    if (units.length && units[0][1] > 0) {
      svg.append("text").attr("x", xs(0) + 4).attr("y", Math.max(m.t + 9, ys(Math.min(units[0][1], ymax)) - 4)).attr("font-size", 10).attr("fill", cssVar("--ink"))
        .text(`${names[units[0][0]]} ${fmt1(units[0][1])}%`);
    }
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 3).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", cssVar("--faint")).text(S().rank_x);
    return el("div", { class: "cd-plotrow" }, el("h4", null, kind === "taxon" ? S().rank_taxa : S().rank_loci), svg.node());
  }
  function flowPlot(D, mode, W, N) {                          // the canonical ranked-category form: a lollipop per taxon, in % of attributed copies
    const groups = copyGroups(D, mode); const total = d3.sum(groups, gr => gr.list.length) || 1;
    const top = groups.slice(0, N); const RH = 21;
    const m = { l: 10, r: 46, t: 4, b: 36 }; const H = m.t + top.length * RH + m.b;
    const nameW = 150;
    const vmax = Math.max(5, 100 * (top[0] ? top[0].list.length : 0) / total * 1.12);
    const xs = d3.scaleLinear().domain([0, vmax]).range([m.l + nameW, W - m.r]);
    const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "cd-rank");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).ticks(4).tickFormat(v => v + "%"));
    top.forEach((gr, i) => {
      const y = m.t + i * RH + RH / 2; const v = 100 * gr.list.length / total;
      svg.append("text").attr("x", m.l + nameW - 8).attr("y", y + 3.5).attr("text-anchor", "end").attr("font-size", 11).attr("fill", cssVar("--ink")).text(gr.label);
      svg.append("line").attr("x1", xs(0)).attr("x2", xs(v)).attr("y1", y).attr("y2", y).attr("stroke", cssVar("--spine")).attr("stroke-width", 1.4);
      svg.append("circle").attr("cx", xs(v)).attr("cy", y).attr("r", 4).attr("fill", D.col.blue);
      svg.append("text").attr("x", xs(v) + 8).attr("y", y + 3.5).attr("font-size", 10.5).attr("fill", cssVar("--muted")).attr("font-variant", "tabular-nums").text(fmtInt(gr.list.length));
    });
    svg.append("text").attr("x", (m.l + nameW + W - m.r) / 2).attr("y", H - 3).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", cssVar("--faint")).text(S().flow_x);
    return el("div", { class: "cd-plotrow half" }, el("h4", null, lmFmt(mode === "donor" ? S().donors : S().recipients, { n: fmtInt(top.length) })), svg.node());
  }
  function drawRibbonInto(cv, C, g, li) {                     // the focal's sites along the locus: red where relatives differ,
    const D = App.D;                                           // blue where the partner agrees at such a site
    let r = -1; for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) { r = k; break; }
    if (r < 0) { cv.remove(); return; }
    const len = C.len; const W = Math.max(300, Math.floor(cv.getBoundingClientRect().width || 700));
    cv.width = W; cv.height = 1; cv.style.height = "12px"; cv.style.imageRendering = "pixelated";
    const img = cv.getContext("2d").createImageData(W, 1); const px = img.data; const base = r * len;
    for (let x = 0; x < W; x++) {
      const s0 = Math.floor(x * len / W), s1 = Math.max(s0, Math.ceil((x + 1) * len / W) - 1);
      let best = -1, code = 255;
      for (let sdx = s0; sdx <= Math.min(s1, len - 1); sdx++) { const c = C.blk[base + sdx]; const sv = codeSeverity(D, c); if (sv > best) { best = sv; code = c; } }
      const rgb = codeRgb(D, code);
      px[x * 4] = rgb[0]; px[x * 4 + 1] = rgb[1]; px[x * 4 + 2] = rgb[2]; px[x * 4 + 3] = 255;
    }
    cv.getContext("2d").putImageData(img, 0, 0);
  }
  function open(kind, a, b) {
    ensure(); const D = App.D;
    const masksAll = maskIndex(D);
    let title, masks, cells, m, spineEl, shape, extra = null, seqState = null;
    if (kind === "dataset") {
      title = D.man.dataset; masks = masksAll; m = D.man.counts.cells;
      spineEl = spine(D, "dataset");
      shape = el("div", null, rankStrip(D, "taxon", 820, 130), rankStrip(D, "locus", 820, 130));
      if (copiesOf(D).length) extra = el("div", { class: "cd-row" }, flowPlot(D, "donor", 400, 10), flowPlot(D, "recipient", 400, 10));
    } else if (kind === "locus" || kind === "taxon") {
      const u = a; title = kind === "locus" ? D.loci.names[u] : D.taxa.names[u];
      masks = masksAll.filter(q => (kind === "locus" ? q.li : q.g) === u);
      cells = []; const N = kind === "locus" ? D.n : D.nloc;
      for (let k = 0; k < N; k++) { const g = kind === "locus" ? k : u, li = kind === "locus" ? u : k; if (D.cols.n_comp[g * D.nloc + li] > 0) cells.push([g, li]); }
      m = cells.length;
      spineEl = spine(D, kind, u);
      const plots = [maskLollipop(D, masks, kind, 400, 10)];
      const cps = copiesOf(D);
      if (kind === "locus") {
        const here = cps.filter(e => e.li === u);
        plots.push(scopedFlowPlot(D, here, e => e.don, S().donors, 400, 10), scopedFlowPlot(D, here, e => e.rec, S().recipients, 400, 10));
      } else {
        plots.push(scopedFlowPlot(D, cps.filter(e => e.rec === u), e => e.don, S().received, 400, 10),
                   scopedFlowPlot(D, cps.filter(e => e.don === u), e => e.rec, S().gave, 400, 10));
      }
      shape = el("div", { class: "cd-row" }, ...plots.filter(Boolean));
    } else {
      const g = a, li = b; title = `${D.taxa.names[g]} · ${D.loci.names[li]}`;
      masks = masksAll.filter(q => q.g === g && q.li === li); m = 1; seqState = D.cellLabel(D.state[g * D.nloc + li], g * D.nloc + li);
      cells = []; for (let gg = 0; gg < D.n; gg++) if (D.cols.n_comp[gg * D.nloc + li] > 0) cells.push([gg, li]);
      spineEl = spine(D, "sequence", null, [g, li]); shape = null;
      const alignSlot = el("div", { class: "cd-plotrow cd-align" }), treesSlot = el("div", { class: "cd-plotrow cd-trees" });
      extra = el("div", null, el("div", { class: "cd-evidence" }, cellHover(g, li, false)), alignSlot, treesSlot);
      (async () => {                                                   // the alignment strip and the trees arrive from the bundle
        const [C, SIM, entry] = await Promise.all([D.codes(li), D.similar(li), D.trees(li)]);
        if (!host || host.style.display === "none") return;
        if (C) {
          alignSlot.append(el("h4", null, S().align_title));
          const cv = el("canvas", { style: "width: 100%" }); alignSlot.append(cv);
          requestAnimationFrame(() => drawRibbonInto(cv, C, g, li));
        }
        let k = -1; if (entry) for (let q = 0; q < entry.cells.length; q++) if (entry.cells[q] === g) { k = q; break; }
        if (k >= 0) {
          treesSlot.append(el("h4", null, S().trees_title));
          const T = cellTests(App.D, g, li, SIM, similarRole(App.D, g, li));
          treesSlot.append(drawCophylo(treesSlot, App.D, entry, k, { g, li }, T));
        }
      })();
    }
    const page = el("div", { class: "cd-page" },
      el("div", { class: "cd-id" }, el("span", { class: "cd-type" }, S()["type_" + kind]), el("b", null, title)),
      verdictBand(D, masks, masks.length, m, seqState),
      kind === "dataset" ? null : extra, spineEl, kind === "dataset" ? extra : null,
      shape,
      el("div", { class: "cd-foot" }, lmFmt(S().footer, { dataset: D.man.dataset, version: D.man.version, date: (D.man.date || "").slice(0, 10), bar: lmP(D.man.test.log10_bar) })));
    const bar = el("div", { class: "cd-bar" }, el("button", { onclick: () => window.print() }, S().print), el("button", { onclick: close }, S().close));
    host.replaceChildren(bar, page); host.style.display = ""; document.body.classList.add("card-open");
  }
  return { open, close };
})();

// ------------------------------------------------------------------ the table view (virtualised)
const Table = (() => {
  let rows = [], keys = [{ k: "p", dir: 1 }], maskedOnly = true, family = "all", root = null;   // keys[0] sorts first, keys[1] breaks its ties
  const T = () => LM_STR.table;
  const COLS = () => [["taxon", T().h_taxon], ["locus", T().h_locus], ["test", T().h_test], ["p", T().h_p], ["obs", T().h_obs], ["exp", T().h_exp], ["excess_m", T().h_excess], ["partner", T().h_partner], ["bm", T().h_best]];
  function build() {
    const D = App.D; rows = [];
    for (let g = 0; g < D.n; g++) for (let li = 0; li < D.nloc; li++) { const i = g * D.nloc + li; if (D.cols.n_comp[i] === 0) continue; const masked = D.cols.side[i] & 3;
      if (maskedOnly && !masked) continue;
      const mp = masked ? maskP(D, g, li) : null; const test = mp ? mp.test.replace(/-/g, "_") : "";
      if (family !== "all" && test !== family) continue;
      const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li); const bm = D.num("best_match", g, li), pt = D.num("partner", g, li);
      rows.push({ g, li, taxon: D.taxa.names[g], locus: D.loci.names[li], test: mp ? mp.test : "", p: mp ? mp.p : null,
        n_fail: D.cols.n_fail[i], n_similar: D.cols.n_similar[i], obs: mi === null ? null : 100 - mi, exp: me === null ? null : 100 - me,
        excess_m: D.num("excess_m", g, li), bm_id: D.num("best_match_identity", g, li), lam: D.num("lam", g, li),
        partner: pt !== null && pt >= 0 ? D.taxa.names[pt] : "", bm: bm !== null && bm >= 0 ? D.taxa.names[bm] : "", code: D.state[i] }); }
    sort();
  }
  function open(fam) { family = fam || "all"; maskedOnly = true; Overview.setView("table"); }
  const cmp = (a, b, k, dir) => { const x = a[k], y = b[k]; const ex = x === null || x === undefined || x === "", ey = y === null || y === undefined || y === ""; if (ex && ey) return 0; if (ex) return 1; if (ey) return -1; return (typeof x === "string" ? x.localeCompare(y) : x - y) * dir; };
  function sort() { rows.sort((a, b) => { for (const { k, dir } of keys) { const c = cmp(a, b, k, dir); if (c) return c; } return 0; }); }
  function render(r) {
    root = r; root.replaceChildren(); build(); const D = App.D;
    const cols = COLS().map(([k, t]) => [k, t, null]);
    const mark = k => { const j = keys.findIndex(x => x.k === k); return j < 0 ? "" : (keys[j].dir < 0 ? " ▾" : " ▴") + (keys.length > 1 ? (j ? "²" : "¹") : ""); };
    const click = (k, ev) => { const j = keys.findIndex(x => x.k === k); const natural = typeof rows.find(r => r[k] !== null && r[k] !== undefined && r[k] !== "")?.[k] === "string" ? 1 : -1;
      if (ev.shiftKey) { if (j === 0 && keys.length === 1) return; if (j === 1) keys[1].dir = -keys[1].dir; else keys = [keys[0], { k, dir: natural }]; }     // shift-click: the secondary key
      else if (j === 0) keys[0].dir = -keys[0].dir; else keys = [{ k, dir: natural }, ...(keys[0].k === k ? [] : [keys[0]])].slice(0, 2);                 // click: the primary key; the old primary breaks ties
      sort(); paint(); head.replaceChildren(...cols.map(([k2, t2, h2]) => el("span", { onclick: e => click(k2, e), title: h2 }, t2 + mark(k2)))); };
    const head = el("div", { class: "thead" }, ...cols.map(([k, t, h]) => el("span", { onclick: e => click(k, e), title: h }, t + mark(k))));
    const body = el("div", { class: "tbody" }); const inner = el("div"); body.append(inner);
    const chips = el("span", { class: "chips" }, ...[["all", T().chips.all], ["close_relative", T().chips.close_relative], ["closest_match", T().chips.closest_match], ["broad_sample", T().chips.broad_sample]].map(([k, t]) =>
      el("button", { class: "chip" + (family === k ? " on" : ""), "data-fam": k, onclick: () => { family = k; render(root); } }, t)));
    const foot = el("div", { class: "tfoot" }, chips, el("label", null, el("input", { type: "checkbox", checked: maskedOnly ? "" : null, onchange: e => { maskedOnly = e.target.checked; build(); paint(); } }), " " + T().masked_only), el("span", { class: "n" }));
    root.append(head, body, foot);
    const RH = 22;
    const paint = () => { inner.style.height = rows.length * RH + "px"; inner.style.position = "relative"; inner.replaceChildren(); const top = body.scrollTop; const i0 = Math.max(0, Math.floor(top / RH) - 5), i1 = Math.min(rows.length, i0 + Math.ceil(body.clientHeight / RH) + 10);
      for (let i = i0; i < i1; i++) { const r = rows[i]; const st = starsOf(D, r.p); inner.append(el("div", { class: "trow", style: `position:absolute; top:${i * RH}px; left:0; right:0`, onclick: () => Sequence.open(r.g, r.li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(r.g * D.nloc + r.li); } },
        el("span", null, sw(D.colorOf(r.code)), r.taxon), el("span", null, r.locus), el("span", null, r.test), el("span", { class: "pcol" }, r.p === null ? "" : lmP(r.p), st ? el("i", { class: "stars" }, " " + st) : null),
        el("span", null, r.obs === null ? "" : fmt2(r.obs)), el("span", null, r.exp === null ? "" : fmt2(r.exp)), el("span", null, r.excess_m === null ? "" : fmt2(r.excess_m)), el("span", null, r.partner), el("span", null, r.bm))); }
      $(".n", foot).textContent = lmFmt(T().foot, { n: fmtInt(rows.length) }); };
    body.addEventListener("scroll", paint); paint();
  }
  function download() {
    const D = App.D; const keepM = maskedOnly, keepF = family; maskedOnly = false; family = "all"; build(); maskedOnly = keepM; family = keepF;
    const v4 = D.format >= 4;                                 // from format 4: the cell's role in its reported pair and the tests' log10 p, in the tests' names
    const p6 = v => v === null || v === undefined || v === "" || Number.isNaN(v) ? "" : Math.min(Math.pow(10, v), 1).toFixed(6);
    const lines = ["taxon\tlocus\tstate\tmask_test\tmask_p\tn_fail\tn_similar\tobserved_divergence\texpected_divergence\texcess_m\tbest_match\tbest_match_identity\tpartner\tdirection_lrt" + (v4 ? "\tclosest_match_role\tclose_relative_p\tclosest_match_p\tbroad_sample_p" : "")];
    for (const r of rows) { const pt = D.num("partner", r.g, r.li); const bm = D.num("best_match", r.g, r.li); const sr = v4 ? similarRole(D, r.g, r.li) : null;
      const i = r.g * D.nloc + r.li; const rcr = D.rcAt ? D.rcAt(i) : null;
      const bsp = rcr && !Number.isNaN(rcr.var_p) ? rcr.var_p : (D.cols.bs_log10p && !Number.isNaN(D.cols.bs_log10p[i]) ? D.cols.bs_log10p[i] : null);
      lines.push([r.taxon, r.locus, D.labelOf(r.code), r.test, p6(r.p), r.n_fail, r.n_similar, r.obs ?? "", r.exp ?? "", r.excess_m ?? "", bm !== null && bm >= 0 ? D.taxa.names[bm] : "", r.bm_id ?? "", pt !== null && pt >= 0 ? D.taxa.names[pt] : "", r.lam ?? "", ...(v4 ? [sr ? sr.role : "", p6(D.num("div_log10p", r.g, r.li)), p6(D.num("sim_log10p", r.g, r.li)), p6(bsp)] : [])].join("\t")); }
    const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/tab-separated-values" })); a.download = `${D.man.dataset || "report"}_cells.tsv`; a.click(); build();
  }
  return { render, download, open };
})();

// ------------------------------------------------------------------ URL hash: view, ordering, viewport, open panels
let hashTimer = null;
function writeHash() {
  clearTimeout(hashTimer); hashTimer = setTimeout(() => {
    const D = App.D; if (!D) return; const V = Overview.V; const p = new URLSearchParams();
    p.set("d", D.man.input.sha256.slice(0, 10)); p.set("view", App.view); p.set("o", `${App.order.taxaMode},${App.order.lociMode}`);
    p.set("v", [V.x0, V.y0, V.cw, V.ch].map(v => +v.toFixed(3)).join(","));
    const dc = Detail.cur(); if (dc && App.panels.detail) { if (dc.kind === "copies") p.set("cp", dc.mode + (dc.key !== null ? ":" + dc.key : "")); else { p.set(dc.kind === "taxon" ? "t" : "l", dc.i); const ds = App.panels.detail.P.el.querySelector(".crumb select"); if (ds) p.set("dv", ds.value); } }
    const sc = Sequence.cur(); if (sc && App.panels.sequence) p.set("c", `${sc.g},${sc.li}`);
    try { history.replaceState(null, "", "#" + p.toString()); } catch (e) { /* a sandboxed host may refuse */ }
  }, 150);
}
function readHash() {
  const D = App.D; const h = location.hash.slice(1); if (!h) return; const p = new URLSearchParams(h);
  if (p.get("d") && p.get("d") !== D.man.input.sha256.slice(0, 10)) return;
  if (p.get("o")) { const [a, b] = p.get("o").split(","); const ve = Overview.viewEl(); if (a) { Overview.setOrder("taxa", a); ve.taxaSel.value = a; } if (b) { Overview.setOrder("loci", b); ve.lociSel.value = b; } }
  if (p.get("v")) { const [x0, y0, cw, ch] = p.get("v").split(",").map(Number); if ([x0, y0, cw, ch].every(Number.isFinite)) { Overview.V.x0 = x0; Overview.V.y0 = y0; Overview.V.cw = cw; Overview.V.ch = ch; Overview.V.user = true; Overview.labelSpace(); App._viewFromHash = true; } }
  if (p.get("view")) App.view = p.get("view");
  if (p.get("dv")) Detail.setView(p.get("dv"));
  if (p.get("t")) Detail.open("taxon", +p.get("t")); else if (p.get("l")) Detail.open("locus", +p.get("l")); else if (p.get("cp")) { const [m, k] = p.get("cp").split(":"); Detail.openCopies(m, k === undefined ? null : k); }
  if (p.get("c")) { const [g, li] = p.get("c").split(",").map(Number); Sequence.open(g, li); }
}
const _render = Overview.render;
Overview.render = function () { if (App._viewFromHash) { App._viewFromHash = false; const V = Overview.V; const keep = { ...V }; _render(); Object.assign(V, keep); Overview.labelSpace(); Overview.redraw(); return; } _render(); };
// the alignment panel redraws when its row is resized
new ResizeObserver(debounce(() => { const p = App.panels.sequence; if (p && p.redraw) p.redraw(); }, 120)).observe(document.body);

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
