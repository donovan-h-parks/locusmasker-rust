/* LocusMasker report viewer. One classic script: reads the bundle (a ZIP of .npy members and JSON) from the page,
   a sidecar global, a dropped file, or the website's worker, and draws everything. Nothing here computes a result:
   every number comes from the bundle the tool wrote. */
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
function ramp(cols, n) {                                   // n steps over the 3 anchor colours, as plot_panels.ramp
  if (n <= 0) return []; if (n === 1) return [cols[0]]; if (n === 2) return [cols[0], cols[2]]; if (n === 3) return cols.slice();
  return d3.range(n).map(i => { const u = 2 * i / (n - 1), j = Math.min(Math.floor(u), 1); return mixHex(cols[j], cols[j + 1], u - j); });
}
function isDark() { const t = document.documentElement.dataset.theme; return t === "dark" || (!t && matchMedia("(prefers-color-scheme: dark)").matches); }

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
  if (man.format_version > 3) throw new Error(`this report was written in bundle format ${man.format_version}; this viewer reads formats 1 to 3. Please use the viewer that came with the tool version ${man.version}.`);
  const taxa = await json("taxa.json"), loci = await json("loci.json");
  const D = { zip, man, format: man.format_version || 1, taxa, loci, n: taxa.names.length, nloc: loci.names.length, K: man.options.K_eff, M: man.options.M_eff, cols: {}, cache: { codes: new Map(), seq: new Map() } };
  const byteCols = ["n_comp", "n_fail", "n_similar", "side", "focal_sites"];
  const numCols = ["med_identity", "med_expected", "excess_m", "closest", "closest_overlap", "closest_identity", "closest_expected", "best_match", "best_match_identity", "partner", "lam", "lam_p", "row_factor"];
  for (const k of byteCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  try { D.cols.class = (await npy("cells/class.npy")).data; } catch (e) { D.cols.class = new Uint8Array(D.n * D.nloc); }   // format 1 bundles carry no class
  D.dense = man.carried.numeric_columns === "dense";
  for (const k of numCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  for (const k of ["attract", "attract_ratio", "attract_logp", "away_frac"]) { try { D.cols[k] = (await npy(`cells/${k}.npy`)).data; } catch (e) { D.cols[k] = null; } }
  D.CLASS = ["", "drift", "towards", "similar"];
  // the taxon a masked cell has moved towards: the attractor (towards), the partner (similar), else nothing
  D.movedTo = (g, li) => { const c = D.cols.class[g * D.nloc + li]; if (c === 3) { const p = D.num("partner", g, li); return p !== null && p >= 0 ? p : -1; } if (c === 2) { const a = D.num("attract", g, li); return a !== null && a >= 0 ? a : -1; } return -1; };
  D.comparators = (await npy("cells/comparators.npy")).data;
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
// 40 + s masked similar with s pairs (1..3); a cell masked on both sides draws divergent and the hover says both
const ST = { ABSENT: 0, UNTESTED: 1, GREY: 2, RED: 20, BLUE: 40 };
function deriveStates(D) {
  const { n, nloc, K, M } = D; const N = n * nloc; const st = new Uint8Array(N), sev = new Uint8Array(N);
  const fs = D.cols.focal_sites, nc = D.cols.n_comp, nf = D.cols.n_fail, ns = D.cols.n_similar, sd = D.cols.side;
  for (let i = 0; i < N; i++) {
    if (fs[i] === 0) { st[i] = ST.ABSENT; sev[i] = 0; continue; }
    if (nc[i] === 0) { st[i] = ST.UNTESTED; sev[i] = 1; continue; }
    const s = sd[i];
    if (s & 1 || s & 16) { const f = clamp(nf[i], M, K); st[i] = ST.RED + f; sev[i] = 20 + f; }
    else if (s & 2) { const k = clamp(ns[i], 1, 3); st[i] = ST.BLUE + k; sev[i] = 40 + k; }
    else { st[i] = ST.GREY; sev[i] = 2; }
  }
  D.state = st; D.sev = sev;
  const greys = ramp([cssVar("--grey-1"), cssVar("--grey-2"), cssVar("--grey-3")], M);
  const reds = ramp([cssVar("--red-1"), cssVar("--red-2"), cssVar("--red-3")], K - M + 1);
  const blues = [cssVar("--blue-1"), cssVar("--blue-2"), cssVar("--blue-3")];
  D.colorOf = code => code === ST.ABSENT ? null : code === ST.UNTESTED ? "hatch" : code < ST.RED ? greys[code - ST.GREY] : code < ST.BLUE ? reds[code - ST.RED - M] : blues[code - ST.BLUE - 1];
  D.labelOf = code => code === ST.ABSENT ? "absent" : code === ST.UNTESTED ? "present, not tested" : code < ST.RED ? "not masked" : code < ST.BLUE ? `masked, too divergent · fails ${code - ST.RED} of ${K}` : `masked, too similar · ${code - ST.BLUE}${code - ST.BLUE === 3 ? "+" : ""} pair${code - ST.BLUE > 1 ? "s" : ""}`;
  D.greys = greys; D.reds = reds; D.blues = blues;
  const surf = hexToRgb(cssVar("--surface")), hatch = hexToRgb(cssVar("--surface-3"));
  const lut = new Uint8ClampedArray(64 * 4);
  for (let c = 0; c < 64; c++) { let rgb = surf; const col = D.colorOf(c); if (col === "hatch") rgb = hatch; else if (col) rgb = hexToRgb(col); lut.set([rgb[0], rgb[1], rgb[2], 255], c * 4); }
  D.lut = lut;
  D.taxCount = { masked: new Int32Array(n), div: new Int32Array(n), sim: new Int32Array(n), tested: new Int32Array(n), present: new Int32Array(n), near: new Int32Array(n), drift: new Int32Array(n), towards: new Int32Array(n), match: new Int32Array(n).fill(-1), matchN: new Int32Array(n), matchTo: new Map() };
  D.locCount = { masked: new Int32Array(nloc), div: new Int32Array(nloc), sim: new Int32Array(nloc), tested: new Int32Array(nloc), present: new Int32Array(nloc), near: new Int32Array(nloc), drift: new Int32Array(nloc), towards: new Int32Array(nloc), match: new Int32Array(nloc).fill(-1), matchN: new Int32Array(nloc), matchTo: new Map() };
  const cls = D.cols.class;
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
    if (m) { const c = cls[i]; if (c === 1) { D.taxCount.drift[g]++; D.locCount.drift[li]++; } else if (c === 2) { D.taxCount.towards[g]++; D.locCount.towards[li]++; }
      let to = D.movedTo(g, li); if (to < 0) { const bm = D.num("best_match", g, li); to = bm !== null && bm >= 0 ? bm : -1; }      // drift cells: the best match is the least different taxon, not a relative
      if (to >= 0) { const mt = D.taxCount.matchTo, ml = D.locCount.matchTo; mt.set(g * n + to, (mt.get(g * n + to) || 0) + 1); ml.set(li * n + to, (ml.get(li * n + to) || 0) + 1); } }
  }
  for (const [cnt, N] of [[D.taxCount, n], [D.locCount, nloc]]) { for (const [key, v] of cnt.matchTo) { const u = Math.floor(key / n), to = key % n; if (v > cnt.matchN[u] || (v === cnt.matchN[u] && to < cnt.match[u])) { cnt.matchN[u] = v; cnt.match[u] = to; } } cnt.matchTo = null; }
  // the dominant class of a unit's masked cells: towards (something specific) before similar before drift
  for (const cnt of [D.taxCount, D.locCount]) { const N = cnt.masked.length; cnt.cls = new Uint8Array(N); for (let u = 0; u < N; u++) { if (!cnt.masked[u]) continue; const sim = cnt.sim[u]; cnt.cls[u] = cnt.towards[u] >= Math.max(sim, cnt.drift[u]) ? 2 : sim >= cnt.drift[u] ? 3 : 1; } }
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

// ------------------------------------------------------------------ the app
const App = { D: null, panels: {}, hover: null, sel: { rows: new Set(), cols: new Set(), cells: new Set() }, order: { taxa: null, loci: null, taxaMode: "input", lociMode: "input" }, hidePass: false, view: "summary" };

function boot() {
  const tag = document.getElementById("lm-bundle"); let b64 = (tag && tag.textContent.trim()) || (window.LM_BUNDLE_B64 || "").trim();
  if (b64 && b64 !== "@@" + "BUNDLE@@") { start(b64ToBytes(b64)).catch(showError); return; }   // the marker the tool replaces, split so this file never carries it
  if (window.LM_BUNDLE_BYTES) { start(window.LM_BUNDLE_BYTES).catch(showError); return; }
  showPicker("This page carries no data. Open a report.npz, a report.html or its report_data.js here:");
}
function showError(e) { console.error(e); $("#app").replaceChildren(el("div", { class: "error" }, "The report could not be read.\n\n" + (e && e.stack || e))); }
function showPicker(msg) {
  const inp = el("input", { type: "file", accept: ".npz,.html,.js,.zip" });
  const drop = el("div", { class: "drop" }, "drop a file here, or ", inp);
  const box = el("div", { class: "picker" }, el("h1", null, "LocusMasker report"), el("p", null, msg), drop);
  $("#app").replaceChildren(box);
  const take = async f => {
    const buf = new Uint8Array(await f.arrayBuffer());
    if (f.name.endsWith(".npz") || f.name.endsWith(".zip")) return start(buf).catch(showError);
    const text = new TextDecoder().decode(buf);
    const m = /id="lm-bundle"[^>]*>([^<]*)</.exec(text) || /LM_BUNDLE_B64\s*=\s*"([^"]+)"/.exec(text);
    if (m && m[1].trim()) return start(b64ToBytes(m[1])).catch(showError);
    showError(new Error("no bundle found in " + f.name));
  };
  inp.addEventListener("change", () => inp.files[0] && take(inp.files[0]));
  drop.addEventListener("dragover", e => { e.preventDefault(); drop.classList.add("over"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("over"));
  drop.addEventListener("drop", e => { e.preventDefault(); drop.classList.remove("over"); if (e.dataTransfer.files[0]) take(e.dataTransfer.files[0]); });
}
async function start(bytes) {
  $("#app").replaceChildren(el("div", { class: "loading" }, "Reading the report…"));
  const D = await loadBundle(bytes); App.D = D;
  App.order.taxa = d3.range(D.n); App.order.loci = d3.range(D.nloc);
  buildShell(); readHash(); renderAll();
  window.addEventListener("resize", () => { Overview.redraw(); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); });
}
window.LM = { start, App, openCell: (g, li) => Sequence.open(g, li), openUnit: (k, i) => Detail.open(k, i), setView: v => Overview.setView(v) };   // the website hands the bundle over through LM.start(bytes); the rest is for tests

// ------------------------------------------------------------------ shell: topbar, the two panel rows, hover card
function buildShell() {
  const D = App.D, m = D.man; const app = $("#app"); app.replaceChildren();
  const c = m.counts;
  const theme = el("button", { title: "light / dark", onclick: () => { const r = document.documentElement; r.dataset.theme = isDark() ? "light" : "dark"; deriveStates(D); renderAll(); } }, "◐");
  const about = el("details", { class: "meta" }, el("summary", null, `${m.tool} ${m.version} · ${m.date} · ${m.platform}`),
    el("div", null, `commit ${m.commit} · python ${m.python} · numpy ${m.numpy} · ${m.system}`),
    el("div", null, `input ${m.input.file} (${(m.input.bytes / 1e6).toFixed(1)} MB, sha256 ${m.input.sha256.slice(0, 12)}…)${m.input.partitions ? ", partitions " + m.input.partitions : ""}`),
    el("div", null, `options: datatype ${m.options.datatype}, n-comp ${m.options.n_comp} (K ${m.options.K_eff}, majority ${m.options.M_eff}), fit-overlap ${m.options.fit_overlap}, alpha ${m.options.alpha}${m.options.mask_similar_only ? ", similar side only" : ""}${m.options.mask_different_only ? ", divergent side only" : ""}`),
    el("div", null, `test: ${fmtInt(m.test.comparisons)} distinct comparisons, Bonferroni k ${m.test.bonferroni_k}, dataset dispersion ${m.test.phi_dataset}; carried: ${m.carried.sequences} sequences, codes for ${m.carried.codes} cells, best matches for ${m.carried.best_match_cells} cells`),
    el("div", null, `command: ${m.command}`));
  app.append(el("div", { class: "topbar" },
    el("h1", null, m.dataset || m.input.file), el("span", { class: "sub" }, `${fmtInt(c.masked)} of ${fmtInt(c.cells)} sequences masked (${pct(c.masked, c.cells)}) · ${fmtInt(c.taxa)} taxa × ${fmtInt(c.loci)} loci`),
    el("span", { class: "spacer" }), about, theme));
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
  let P, matrix, ctx, rowCanvas, colCanvas, strip, stripCtx, wrap, viewEl = {};
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
      el("button", { title: "export the matrix as PNG", onclick: exportPng }, "⤓ PNG"));
    viewEl.matrixCtl = [el("label", null, "taxa ", taxaSel), el("label", null, "loci ", lociSel), search];
    P.body.append(controls);
    // the summary view: stat tiles and the worst-ten lists
    viewEl.summary = el("div", { class: "summary-view", style: "display:none" });
    viewEl.tiles = el("div", { class: "tiles" }); viewEl.worst = el("div", { class: "worst" }); viewEl.summary.append(viewEl.tiles, viewEl.worst); P.body.append(viewEl.summary);
    // the matrix
    wrap = el("div", { class: "matrix-wrap" });
    rowCanvas = el("canvas", { class: "rowlabels" }); colCanvas = el("canvas", { class: "collabels" }); matrix = el("canvas", { class: "matrix", tabindex: 0 });
    wrap.append(colCanvas, rowCanvas, matrix); P.body.append(wrap); ctx = matrix.getContext("2d");
    const legend = el("div", { class: "legend" });
    legend.append(el("span", { class: "fam", title: "tested and not masked" }, sw(D.greys[0]), "not masked"),
      el("span", { class: "fam" }, el("span", { class: "lab" }, "masked, too divergent, fails"), ...d3.range(D.M, D.K + 1).map(f => el("span", { class: "fam", title: D.labelOf(ST.RED + f) }, sw(D.reds[f - D.M]), String(f)))),
      el("span", { class: "fam" }, el("span", { class: "lab" }, "masked, too similar, pairs"), ...[1, 2, 3].map(k => el("span", { class: "fam", title: D.labelOf(ST.BLUE + k) }, sw(D.blues[k - 1]), k === 3 ? "3+" : String(k)))),
      el("span", { class: "fam", title: "present, not tested" }, sw("hatch"), "not tested"), el("span", { class: "fam", title: "no residue" }, sw(null), "absent"), el("span", { class: "lab" }, `of ${D.K} comparators`));
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
    const D = App.D, c = D.man.counts; const t = viewEl.tiles; t.replaceChildren();
    const tile = (v, k, title) => el("div", { class: "tile", title: title || "" }, el("div", { class: "v" }, v), el("div", { class: "k" }, k));
    const taxaWith = d3.sum(D.taxCount.masked, v => v > 0), lociWith = d3.sum(D.locCount.masked, v => v > 0);
    t.append(tile(fmtInt(c.masked), `sequences masked · ${pct(c.masked, c.cells)} of ${fmtInt(c.cells)} tested`),
      tile(fmtInt(c.masked_divergent), "too divergent" + (() => { const dr = d3.sum(D.taxCount.drift), tw = d3.sum(D.taxCount.towards); return dr + tw ? ` · ${fmtInt(dr)} drift, ${fmtInt(tw)} towards a taxon` : ""; })() + (c.ambiguous ? ` · ${c.ambiguous} with fewer comparators` : "")),
      tile(fmtInt(c.masked_similar), "too similar" + (c.masked_both ? ` · ${c.masked_both} both sides` : "")),
      tile(`${fmtInt(taxaWith)} / ${fmtInt(c.taxa)}`, "taxa with a masked sequence" + (c.excluded_taxa ? ` · ${c.excluded_taxa} excluded` : "")),
      tile(`${fmtInt(lociWith)} / ${fmtInt(c.loci)}`, "loci with a masked sequence" + (c.skipped_loci ? ` · ${c.skipped_loci} skipped` : "")));
  }
  function renderWorst() {
    const D = App.D; const w = viewEl.worst; w.replaceChildren();
    const col = (title, items, open) => {
      const c = el("div", { class: "col" }, el("h3", null, title));
      const max = d3.max(items, d => d.v) || 1;
      for (const it of items) {
        const bar = el("div", { class: "bar" });
        const div = it.div / (it.div + it.sim || 1);
        bar.append(el("i", { style: `width:${100 * it.v / max}%; background: linear-gradient(90deg, ${D.reds[D.reds.length - 2] || D.reds[0]} ${100 * div}%, ${D.blues[1]} ${100 * div}%)` }));
        c.append(el("div", { class: "row", title: it.title, onclick: () => open(it.i), onmouseenter: () => highlight(it.kind, it.i), onmouseleave: () => highlight(null) },
          el("span", { class: "nm" }, it.name), bar, el("span", { class: "num" }, it.label)));
      }
      return c;
    };
    const taxa = d3.range(D.n).filter(g => D.taxCount.masked[g] > 0).map(g => ({ i: g, kind: "row", name: D.taxa.names[g], v: D.taxCount.masked[g] / Math.max(D.taxCount.tested[g], 1), div: D.taxCount.div[g], sim: D.taxCount.sim[g], label: pct(D.taxCount.masked[g], D.taxCount.tested[g]), title: `${D.taxCount.masked[g]} of ${D.taxCount.tested[g]} tested loci masked` }))
      .sort((a, b) => b.v - a.v || b.div + b.sim - a.div - a.sim).slice(0, 10);
    const loci = d3.range(D.nloc).filter(l => D.locCount.masked[l] > 0).map(l => ({ i: l, kind: "col", name: D.loci.names[l], v: D.locCount.masked[l], div: D.locCount.div[l], sim: D.locCount.sim[l], label: fmtInt(D.locCount.masked[l]), title: `${D.locCount.masked[l]} of ${D.locCount.tested[l]} tested taxa masked` }))
      .sort((a, b) => b.v - a.v).slice(0, 10);
    const ct = col("most masked taxa (share of tested loci)", taxa, g => Detail.open("taxon", g)), cl = col("most masked loci (sequences)", loci, l => Detail.open("locus", l));
    w.append(ct, cl);
    // the histograms under the lists, sized to the room left in the panel
    const room = P.body.getBoundingClientRect().height - viewEl.tiles.getBoundingClientRect().height - w.getBoundingClientRect().height - 60;
    const hh = clamp(room, 190, 420), hw = (w.getBoundingClientRect().width || 800) / 2 - 20;
    const ht = el("div", { class: "hist-wrap" }), hl = el("div", { class: "hist-wrap" }); ct.append(ht); cl.append(hl);
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
    const img = ctx.createImageData(matrix.width, matrix.height); const px = img.data; const lut = D.lut; const sev = D.sev, st = D.state;
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
      }
    }
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
    const zoomAt = (x, y, fx, fy) => { if (fx !== 1) { const ncw = clamp(V.cw * fx, 0.05, 40); V.x0 = V.x0 + x / V.cw - x / ncw; V.cw = ncw; } if (fy !== 1) { const nch = clamp(V.ch * fy, 0.05, 40); V.y0 = V.y0 + y / V.ch - y / nch; V.ch = nch; } V.user = true; labelSpace(); redraw(); writeHash(); };
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
    const D = App.D; App.sel.rows.clear(); App.sel.cols.clear(); q = q.trim().toLowerCase(); if (!q) { redraw(); return; }
    let first = null;
    D.taxa.names.forEach((nm, g) => { if (nm.toLowerCase().includes(q)) { App.sel.rows.add(g); if (!first) first = ["row", g]; } });
    D.loci.names.forEach((nm, l) => { if (nm.toLowerCase().includes(q)) { App.sel.cols.add(l); if (!first) first = ["col", l]; } });
    if (first) { const { w, h } = canvasArea(); if (first[0] === "row") { const r = App.order.taxa.indexOf(first[1]); V.y0 = r - h / V.ch / 2; } else { const c = App.order.loci.indexOf(first[1]); V.x0 = c - w / V.cw / 2; } }
    redraw();
  }
  return { create, render, redraw, setOrder, setView, refresh, V, highlight, labelSpace, fitAll, viewEl: () => viewEl };
})();
function debounce(fn, ms) { let t = null; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; }

// hover content for a cell, shared by the matrix, the ribbons and the table
function cellHover(g, li, compact) {
  const D = App.D; const i = g * D.nloc + li; const code = D.state[i]; const col = D.colorOf(code);
  const t = el("table"); const row = (k, v) => { if (v === null || v === undefined || v === "") return; t.append(el("tr", null, el("td", { class: "k" }, k), el("td", null, v))); };
  const head = el("div", { class: "state" }, sw(col), `${D.taxa.names[g]} · ${D.loci.names[li]}`);
  const box = el("div", null, head);
  const reason = stateReason(D, g, li);
  const sd = D.cols.side[i], nf = D.cols.n_fail[i], ns = D.cols.n_similar[i], nc = D.cols.n_comp[i];
  let label = D.labelOf(code); if ((sd & 3) === 3) label += " · also masked as too similar"; if (sd & 4) label += " · fewer comparators than requested"; if (sd & 16) label += " (flagged; divergent masking was off)"; if (sd & 32) label += " · a similar-side fail (similar masking was off)";
  row("state", reason ? `${D.labelOf(code)}: ${reason}` : label);
  if (nc > 0) {
    row("residues", fmtInt(D.cols.focal_sites[i]));
    row("comparators", `${nc} tested, ${nf} failed as too divergent${ns ? `, ${ns} as too similar` : ""}`);
    const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li), xm = D.num("excess_m", g, li);
    if (mi !== null) row("divergence, median over comparators", `${fmt2(100 - mi)}% observed, ${fmt2(100 - me)}% expected`);
    const cl = D.num("closest", g, li);
    if (cl !== null && cl >= 0 && !compact) row("closest comparator", `${D.taxa.names[cl]}: ${fmt2(100 - D.num("closest_identity", g, li))}% observed vs ${fmt2(100 - D.num("closest_expected", g, li))}% expected over ${fmtInt(D.num("closest_overlap", g, li))} shared sites`);
    if (xm !== null) row("excess beyond expectation", `${fmt2(xm)} m`);
    const bm = D.num("best_match", g, li);
    if (bm !== null && bm >= 0 && !compact) { const close = Array.from(D.closest.subarray(g * 5, g * 5 + 5)).includes(bm); row("best match", `${D.taxa.names[bm]} (${fmt1(D.num("best_match_identity", g, li))}% identity over shared sites${close ? ", one of the 5 closest relatives" : ""})`); }
    const cc = D.cols.class[i]; const aw = D.num("away_frac", g, li);
    if (cc === 1) row("class", `drift: further than expected from ${aw === null ? "its" : pct(aw, 1)} of the taxa it shares sites with, and not significantly closer to any${(() => { const a = D.num("attract", g, li); return a !== null && a >= 0 ? ` (closest to a deficit: ${D.taxa.names[a]}, ${fmt2(D.num("attract_ratio", g, li))} × expected, log p ${fmt1(D.num("attract_logp", g, li))})` : ""; })()}`);
    else if (cc === 2) { const a = D.num("attract", g, li); row("class", `towards ${D.taxa.names[a]}: ${fmt2(D.num("attract_ratio", g, li))} × the divergence the genome predicts for the pair (log p ${fmt1(D.num("attract_logp", g, li))}, past the bar); further than expected from ${aw === null ? "most" : pct(aw, 1)} of the other taxa`); }
    const pt = D.num("partner", g, li);
    if (pt !== null && pt >= 0) { const lam = D.num("lam", g, li), lp = D.num("lam_p", g, li); row("too similar to", `${D.taxa.names[pt]}${sd & 8 ? " (identified as the recipient)" : " (direction not identified)"}`); if (lam !== null) row("direction LRT", `|2 ln LR| = ${fmt1(lam)}, p = ${lp < 1e-4 ? lp.toExponential(1) : fmt3(lp)} (χ², 1 df)`); const rf = D.num("row_factor", g, li); if (rf !== null && rf < 1) row("row factor", `${fmt3(rf)} (the pair was tested against ${fmt3(rf)} × its expectation: the sequence sits below expectation against its other comparators too)`); }
    if (!compact) row("", el("span", { class: "more" }, "click: open the alignment"));
  }
  box.append(t); return box;
}
function unitHover(kind, i) {
  const D = App.D; const c = kind === "taxon" ? D.taxCount : D.locCount; const nm = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
  const t = el("table"); const row = (k, v) => t.append(el("tr", null, el("td", { class: "k" }, k), el("td", null, v)));
  row(kind, nm); row("status", kind === "taxon" ? D.taxa.status[i] : D.loci.status[i]);
  row("sequences", `${fmtInt(c.present[i])} present, ${fmtInt(c.tested[i])} tested`);
  row("masked", `${fmtInt(c.masked[i])} (${pct(c.masked[i], c.tested[i])}): ${c.div[i]} too divergent, ${c.sim[i]} too similar; ${c.near[i]} near misses`);
  if (c.masked[i]) { row("class", `${c.drift[i]} drift, ${c.towards[i]} towards, ${c.sim[i]} similar`); if (c.match[i] >= 0) row("closest match", `${D.taxa.names[c.match[i]]} for ${c.matchN[i]} of the masked sequences`); }
  if (kind === "locus") { row("length", `${fmtInt(D.loci.length[i])} sites`); if (D.loci.rate[i] !== null) row("rate", `${fmt3(D.loci.rate[i])} × the typical locus`); row("plateau", `${fmt1(D.loci.plateau[i])}%`); if (D.loci.phi[i] !== null) row("dispersion φ", fmt3(D.loci.phi[i])); }
  else { const cl = Array.from(D.closest.subarray(i * 5, i * 5 + 5)).filter(x => x >= 0).map(x => D.taxa.names[x]); if (cl.length) row("closest by m", cl.join(", ")); }
  return el("div", null, t, el("div", { class: "more" }, "click: open the panel"));
}

// ------------------------------------------------------------------ residual ribbons: the cells of one unit on an asinh axis
// draws the residuals of the cells of a locus (by="locus") or a taxon into a canvas row; returns the drawn points for hover
function drawResidualRow(cv, kind, u, W, H, opts) {
  const D = App.D; const dpr = devicePixelRatio || 1; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + "px"; cv.style.height = H + "px";
  const c = cv.getContext("2d"); c.scale(dpr, dpr); c.clearRect(0, 0, W, H);
  const xs = d3.scaleSymlog().constant(1).domain(D.exDomain).range([opts.padL || 6, W - (opts.padR || 6)]);
  const pts = []; let identical = 0; const N = kind === "locus" ? D.n : D.nloc;
  for (let k = 0; k < N; k++) {
    const g = kind === "locus" ? k : u, li = kind === "locus" ? u : k; const i = g * D.nloc + li; if (D.cols.n_comp[i] === 0) continue;
    const x = D.num("excess_m", g, li); if (x === null || !Number.isFinite(x)) continue;
    const mi = D.num("med_identity", g, li); if (mi !== null && 100 - mi <= 0) { identical++; continue; }
    pts.push({ g, li, x, code: D.state[i], px: xs(x), py: 0 });
  }
  // the saturation band for a locus: where a saturated sequence would land, from the plateau, the rate and the expectation
  if (kind === "locus" && D.loci.rate[u]) {
    const P = D.loci.plateau[u], rate = D.loci.rate[u]; const corr = p => -P * Math.log(1 - Math.min(Math.max(p, 0), 0.98 * P) / P);
    const ce = pts.map(p => { const me = D.num("med_expected", p.g, p.li); return me === null ? NaN : (corr(0.98 * P) - corr(100 - me)) / rate; }).filter(Number.isFinite).sort(d3.ascending);
    if (ce.length) { const c10 = d3.quantileSorted(ce, 0.1); c.fillStyle = cssVar("--surface-3"); c.fillRect(xs(c10), 1, Math.max(0, W - xs(c10)), H - 2); }
  }
  c.strokeStyle = cssVar("--model"); c.lineWidth = 1; c.beginPath(); c.moveTo(xs(0) + 0.5, 1); c.lineTo(xs(0) + 0.5, H - 1); c.stroke();
  const jit = d3.randomUniform.source(d3.randomLcg(0.31 + (u % 97) / 100))(-(H / 2 - 4), H / 2 - 4); const cy = H / 2;
  pts.sort((a, b) => D.sev[a.g * D.nloc + a.li] - D.sev[b.g * D.nloc + b.li]);
  const many = pts.length > 400;
  for (const p of pts) {
    p.py = cy + jit(); const col = D.colorOf(p.code) || cssVar("--grey-1");
    if (p.code >= ST.RED) { const r = 4.2; c.fillStyle = col; c.strokeStyle = cssVar("--surface"); c.lineWidth = 0.8; c.beginPath();
      if (p.code >= ST.BLUE) { c.moveTo(p.px, p.py + r); c.lineTo(p.px - r, p.py - r * 0.7); c.lineTo(p.px + r, p.py - r * 0.7); } else { c.moveTo(p.px, p.py - r); c.lineTo(p.px - r, p.py + r * 0.7); c.lineTo(p.px + r, p.py + r * 0.7); }
      c.closePath(); c.fill(); c.stroke(); }
    else { c.globalAlpha = many ? 0.35 : 0.7; c.fillStyle = col; c.beginPath(); c.arc(p.px, p.py, many ? 1.4 : 2, 0, 2 * Math.PI); c.fill(); c.globalAlpha = 1; }
  }
  return { pts, identical, xs };
}
function residualAxis(W, padL, padR) {
  const D = App.D; const xs = d3.scaleSymlog().constant(1).domain(D.exDomain).range([padL, W - padR]);
  const svg = d3.create("svg").attr("width", W).attr("height", 30).attr("class", "resid-axis");
  const tv = [-1000, -100, -10, -1, 0, 1, 10, 100, 1000].filter(v => v >= D.exDomain[0] && v <= D.exDomain[1]);
  svg.append("g").attr("class", "axis").attr("transform", "translate(0,16)").call(d3.axisTop(xs).tickValues(tv).tickFormat(d3.format(",")));
  svg.append("text").attr("x", xs(0)).attr("y", 29).attr("text-anchor", "middle").attr("fill", cssVar("--muted")).attr("font-size", 10.5).text("◀ too similar · excess divergence in units of m · too divergent ▶");
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
  let root, kind, sortMode = "masked", units = [], list, inner, W = 600, axisBox;
  const RH = 22, NAME_W = 190 + 52 + 150 + 170 + 16;                            // name, count, class, closest match, the grid gaps
  const CLS_ORDER = { 2: 0, 3: 1, 1: 2, 0: 3 };                                  // towards, similar, drift, none
  let ctlHost = null;
  function render(container, k, host) { root = container; kind = k; ctlHost = host; root.replaceChildren(); build(); }
  function build() {
    const D = App.D; const cnt = kind === "locus" ? D.locCount : D.taxCount; const names = kind === "locus" ? D.loci.names : D.taxa.names; const N = kind === "locus" ? D.nloc : D.n;
    units = d3.range(N).filter(u => cnt.tested[u] > 0);
    const T = D.taxa.names; const matchName = u => cnt.match[u] >= 0 ? T[cnt.match[u]] : "\uffff";
    const byClass = (a, b) => CLS_ORDER[cnt.cls[a]] - CLS_ORDER[cnt.cls[b]] || cnt.towards[b] - cnt.towards[a] || cnt.masked[b] - cnt.masked[a];
    const byMatch = (a, b) => matchName(a).localeCompare(matchName(b)) || cnt.matchN[b] - cnt.matchN[a] || cnt.masked[b] - cnt.masked[a];
    if (sortMode === "masked") units.sort((a, b) => cnt.masked[b] - cnt.masked[a] || cnt.near[b] - cnt.near[a] || a - b); else if (sortMode === "name") units.sort((a, b) => names[a].localeCompare(names[b]));
    else if (sortMode === "class") units.sort((a, b) => byClass(a, b) || byMatch(a, b) || a - b); else if (sortMode === "match") units.sort((a, b) => byMatch(a, b) || byClass(a, b) || a - b);
    const sortSel = el("select", { onchange: e => { sortMode = e.target.value; build(); } }, ...[["masked", "most masked first"], ["class", "class, then closest match"], ["match", "closest match, then class"], ["name", "name"], ["input", "input order"]].map(([v, t]) => el("option", { value: v, selected: v === sortMode ? "" : null }, t)));
    if (ctlHost) { ctlHost.replaceChildren(el("label", null, "sort ", sortSel), el("span", { class: "k", title: `each row is one ${kind}, each mark one of its sequences: ▲ masked too divergent, ▼ too similar, dots pass${kind === "locus" ? "; the shaded band is where a saturated sequence lands" : ""}. Click a name for its panel, a mark for the alignment.` }, `${fmtInt(units.length)} ${kind === "locus" ? "loci" : "taxa"}`)); }
    axisBox = el("div", { class: "rib-axis" }); list = el("div", { class: "rib-list" }); inner = el("div", { style: "position:relative" }); list.append(inner);
    root.replaceChildren(axisBox, list);
    W = Math.max(300, root.getBoundingClientRect().width - NAME_W - 24);
    axisBox.replaceChildren(el("div", { class: "rib-head", style: `width:${NAME_W}px; flex: 0 0 auto` }, el("span", null, kind === "locus" ? "locus" : "taxon"), el("span", null, "masked"), el("span", { title: "the classes of the unit's masked sequences: drift (further from everything), towards (closer than the genome predicts to one taxon), similar (too similar to a partner)" }, "class"), el("span", { title: "the taxon most of the unit's masked sequences moved towards (or, for drift, their least different taxon)" }, "closest match")), residualAxis(W, 8, 8));
    inner.style.height = units.length * RH + "px";
    list.addEventListener("scroll", paint); paint();
  }
  function paint() {
    const D = App.D; const cnt = kind === "locus" ? D.locCount : D.taxCount; const names = kind === "locus" ? D.loci.names : D.taxa.names;
    inner.replaceChildren(); const top = list.scrollTop; const i0 = Math.max(0, Math.floor(top / RH) - 3), i1 = Math.min(units.length, i0 + Math.ceil(list.clientHeight / RH) + 6);
    for (let i = i0; i < i1; i++) {
      const u = units[i]; const cv = el("canvas"); let drawn = null;
      const nm = el("span", { class: "nm", title: names[u], onclick: () => Detail.open(kind, u) }, sw(cnt.masked[u] ? (cnt.div[u] >= cnt.sim[u] ? D.reds[1] || D.reds[0] : D.blues[1]) : null), names[u]);
      const clsCell = el("span", { class: "cls" }); if (cnt.masked[u]) { const parts = []; if (cnt.towards[u]) parts.push(el("span", { style: `color:${D.reds[D.reds.length - 1]}` }, `towards ${cnt.towards[u]}`)); if (cnt.drift[u]) parts.push(el("span", null, `drift ${cnt.drift[u]}`)); if (cnt.sim[u]) parts.push(el("span", { style: `color:${D.blues[1]}` }, `similar ${cnt.sim[u]}`)); parts.forEach((p, j) => { if (j) clsCell.append(" · "); clsCell.append(p); }); }
      const mt = el("span", { class: "nm match", title: cnt.match[u] >= 0 ? `${D.taxa.names[cnt.match[u]]}: ${cnt.matchN[u]} of ${cnt.masked[u]} masked sequences` : "" }, cnt.match[u] >= 0 ? `${D.taxa.names[cnt.match[u]]} (${cnt.matchN[u]})` : "");
      if (cnt.match[u] >= 0) mt.onclick = () => Detail.open("taxon", cnt.match[u]);
      const row = el("div", { class: "rib-row", style: `position:absolute; top:${i * RH}px; left:0; right:0; height:${RH}px` }, nm, el("span", { class: "cnt" }, cnt.masked[u] ? fmtInt(cnt.masked[u]) : ""), clsCell, mt, cv);
      inner.append(row);
      drawn = drawResidualRow(cv, kind, u, W, RH - 2, { padL: 8, padR: 8 });
      if (drawn.identical) row.title = `${drawn.identical} sequences identical to their comparators are not drawn`;
      bindResidualHover(cv, () => drawn.pts);
    }
  }
  return { render };
})();

// ------------------------------------------------------------------ the taxon / locus panel: one plot at a time, sized to the panel
const Detail = (() => {
  let cur = null, view = "scatter", viewSel = null, body = null;
  function open(kind, i) {
    cur = { kind, i };
    if (!App.panels.detail) {
      const P = makePanel("detail", [], () => closePanel("detail")); App.panels.detail = { el: P.el, P, render, redraw: render };
      App.topEl.append(App.vsplitEl, P.el); App.vsplitEl.style.display = "";
      new ResizeObserver(debounce(() => { if (App.panels.detail && cur) draw(); }, 120)).observe(P.body);
    }
    if (kind === "locus" && view === "matches") view = "scatter";
    render(); writeHash(); Overview.highlight(kind === "taxon" ? "row" : "col", i);
  }
  function cellsOf() {
    const D = App.D; const out = [];
    if (cur.kind === "taxon") { const g = cur.i; for (let li = 0; li < D.nloc; li++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    else { const li = cur.i; for (let g = 0; g < D.n; g++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    return out;
  }
  function area() { const r = body.getBoundingClientRect(); return { W: Math.max(300, r.width - 20), H: Math.max(220, r.height - 20) }; }
  function render() {
    if (!cur || !App.panels.detail) return;
    const D = App.D; const { kind, i } = cur; const P = App.panels.detail.P; const name = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
    const cnt = kind === "taxon" ? D.taxCount : D.locCount;
    const views = [["scatter", "observed vs expected"], ...(kind === "taxon" ? [["matches", "best matches"]] : []), ["ribbons", "site ribbons"]];
    viewSel = el("select", { onchange: e => { view = e.target.value; draw(); } }, ...views.map(([v, t]) => el("option", { value: v, selected: v === view ? "" : null }, t)));
    P.setCrumb([el("span", null, kind, " ", el("b", null, name)), el("span", { class: "sub" }, `${fmtInt(cnt.masked[i])} of ${fmtInt(cnt.tested[i])} sequences masked (${pct(cnt.masked[i], cnt.tested[i])})`), el("label", { class: "crumb-ctl" }, "view ", viewSel)]);
    P.body.replaceChildren(); body = P.body; body.classList.add("detail-body");
    draw();
  }
  function draw() {
    if (!cur || !App.panels.detail) return; body.replaceChildren();
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
    const line = (f, dash, col) => { const [x0, x1] = xs.domain(); const p0 = [Math.max(x0, ys.domain()[0] / f), Math.min(x1, ys.domain()[1] / f)]; if (p0[0] < p0[1]) svg.append("line").attr("x1", xs(p0[0])).attr("y1", ys(p0[0] * f)).attr("x2", xs(p0[1])).attr("y2", ys(p0[1] * f)).attr("stroke", col).attr("stroke-width", f === 1 ? 1.2 : 0.9).attr("stroke-dasharray", dash); };
    line(1, null, cssVar("--model"));
    wrap.append(svg.node());
    const cv = el("canvas"); const dpr = devicePixelRatio || 1; cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px"; wrap.append(cv);
    const c = cv.getContext("2d"); c.scale(dpr, dpr);
    const drawPt = (p, r) => { const x = xs(X(p)), y = ys(Y(p)); const col = D.colorOf(p.code) || cssVar("--grey-1"); c.fillStyle = col; c.strokeStyle = cssVar("--surface"); c.lineWidth = 0.8; c.beginPath();
      if (p.code >= ST.BLUE) { c.moveTo(x, y + r); c.lineTo(x - r, y - r * 0.7); c.lineTo(x + r, y - r * 0.7); c.closePath(); } else if (p.code >= ST.RED) { c.moveTo(x, y - r); c.lineTo(x - r, y + r * 0.7); c.lineTo(x + r, y + r * 0.7); c.closePath(); } else c.arc(x, y, r * 0.75, 0, 2 * Math.PI);
      c.fill(); c.stroke(); };
    const order = pts.slice().sort((a, b) => D.sev[a.g * D.nloc + a.li] - D.sev[b.g * D.nloc + b.li]);
    const r = pts.length > 2000 ? 2.4 : pts.length > 500 ? 3.2 : 4.2;
    for (const p of order) drawPt(p, p.code >= ST.RED ? r + 1.5 : r);
    svg.style("pointer-events", "all");
    const nearest = ev => { const [mx, my] = d3.pointer(ev, svg.node()); let best = null, bd = 81; for (const p of pts) { const dx = xs(X(p)) - mx, dy = ys(Y(p)) - my; const d = dx * dx + dy * dy; if (d < bd) { bd = d; best = p; } } return best; };
    svg.on("mousemove", ev => { const p = nearest(ev); if (p) { showHover(cellHover(p.g, p.li), ev.clientX, ev.clientY); App.sel.cells.clear(); App.sel.cells.add(p.g * D.nloc + p.li); Overview.redraw(); svg.style("cursor", "pointer"); } else { hideHover(); svg.style("cursor", null); } });
    svg.on("mouseleave", () => { hideHover(); App.sel.cells.clear(); Overview.redraw(); });
    svg.on("click", ev => { const p = nearest(ev); if (p) Sequence.open(p.g, p.li); });
    box.append(el("div", { class: "note" }, `each point is one sequence. Diagonal: observed = expected. ● not masked, ▲ masked too divergent, ▼ masked too similar. Hover a point for its numbers; click it for the alignment.`));
    return box;
  }
  // best matches of the masked sequences, with the count expected from the unmasked ones
  function matchesPlot(cells) {
    const D = App.D; const g0 = cur.i; const box = el("div", { class: "plot" });
    const close = new Set(Array.from(D.closest.subarray(g0 * 5, g0 * 5 + 5)).filter(x => x >= 0));
    const masked = [], unmasked = [];
    for (const [g, li] of cells) { const bm = D.num("best_match", g, li); if (bm === null || bm < 0) continue; ((D.cols.side[g * D.nloc + li] & 3) ? masked : unmasked).push({ g, li, bm, id: D.num("best_match_identity", g, li), code: D.state[g * D.nloc + li] }); }
    if (!masked.length) { box.append(el("div", { class: "note" }, unmasked.length ? "no masked sequence for this taxon" : "best matches were not carried for this taxon (the dataset is large: they are carried for cells with a fail only)")); return box; }
    const cm = d3.rollup(masked, v => v.length, d => d.bm), cu = d3.rollup(unmasked, v => v.length, d => d.bm); const nm = masked.length, nu = unmasked.length;
    let taxaSorted = Array.from(new Set([...cm.keys(), ...Array.from(cu.keys()).filter(u => (cu.get(u) / Math.max(nu, 1)) * nm >= 1)])).sort((a, b) => (cm.get(b) || 0) - (cm.get(a) || 0) || (cu.get(b) || 0) - (cu.get(a) || 0));
    const { W, H: H0 } = area(); const H = H0 - 60; const MAXB = Math.max(6, Math.floor((W - 60) / 34));
    let other = null; if (taxaSorted.length > MAXB) { const rest = taxaSorted.slice(MAXB - 1); taxaSorted = taxaSorted.slice(0, MAXB - 1); other = { n: rest.length, cm: d3.sum(rest, u => cm.get(u) || 0), cu: d3.sum(rest, u => cu.get(u) || 0), set: new Set(rest) }; }
    const m = { l: 44, r: 10, t: 14, b: Math.min(130, H * 0.42) };
    const cats = taxaSorted.map(u => ({ key: u, name: D.taxa.names[u], cm: cm.get(u) || 0, cu: cu.get(u) || 0, rows: masked.filter(d => d.bm === u) })); if (other) cats.push({ key: -1, name: `other (${other.n} taxa)`, cm: other.cm, cu: other.cu, rows: masked.filter(d => other.set.has(d.bm)) });
    const xs = d3.scaleBand().domain(cats.map(c => c.name)).range([m.l, W - m.r]).padding(0.18); const top = Math.max(d3.max(cats, c => Math.max(c.cm, c.cu / Math.max(nu, 1) * nm)), 1);
    const ys = d3.scaleLinear().domain([0, top * 1.15]).range([H - m.b, m.t]);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).tickValues(d3.range(0, top + 1, Math.max(1, Math.ceil(top / 4)))).tickFormat(d3.format("d")));
    const ax = svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs));
    ax.selectAll("text").attr("transform", "rotate(-60)").attr("text-anchor", "end").attr("dx", "-0.6em").attr("dy", "0.2em").style("font-weight", d => cats.find(c => c.name === d) && close.has(cats.find(c => c.name === d).key) ? "700" : "400").style("fill", cssVar("--ink"));
    for (const c of cats) {
      let y0 = H - m.b; const byCode = d3.rollup(c.rows, v => v.length, d => d.code); const codes = Array.from(byCode.keys()).sort(d3.ascending);
      for (const code of codes) { const h = ys(0) - ys(byCode.get(code)); svg.append("rect").attr("x", xs(c.name)).attr("y", y0 - h).attr("width", xs.bandwidth()).attr("height", Math.max(h, 0)).attr("fill", D.colorOf(code)).attr("stroke", cssVar("--surface")).attr("stroke-width", 0.8).style("cursor", "pointer")
        .on("mousemove", ev => showHover(el("div", null, el("b", null, c.name), el("div", null, `${byCode.get(code)} masked sequence${byCode.get(code) > 1 ? "s" : ""} match it best (${D.labelOf(code)})`), el("div", { class: "more" }, `${c.cm} masked in all; ${fmt1(c.cu / Math.max(nu, 1) * nm)} expected from the unmasked sequences`)), ev.clientX, ev.clientY)).on("mouseleave", hideHover)
        .on("click", () => { App.sel.cells.clear(); for (const d of c.rows) App.sel.cells.add(d.g * D.nloc + d.li); Overview.redraw(); }); y0 -= h; }
      const e = c.cu / Math.max(nu, 1) * nm; svg.append("line").attr("x1", xs(c.name) - 1).attr("x2", xs(c.name) + xs.bandwidth() + 1).attr("y1", ys(e)).attr("y2", ys(e)).attr("stroke", cssVar("--ink")).attr("stroke-width", 1.4);
      if (c.cm) svg.append("text").attr("x", xs(c.name) + xs.bandwidth() / 2).attr("y", ys(c.cm) - 3).attr("text-anchor", "middle").style("fill", cssVar("--ink")).text(c.cm);
    }
    box.append(svg.node());
    const pm = 100 * d3.mean(masked, d => close.has(d.bm) ? 1 : 0), pu = nu ? 100 * d3.mean(unmasked, d => close.has(d.bm) ? 1 : 0) : NaN;
    box.append(el("div", { class: "note" }, `masked sequences by the taxon each matches best at its locus, stacked by mask kind; the tick is the count expected if they matched like the taxon's ${fmtInt(nu)} unmasked sequences; bold names are the 5 closest relatives genome-wide. Best match is one of the 5 closest for ${fmt1(pm)}% of masked and ${Number.isFinite(pu) ? fmt1(pu) : "–"}% of unmasked sequences; median identity to it ${fmt1(d3.median(masked, d => d.id))}% masked vs ${nu ? fmt1(d3.median(unmasked, d => d.id)) : "–"}% unmasked. Bars on their ticks are damaged copies of the right sequence (a wrong exon); bars shifted onto a distant taxon point to contamination or a shared paralogue.`));
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
    for (const [g, li] of masked) {
      const cv = el("canvas"); const row = el("div", { class: "rib-row2", onclick: () => Sequence.open(g, li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(g * D.nloc + li); Overview.redraw(); }, onmouseleave: () => { App.sel.cells.clear(); Overview.redraw(); } },
        el("span", { class: "nm", title: cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g] }, sw(D.colorOf(D.state[g * D.nloc + li])), cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g]), cv);
      wrap.append(row);
      D.codes(li).then(C => drawRibbon(cv, C, g, li, maxLen));
    }
    box.append(el("div", { class: "note" }, `${fmtInt(masked.length)} masked sequences, sites along x: red where a majority of the comparators differ, blue where a pair shares a state against the others. ` + (cur.kind === "locus" ? "Rows share the alignment's columns: a vertical band is a column problem (misalignment), a horizontal block in one row is a taxon problem (a wrong exon, a frameshift)." : "Each row is a different locus drawn to its own length on a shared axis.") + " Hover a row for the site; click it for the alignment."));
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
  return { open, render, cur: () => cur, setView: v => { view = v; } };
})();

// site-code helpers shared by the ribbons and the alignment panel
function codeSeverity(D, c) { if (c === 255) return 0; if (c === 254) return 5; if (c >= 128) return 30 + (c - 128); if (c >= D.M) return 20 + c; return 2 + c; }
function codeHex(D, c) {
  if (c === 255) return null; if (c === 254) return cssVar("--muted");
  if (c >= 128) { const j = c - 128; return D.blues[clamp(j - 2, 0, 2)]; }
  if (c >= D.M) return D.reds[Math.min(c - D.M, D.reds.length - 1)];
  return D.greys[Math.min(c, D.greys.length - 1)];
}
function codeRgb(D, c) { const h = codeHex(D, c); return hexToRgb(h || cssVar("--surface")); }
function codeLabel(D, c) { if (c === 255) return "no residue in the focal sequence"; if (c === 254) return "the pair differs here"; if (c >= 128) return `the pair shares this state while ${c - 128} other comparators differ`; return `${c} of the comparators present differ${c >= D.M ? " (a majority)" : ""}`; }
function alignmentPos(D, li, s) { let left = s; for (const [a, b, step] of D.loci.ranges[li]) { const n = Math.floor((b - a) / step) + 1; if (left < n) return a + left * step; left -= n; } return s; }

// ------------------------------------------------------------------ residue colour schemes (Jalview's standard tables: www.jalview.org/help/html/colourSchemes)
const RESIDUE_SCHEMES = {
  zappo: { I: "#ffafaf", L: "#ffafaf", V: "#ffafaf", A: "#ffafaf", M: "#ffafaf", F: "#ffc800", W: "#ffc800", Y: "#ffc800", K: "#6464ff", R: "#6464ff", H: "#6464ff", D: "#ff0000", E: "#ff0000", S: "#00ff00", T: "#00ff00", N: "#00ff00", Q: "#00ff00", P: "#ff00ff", G: "#ff00ff", C: "#ffff00" },
  taylor: { A: "#ccff00", V: "#99ff00", I: "#66ff00", L: "#33ff00", M: "#00ff00", F: "#00ff66", Y: "#00ffcc", W: "#00ccff", H: "#0066ff", R: "#0000ff", K: "#6600ff", N: "#cc00ff", Q: "#ff00cc", E: "#ff0066", D: "#ff0000", S: "#ff3300", T: "#ff6600", G: "#ff9900", P: "#ffcc00", C: "#ffff00" },
  hydrophobic: { I: "#ff0000", V: "#f60009", L: "#ea0015", F: "#cb0034", C: "#c2003d", M: "#b0004f", A: "#ad0052", G: "#6a0095", X: "#680097", T: "#61009e", S: "#5e00a1", W: "#5b00a4", Y: "#4f00b0", P: "#4600b9", H: "#1500ea", E: "#0c00f3", Z: "#0c00f3", Q: "#0c00f3", D: "#0c00f3", B: "#0c00f3", N: "#0c00f3", K: "#0000ff", R: "#0000ff" },
  helix: { E: "#ff00ff", M: "#ef10ef", A: "#e718e7", Z: "#c936c9", L: "#ae51ae", K: "#a05fa0", F: "#986798", Q: "#926d92", I: "#8a758a", W: "#8a758a", V: "#857a85", D: "#778877", X: "#758a75", H: "#758a75", R: "#6f906f", B: "#49b649", T: "#47b847", S: "#36c936", C: "#23dc23", Y: "#21de21", G: "#00ff00", N: "#1be41b", P: "#00ff00" },
  strand: { V: "#ffff00", I: "#ecec13", Y: "#d3d32c", F: "#c2c23d", W: "#c0c03f", L: "#b2b24d", T: "#9d9d62", C: "#9d9d62", Q: "#8c8c73", M: "#82827d", X: "#797986", R: "#6b6b94", N: "#64649b", H: "#60609f", A: "#5858a7", S: "#4949b6", G: "#4949b6", Z: "#4747b8", K: "#4747b8", B: "#4343bc", P: "#2323dc", D: "#2121de", E: "#0c0cf3" },
  turn: { N: "#ff0000", G: "#ff0000", P: "#f60009", B: "#f30c0c", D: "#e81717", S: "#e81717", C: "#a85757", Y: "#9d6262", K: "#7e8181", X: "#7c8383", Q: "#778888", W: "#738c8c", T: "#738c8c", R: "#708f8f", H: "#708f8f", Z: "#5ba4a4", E: "#3fc0c0", A: "#2cd3d3", F: "#1ee1e1", M: "#1ee1e1", L: "#1ce3e3", V: "#07f8f8", I: "#00ffff" },
  buried: { C: "#0000ff", I: "#0054ab", V: "#005fa0", L: "#007b84", F: "#008778", M: "#009768", G: "#009d62", A: "#00a35c", W: "#00a857", X: "#00b649", S: "#00d52a", H: "#00d52a", T: "#00db24", P: "#00e619", Y: "#00e619", N: "#00eb14", B: "#00eb14", D: "#00eb14", Q: "#00f10e", Z: "#00f10e", E: "#00f10e", R: "#00fc03", K: "#00ff00" },
  nucleotide: { A: "#64f73f", C: "#ffb340", G: "#eb413c", T: "#3c88ee", U: "#3c88ee" },
  purine: { A: "#ff83fa", G: "#ff83fa", R: "#ff83fa", C: "#40e0d0", T: "#40e0d0", U: "#40e0d0", Y: "#40e0d0" },
};
const CLUSTAL = { BLUE: "#80a0f0", RED: "#f01505", MAGENTA: "#c048c0", GREEN: "#15c015", PINK: "#f08080", ORANGE: "#f09048", YELLOW: "#c0c000", CYAN: "#15a4a4" };
// Clustal X colouring: a residue is coloured when its column meets the thresholds (Jalview's description of the Clustal X rules)
function clustalColour(res, pc) {                           // pc: percentage of each residue among the residues present in the column
  const sum = s => d3.sum(s, ch => pc[ch] || 0); const any = (s, t) => s.split("").some(ch => (pc[ch] || 0) > t);
  const hyd = sum("WLVIMAFCYHP") > 60;
  switch (res) {
    case "A": case "I": case "L": case "M": case "F": case "W": case "V": return hyd ? CLUSTAL.BLUE : null;
    case "C": return any("C", 85) ? CLUSTAL.PINK : hyd ? CLUSTAL.BLUE : null;
    case "K": case "R": return sum("KR") > 60 || any("KRQ", 85) ? CLUSTAL.RED : null;
    case "E": return sum("KR") > 60 || sum("QE") > 50 || sum("ED") > 50 || any("EQD", 85) ? CLUSTAL.MAGENTA : null;
    case "D": return sum("KR") > 60 || any("DEN", 85) || sum("ED") > 50 ? CLUSTAL.MAGENTA : null;
    case "N": return (pc.N || 0) > 50 || any("ND", 85) ? CLUSTAL.GREEN : null;
    case "Q": return sum("KR") > 60 || sum("QE") > 50 || any("QTKR", 85) ? CLUSTAL.GREEN : null;
    case "S": case "T": return hyd || sum("TS") > 50 || any("ST", 85) ? CLUSTAL.GREEN : null;
    case "G": return CLUSTAL.ORANGE;
    case "P": return CLUSTAL.YELLOW;
    case "H": case "Y": return hyd || any("WYACPQFHILMV", 85) ? CLUSTAL.CYAN : null;
    default: return null;
  }
}
const SCHEME_LIST = {
  aa: [["codes", "site codes"], ["none", "plain letters"], ["clustal", "Clustal X"], ["zappo", "Zappo"], ["taylor", "Taylor"], ["hydrophobic", "hydrophobicity"], ["helix", "helix propensity"], ["strand", "strand propensity"], ["turn", "turn propensity"], ["buried", "buried index"]],
  dna: [["codes", "site codes"], ["none", "plain letters"], ["nucleotide", "nucleotide"], ["purine", "purine / pyrimidine"]],
};

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
function drawCophylo(host, D, entry, k, cell) {
  const { g, li } = cell; const exp = parseNewick(entry.expected[k]), obs = parseNewick(entry.observed[k]);
  const name = x => D.taxa.names[+x]; const partner = D.num("partner", g, li), bm = D.num("best_match", g, li), bmId = D.num("best_match_identity", g, li); const comps = new Set(D.comps(g, li));
  // ladderise the expected tree, then rotate the observed tree to follow it; midpoint roots come from the tool
  const ladder = n => { n.children.forEach(ladder); n.children.sort((a, b) => tipsOf(a).length - tipsOf(b).length); }; ladder(exp);
  const rank = new Map(); tipsOf(exp).forEach((t, i) => rank.set(t.name, i)); orderLike(obs, rank);
  splitsOf(exp); splitsOf(obs); const expSplits = new Set(), obsSplits = new Set();
  const collect = (n, out) => { if (n.children.length && n.split && n.set.size > 1 && n.set.size < tipsOf(exp).length - 1) out.add(n.split); n.children.forEach(c => collect(c, out)); }; collect(exp, expSplits); collect(obs, obsSplits);
  const labelOf = t => { const i = +t; return name(t) + (i === bm && !comps.has(i) ? "*" : ""); };          // * : the best match, brought in from outside the comparators
  const meas = document.createElement("canvas").getContext("2d"); meas.font = `bold 11.5px ${cssVar("--sans")}`; const LAB = Math.ceil(d3.max(tipsOf(exp), t => meas.measureText(labelOf(t.name)).width) + 10);
  const nTip = tipsOf(exp).length; const RH = 20, GAP = 44, W = Math.min(820, Math.max(560, 2 * LAB + GAP + 240)), TW = (W - 2 * LAB - GAP - 28) / 2, top = 10, H = top + nTip * RH + 40;
  const maxD = Math.max(depthOf(exp), depthOf(obs), 1e-6); const sx = d3.scaleLinear().domain([0, maxD]).range([0, TW]);
  assignY(exp, { y: 0 }); assignY(obs, { y: 0 });
  const svg = d3.create("svg").attr("width", W).attr("height", H).attr("class", "tangle");
  const yOf = n => top + n.y * RH + RH / 2;
  const focalCol = D.colorOf(D.state[g * D.nloc + li]) || cssVar("--ink"); const partnerCol = D.colorOf(ST.BLUE + 1) || cssVar("--ink");
  const draw = (root, x0, dir, otherSplits) => {                  // dir +1 draws left to right, -1 right to left
    const gg = svg.append("g"); const walk = n => { for (const c of n.children) { const xa = x0 + dir * sx(n.depth), xb = x0 + dir * sx(c.depth);
      const conflict = c.children.length && c.split && c.set.size > 1 && c.set.size < nTip - 1 && !otherSplits.has(c.split);
      gg.append("path").attr("d", `M${xa},${yOf(n)}V${yOf(c)}H${xb}`).attr("fill", "none").attr("stroke", conflict ? cssVar("--ink") : cssVar("--muted")).attr("stroke-width", conflict ? 2.6 : 1.2).append("title").text(conflict ? "this split is not in the other tree" : `branch length ${fmt1(c.length)}%`);
      if (!c.children.length) { const isF = +c.name === g, isBM = +c.name === bm, isP = +c.name === partner; gg.append("circle").attr("cx", xb).attr("cy", yOf(c)).attr("r", 3).attr("fill", isF ? focalCol : cssVar("--surface")).attr("stroke", isF ? focalCol : isP ? partnerCol : cssVar("--muted")).attr("stroke-width", 1.2);
        const lab = gg.append("text").attr("x", xb + dir * 6).attr("y", yOf(c) + 4).attr("text-anchor", dir > 0 ? "start" : "end").attr("font-size", 11.5).attr("font-weight", isF || isBM ? 700 : 400).attr("fill", isF ? focalCol : isP ? partnerCol : cssVar("--ink")).text(labelOf(c.name));
        lab.append("title").text(`${name(c.name)}: ${fmt1(c.depth)}% ${D.format >= 3 ? "corrected divergence" : "divergence"} from the root of this tree`); }
      walk(c); } }; walk(root); return gg; };
  draw(exp, 8, 1, obsSplits); draw(obs, W - 8, -1, expSplits);
  const obsTips = new Map(); tipsOf(obs).forEach(t => obsTips.set(t.name, t));
  for (const t of tipsOf(exp)) { const o = obsTips.get(t.name); if (!o) continue;      // the connector runs through the gap between the two label columns
    const x1 = 8 + sx(t.depth) + 6 + meas.measureText(labelOf(t.name)).width + 4, x2 = W - 8 - sx(o.depth) - 6 - meas.measureText(labelOf(t.name)).width - 4;
    svg.append("line").attr("x1", x1).attr("x2", x2).attr("y1", yOf(t)).attr("y2", yOf(o)).attr("stroke", yOf(t) === yOf(o) ? cssVar("--grid") : cssVar("--floor")).attr("stroke-dasharray", "2 3"); }
  const bar = Math.pow(10, Math.floor(Math.log10(maxD))) * (maxD / Math.pow(10, Math.floor(Math.log10(maxD))) >= 5 ? 2 : 1);
  const by = H - 12; svg.append("line").attr("x1", 8).attr("x2", 8 + sx(bar)).attr("y1", by).attr("y2", by).attr("stroke", cssVar("--ink")).attr("stroke-width", 1.5);
  svg.append("text").attr("x", 12 + sx(bar)).attr("y", by + 4).attr("font-size", 10.5).attr("fill", cssVar("--muted")).text(`${bar}% ${D.format >= 3 ? "corrected divergence" : "divergence"}`);
  svg.append("text").attr("x", W - 36).attr("y", by + 4).attr("text-anchor", "end").attr("font-size", 10.5).attr("fill", cssVar("--muted")).text("expected  ·  observed at this locus");
  const q = svg.append("g").attr("transform", `translate(${W - 20},${by})`).style("cursor", "help");
  q.append("circle").attr("r", 8).attr("fill", "none").attr("stroke", cssVar("--muted")).attr("stroke-width", 1.2); q.append("text").attr("y", 4).attr("text-anchor", "middle").attr("font-size", 11).attr("font-weight", 700).attr("fill", cssVar("--muted")).text("?");
  const explain = el("div", { style: "max-width: 380px" },
    el("div", null, el("b", null, "Two neighbour-joining trees"), " over the focal sequence (filled dot, in its mask colour), its comparators, the partner of a too-similar pair (in the too-similar colour), and the best match, the taxon closest to the focal by per-site divergence compared on the same sites (in bold; marked * when it is not one of the comparators" + (bmId !== null ? `; identity ${fmt1(bmId)}%` : "") + ")."),
    el("div", { style: "margin-top: 4px" }, el("b", null, "Left:"), (D.format >= 3 ? " the distances the test expects for every pair: each pair's genome-wide divergence, F81-corrected with this locus's plateau and scaled by its rate. " : " the distances the test expects for every pair, each pair's genome-wide divergence through this locus's expectation. "), el("b", null, "Right:"), (D.format >= 3 ? " the divergence observed at this locus over shared sites, corrected the same way (capped at 0.98 of the plateau). Both midpoint-rooted, one scale in corrected divergence; negative branches drawn as zero." : " the mismatch fractions observed at this locus over shared sites. Both midpoint-rooted, one scale; negative branches drawn as zero.")),
    el("div", { style: "margin-top: 4px" }, "Thick branches are splits absent from the other tree; dotted lines join the same taxon on both sides."),
    el("div", { style: "margin-top: 4px", class: "more" }, "A long focal branch with its neighbours unchanged is a damaged copy of the right sequence. A focal that moves across to its best match points to contamination or a paralogue. A whole clade pushed away from the rest on one long branch is a paralogue shared by the clade. A best match on its own long branch is merely the least different taxon, not a relative."));
  q.on("mouseenter", ev => showHover(explain, ev.clientX, ev.clientY)).on("mouseleave", hideHover);
  return svg.node();
}

// ------------------------------------------------------------------ the alignment panel (full width, bottom)
const Sequence = (() => {
  let cur = null, pps = null, x0 = 0, scheme = "codes", dots = true, drawFn = null, showTrees = false;
  const RULER_H = 18, CONS_H = 26, CODE_H = 9, ROW_H = 18, NAME_W = 190;
  function open(g, li) {
    cur = { g, li }; pps = null; x0 = 0;
    if (!App.panels.sequence) {
      const P = makePanel("sequence", [], () => closePanel("sequence")); App.panels.sequence = { el: P.el, P, render, redraw: () => drawFn && drawFn() };
      App.bottomEl.append(P.el); App.bottomEl.style.display = ""; App.splitEl.style.display = ""; if (!App.bottomEl.style.flexBasis) App.bottomEl.style.flexBasis = "42%";
    }
    App.sel.cells.clear(); App.sel.cells.add(g * App.D.nloc + li); Overview.redraw(); render(); writeHash();
  }
  const isGap = b => b === 45 || b === 63 || b === 46 || b === 42 || b === 32;      // - ? . * space
  function tint(hex) { return mixHex(hex, cssVar("--surface"), isDark() ? 0.5 : 0.62); }
  function inkFor(hex) { return luminance(hex) < 0.35 ? "#ffffff" : "#0b0b0b"; }
  async function render() {
    if (!cur || !App.panels.sequence) return;
    const D = App.D; const { g, li } = cur; const P = App.panels.sequence.P; const i = g * D.nloc + li; const dtype = D.man.options.datatype;
    P.setCrumb([D.man.dataset || D.man.input.file, el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("taxon", g); } }, D.taxa.names[g]), el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("locus", li); } }, D.loci.names[li])]);
    P.body.replaceChildren(el("div", { class: "loading" }, "loading the sequences…"));
    const [C, S] = await Promise.all([D.codes(li), D.seq(li)]);
    if (!cur || cur.g !== g || cur.li !== li) return;
    const comps = D.comps(g, li); const partner = D.num("partner", g, li); const sd = D.cols.side[i];
    const rows = [{ g, role: "focal" }]; if (partner !== null && partner >= 0) rows.push({ g: partner, role: "partner" }); comps.forEach((c, r) => { if (c !== partner) rows.push({ g: c, role: `#${r + 1}` }); });
    const seqRow = gg => { if (!S) return null; for (let k = 0; k < S.taxa.length; k++) if (S.taxa[k] === gg) return S.blk.subarray(k * S.len, (k + 1) * S.len); return null; };
    let codeRow = null; if (C) for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) codeRow = C.blk.subarray(k * C.len, (k + 1) * C.len);
    const len = D.loci.length[li]; for (const r of rows) r.bytes = seqRow(r.g);
    const focal = rows[0].bytes;
    // column composition over the displayed rows: conservation (1 - normalised Shannon entropy) and effective residues
    const alphaN = dtype === "dna" ? 4 : 20; const cons = new Float32Array(len), neff = new Float32Array(len), npres = new Uint8Array(len);
    const colCounts = new Array(len);
    for (let s = 0; s < len; s++) {
      const cnt = {}; let tot = 0; for (const r of rows) { if (!r.bytes) continue; const b = r.bytes[s]; if (isGap(b)) continue; const ch = String.fromCharCode(b).toUpperCase(); cnt[ch] = (cnt[ch] || 0) + 1; tot++; }
      colCounts[s] = cnt; npres[s] = tot; if (tot < 2) { cons[s] = NaN; neff[s] = tot; continue; }
      let H = 0; for (const k in cnt) { const p = cnt[k] / tot; H -= p * Math.log(p); } const smax = Math.min(tot, alphaN); neff[s] = Math.exp(H); cons[s] = smax > 1 ? 1 - H / Math.log(smax) : 1;
    }
    // header: the cell's numbers in words, and the controls
    const head = el("div", { class: "seq-head" });
    const cl = D.num("closest", g, li), co = D.num("closest_overlap", g, li), ci = D.num("closest_identity", g, li), ce = D.num("closest_expected", g, li);
    const parts = [el("span", { class: "state" }, sw(D.colorOf(D.state[i])), D.labelOf(D.state[i]))];
    if (cl !== null && cl >= 0) parts.push(el("span", null, ` · against its closest comparator ${D.taxa.names[cl]}: `, el("span", { class: "num" }, `${fmtInt(Math.round((100 - ci) / 100 * co))} of ${fmtInt(co)} shared sites differ, ${fmtInt(Math.round((100 - ce) / 100 * co))} expected`)));
    if (partner !== null && partner >= 0) { const lam = D.num("lam", g, li), lp = D.num("lam_p", g, li); parts.push(el("span", null, ` · too similar to ${D.taxa.names[partner]}${sd & 8 ? ", identified as the recipient" : ", direction not identified"}${lam !== null ? ` (direction LRT |2 ln LR| = ${fmt1(lam)}, p = ${lp < 1e-4 ? lp.toExponential(1) : fmt3(lp)})` : ""}${(() => { const rf = D.num("row_factor", g, li); return rf !== null && rf < 1 ? `, row factor ${fmt3(rf)}` : ""; })()}`)); }
    const xm = D.num("excess_m", g, li); if (xm !== null) parts.push(el("span", null, ` · excess ${fmt2(xm)} m`));
    const schemeSel = el("select", { onchange: e => { scheme = e.target.value; draw(); } }, ...SCHEME_LIST[dtype === "dna" ? "dna" : "aa"].map(([v, t]) => el("option", { value: v, selected: v === scheme ? "" : null }, t)));
    const dotsBox = el("input", { type: "checkbox", checked: dots ? "" : null, onchange: e => { dots = e.target.checked; draw(); } });
    const treesBtn = el("button", { title: "the expected and the observed distances among these sequences as two neighbour-joining trees, face to face", onclick: () => toggleTrees() }, "trees");
    const ctl = el("span", { class: "seq-ctl" }, el("label", null, "colours ", schemeSel), el("label", null, dotsBox, " dots where a comparator matches the focal residue"),
      el("button", { onclick: () => { pps = fitPps(); x0 = 0; draw(); } }, "whole locus"), el("button", { onclick: () => { pps = 12; draw(); } }, "letters"), treesBtn);
    head.append(el("div", { class: "seq-line" }, ...parts), ctl);
    const body = el("div", { class: "seq-body" }); const names = el("div", { class: "seq-names" }); const cwrap = el("div", { class: "seq-canvas-wrap" }); const cv = el("canvas"); cwrap.append(cv); body.append(names, cwrap);
    let treesCard = null;
    const toggleTrees = async () => {
      if (treesCard) { treesCard.remove(); treesCard = null; treesBtn.classList.remove("primary"); showTrees = false; return; }
      showTrees = true; treesBtn.classList.add("primary");
      treesCard = el("div", { class: "trees-card" }, el("div", { class: "loading", style: "padding: 10px" }, "loading the trees…")); cwrap.append(treesCard);
      const entry = await D.trees(li); if (!treesCard) return; treesCard.replaceChildren();
      let k = -1; if (entry) for (let q = 0; q < entry.cells.length; q++) if (entry.cells[q] === g) { k = q; break; }
      if (k < 0) { treesCard.append(el("div", { class: "note", style: "padding: 8px" }, "trees were not carried for this cell")); return; }
      treesCard.append(el("button", { class: "close", title: "close", onclick: toggleTrees }, "✕"), drawCophylo(treesCard, D, entry, k, { g, li }));
    };
    if (showTrees) setTimeout(toggleTrees, 0);
    names.append(el("div", { class: "trk", style: `height:${RULER_H}px` }, "position"), el("div", { class: "trk", style: `height:${CONS_H}px`, title: "1 − normalised Shannon entropy of the column over the rows shown" }, "conservation"), el("div", { class: "trk", style: `height:${CODE_H + 4}px` }, "site codes"));
    for (const r of rows) names.append(el("div", { class: "seqnm" + (r.role === "focal" ? " focal" : ""), style: `height:${ROW_H}px`, title: `${D.taxa.names[r.g]} (${r.role === "#1" ? "comparator 1" : r.role})` }, el("span", { class: "nm" }, D.taxa.names[r.g]), el("span", { class: "role" }, r.role === "focal" ? "" : r.role)));
    const foot = el("div", { class: "seq-foot" }, el("span", null, `${fmtInt(len)} sites · pinch or ctrl+wheel: zoom · two-finger scroll or drag: move · letters appear from 7 px per site`), el("span", { class: "k", style: "margin-left:auto" }, S ? "" : "sequences were not carried for this locus"));
    P.body.replaceChildren(head, body, foot);
    if (!S && !codeRow) { cwrap.append(el("div", { class: "note", style: "padding: 10px" }, "neither sequences nor site codes were carried for this cell in this report")); return; }
    const dpr = devicePixelRatio || 1;
    const fitPps = () => Math.max(0.05, cwrap.getBoundingClientRect().width / len);
    const topOf = ri => RULER_H + CONS_H + CODE_H + 4 + ri * ROW_H;
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
      // conservation track: bars, averaged over the sites a pixel covers
      const cy0 = RULER_H + 2, ch = CONS_H - 5; c.fillStyle = cssVar("--grey-3");
      if (letters) { for (let s = s0; s <= s1; s++) { const v = cons[s]; if (!Number.isFinite(v)) continue; c.fillRect(sx(s) + 1, cy0 + ch * (1 - v), Math.max(pps - 2, 1), ch * v); } }
      else { for (let x = 0; x < W; x++) { const a = Math.floor(x0 + x / pps), b = Math.max(a, Math.ceil(x0 + (x + 1) / pps) - 1); let sum = 0, k = 0; for (let s = a; s <= Math.min(b, len - 1); s++) if (Number.isFinite(cons[s])) { sum += cons[s]; k++; } if (k) { const v = sum / k; c.fillRect(x, cy0 + ch * (1 - v), 1, ch * v); } } }
      c.strokeStyle = cssVar("--grid"); c.beginPath(); c.moveTo(0, RULER_H + CONS_H + 0.5); c.lineTo(W, RULER_H + CONS_H + 0.5); c.stroke();
      // code track: the site codes at full strength
      const ty = RULER_H + CONS_H + 2;
      if (codeRow) {
        if (letters) { for (let s = s0; s <= s1; s++) { const h = codeHex(D, codeRow[s]); if (h) { c.fillStyle = h; c.fillRect(sx(s), ty, pps, CODE_H); } } }
        else { for (let x = 0; x < W; x++) { const a = Math.floor(x0 + x / pps), b = Math.max(a, Math.ceil(x0 + (x + 1) / pps) - 1); let best = -1, code = 255; for (let s = a; s <= Math.min(b, len - 1); s++) { const sv = codeSeverity(D, codeRow[s]); if (sv > best) { best = sv; code = codeRow[s]; } } const h = codeHex(D, code); if (h) { c.fillStyle = h; c.fillRect(x, ty, 1, CODE_H); } } }
      }
      // sequence rows
      const resColour = (ch, s) => {                           // the residue colour scheme, or null
        if (scheme === "codes" || scheme === "none") return null;
        if (scheme === "clustal") { const cnt = colCounts[s]; const tot = npres[s] || 1; const pc = {}; for (const k in cnt) pc[k] = 100 * cnt[k] / tot; return clustalColour(ch, pc); }
        return (RESIDUE_SCHEMES[scheme] || {})[ch] || null;
      };
      const grey1 = cssVar("--grey-1"), ink = cssVar("--ink"), faint = cssVar("--faint"), inkRgb = hexToRgb(cssVar("--muted"));
      c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; c.textBaseline = "middle"; c.textAlign = "center";
      rows.forEach((r, ri) => {
        const y = topOf(ri); const bytes = r.bytes; const isFocal = ri === 0, isPartner = r.role === "partner";
        if (!bytes) { c.fillStyle = faint; c.font = `11px ${cssVar("--sans")}`; c.textAlign = "left"; c.fillText("row not carried in this report", 4, y + ROW_H / 2); c.textAlign = "center"; c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; return; }
        if (letters) {
          for (let s = s0; s <= s1; s++) {
            const b = bytes[s]; if (isGap(b)) continue; const ch = String.fromCharCode(b).toUpperCase(); const x = sx(s);
            let bg = null;
            if (isFocal && codeRow && scheme === "codes") { const h = codeHex(D, codeRow[s]); if (h && codeRow[s] !== 254) bg = tint(h); }
            else if (isPartner && codeRow && codeRow[s] >= 128 && codeRow[s] < 254) bg = tint(codeHex(D, codeRow[s]));   // the matching residue of the too-similar partner
            const same = !isFocal && focal && !isGap(focal[s]) && focal[s] === b;
            const rc = (isFocal || !same || !dots) ? resColour(ch, s) : null; if (rc && (!bg || !isFocal)) bg = rc;
            if (bg) { c.fillStyle = bg; c.fillRect(x, y + 1, pps, ROW_H - 2); }
            if (!isFocal && same && dots) { c.fillStyle = faint; c.fillText("·", x + pps / 2, y + ROW_H / 2); }
            else { c.fillStyle = bg ? inkFor(bg) : (isFocal || same ? ink : ink); if (!isFocal && !same && !bg) c.fillStyle = ink; c.fillText(ch, x + pps / 2, y + ROW_H / 2); }
          }
        } else {
          // block zoom: the focal row wears its codes (tinted); every other row shows presence and, per pixel, the share of its sites that differ from the focal
          for (let x = 0; x < W; x++) {
            const a = Math.floor(x0 + x / pps), b = Math.min(len - 1, Math.max(a, Math.ceil(x0 + (x + 1) / pps) - 1)); if (a > b) break;
            let pres = 0, diff = 0, best = -1, code = 255, blue = false;
            for (let s = a; s <= b; s++) { const bb = bytes[s]; if (isGap(bb)) continue; pres++; if (!isFocal && focal && !isGap(focal[s]) && focal[s] !== bb) diff++;
              if (codeRow) { const cc = codeRow[s]; if (isFocal) { const sv = codeSeverity(D, cc); if (sv > best) { best = sv; code = cc; } } else if (isPartner && cc >= 128 && cc < 254) blue = true; } }
            if (!pres) continue;
            if (isFocal) { const h = codeRow && scheme === "codes" ? codeHex(D, code) : null; c.fillStyle = h && code !== 254 ? tint(h) : grey1; c.fillRect(x, y + 2, 1, ROW_H - 4); }
            else { c.fillStyle = blue ? tint(D.blues[1]) : grey1; c.fillRect(x, y + 4, 1, ROW_H - 8); if (diff) { c.fillStyle = `rgba(${inkRgb[0]},${inkRgb[1]},${inkRgb[2]},${0.15 + 0.75 * diff / pres})`; c.fillRect(x, y + 3, 1, ROW_H - 6); } }
          }
        }
      });
    };
    drawFn = draw; draw();
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
      if (codeRow) t.append(el("tr", null, el("td", { class: "k" }, "code"), el("td", null, codeLabel(D, codeRow[s]))));
      if (npres[s] >= 2) t.append(el("tr", null, el("td", { class: "k" }, "conservation"), el("td", null, `${fmt2(cons[s])} · ${fmt1(neff[s])} effective residues among the ${npres[s]} rows with one`)));
      for (const rr of rows) if (rr.bytes) t.append(el("tr", null, el("td", { class: "k" }, D.taxa.names[rr.g]), el("td", { style: "font-family: var(--mono)" }, String.fromCharCode(rr.bytes[s]))));
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
    const out = [[0, 0]]; for (let lo = 1; lo <= kmax; lo += w) out.push([lo, Math.min(lo + w - 1, kmax)]); return out;
  }
  function facet(by, root, width, height) {
    const D = App.D; const N = by === "locus" ? D.nloc : D.n; const cnt = by === "locus" ? D.locCount : D.taxCount; const names = by === "locus" ? D.loci.names : D.taxa.names;
    const units = d3.range(N).filter(u => cnt.present[u] > 0); const kof = units.map(u => cnt.masked[u]);
    const B = bins(kof); const heights = B.map(([lo, hi]) => units.filter((u, j) => kof[j] >= lo && kof[j] <= hi));
    const rotated = B.length > 14; const W = Math.max(300, width), H = Math.max(170, height), m = { l: 44, r: 12, t: 22, b: rotated ? 76 : 40 };
    const xs = d3.scaleBand().domain(B.map((b, i) => i)).range([m.l, W - m.r]).padding(0.15);
    const ymax = d3.max(heights, h => h.length) || 1; const ys = d3.scaleSymlog().constant(2).domain([0, ymax * 1.15]).range([H - m.b, m.t]);
    const box = el("div", { class: "plot" }, el("h3", null, `${by === "locus" ? "loci" : "taxa"} by the number of masked sequences they hold`, el("span", { class: "sub" }, `${fmtInt(units.filter((u, j) => kof[j] > 0).length)} of ${fmtInt(units.length)}`)));
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).tickValues([0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].filter(v => v <= ymax * 1.15)).tickFormat(d3.format(",")));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).tickFormat(i => { const [lo, hi] = B[i]; return lo === hi ? String(lo) : `${lo}–${hi}`; })).selectAll("text").attr("transform", rotated ? "rotate(-60)" : null).attr("text-anchor", rotated ? "end" : "middle").attr("dx", rotated ? "-0.5em" : null).attr("dy", rotated ? "0.1em" : null);
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 4).attr("text-anchor", "middle").text(`sequences masked in the ${by}`);
    svg.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text(by === "locus" ? "loci" : "taxa");
    B.forEach(([lo, hi], i) => {
      const us = heights[i]; if (!us.length) return;
      const comp = new Map();
      if (lo === 0) comp.set(ST.GREY, 1);
      else { for (const u of us) for (let k = 0; k < (by === "locus" ? D.n : D.nloc); k++) { const idx = by === "locus" ? k * D.nloc + u : u * D.nloc + k; const code = D.state[idx]; if (code >= ST.RED) comp.set(code, (comp.get(code) || 0) + 1); } }
      const tot = d3.sum(comp.values()) || 1; let y0 = ys(0); const codes = Array.from(comp.keys()).sort(d3.ascending);
      for (const code of codes) { const h = (ys(0) - ys(us.length)) * comp.get(code) / tot; svg.append("rect").attr("x", xs(i)).attr("y", y0 - h).attr("width", xs.bandwidth()).attr("height", h).attr("fill", D.colorOf(code)).attr("stroke", cssVar("--surface")).attr("stroke-width", 0.8).style("cursor", "pointer")
        .on("mousemove", ev => showHover(el("div", null, el("b", null, `${fmtInt(us.length)} ${by === "locus" ? "loci" : "taxa"} with ${lo === hi ? lo : `${lo} to ${hi}`} masked sequence${hi !== 1 ? "s" : ""}`), lo === 0 ? null : el("div", null, `${D.labelOf(code)}: ${fmtInt(comp.get(code))} sequences (${pct(comp.get(code), tot)})`), el("div", { class: "more" }, us.slice(0, 8).map(u => names[u]).join(", ") + (us.length > 8 ? ", …" : "") + " · click to highlight them in the matrix")), ev.clientX, ev.clientY))
        .on("mouseleave", hideHover).on("click", () => { if (by === "locus") { App.sel.cols = new Set(us); App.sel.rows.clear(); } else { App.sel.rows = new Set(us); App.sel.cols.clear(); } Overview.setView("matrix"); }); y0 -= h; }
      svg.append("text").attr("x", xs(i) + xs.bandwidth() / 2).attr("y", ys(us.length) - 3).attr("text-anchor", "middle").style("fill", cssVar("--ink")).text(fmtInt(us.length));
    });
    box.append(svg.node()); svg.append("title").text("zero is its own bar; the other bars are coloured by the kinds of mask they hold. The y axis is linear near zero and logarithmic above. Click a bar to highlight its units in the matrix.");
    return box;
  }
  return { render };
})();

// ------------------------------------------------------------------ the table view (virtualised)
const Table = (() => {
  let rows = [], keys = [{ k: "excess_m", dir: -1 }], maskedOnly = true, root = null;          // keys[0] sorts first, keys[1] breaks its ties
  const COLS = [["taxon", "taxon"], ["locus", "locus"], ["cls", "class"], ["to", "moved towards"], ["bm", "best match"], ["n_fail", "fails"], ["n_similar", "similar"], ["obs", "observed %"], ["exp", "expected %"], ["excess_m", "excess (m)"], ["bm_id", "best match %"], ["away", "away %"], ["lam", "|2 ln LR|"]];
  function build() {
    const D = App.D; rows = [];
    for (let g = 0; g < D.n; g++) for (let li = 0; li < D.nloc; li++) { const i = g * D.nloc + li; if (D.cols.n_comp[i] === 0) continue; if (maskedOnly && !(D.cols.side[i] & 3)) continue;
      const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li); const bm = D.num("best_match", g, li), pt = D.num("partner", g, li);
      const cc = D.cols.class[i]; const to = D.movedTo(g, li); const aw = D.num("away_frac", g, li);
      rows.push({ g, li, taxon: D.taxa.names[g], locus: D.loci.names[li], cls: D.CLASS[cc] || "", to: to >= 0 ? D.taxa.names[to] : "", n_fail: D.cols.n_fail[i], n_similar: D.cols.n_similar[i], obs: mi === null ? null : 100 - mi, exp: me === null ? null : 100 - me, excess_m: D.num("excess_m", g, li), bm_id: D.num("best_match_identity", g, li), away: aw === null ? null : 100 * aw, lam: D.num("lam", g, li), bm: bm !== null && bm >= 0 ? D.taxa.names[bm] : "", code: D.state[i] }); }
    sort();
  }
  const cmp = (a, b, k, dir) => { const x = a[k], y = b[k]; const ex = x === null || x === undefined || x === "", ey = y === null || y === undefined || y === ""; if (ex && ey) return 0; if (ex) return 1; if (ey) return -1; return (typeof x === "string" ? x.localeCompare(y) : x - y) * dir; };
  function sort() { rows.sort((a, b) => { for (const { k, dir } of keys) { const c = cmp(a, b, k, dir); if (c) return c; } return 0; }); }
  function render(r) {
    root = r; root.replaceChildren(); build(); const D = App.D;
    const mark = k => { const j = keys.findIndex(x => x.k === k); return j < 0 ? "" : (keys[j].dir < 0 ? " ▾" : " ▴") + (keys.length > 1 ? (j ? "²" : "¹") : ""); };
    const click = (k, ev) => { const j = keys.findIndex(x => x.k === k); const natural = typeof rows.find(r => r[k] !== null && r[k] !== undefined && r[k] !== "")?.[k] === "string" ? 1 : -1;
      if (ev.shiftKey) { if (j === 0 && keys.length === 1) return; if (j === 1) keys[1].dir = -keys[1].dir; else keys = [keys[0], { k, dir: natural }]; }     // shift-click: the secondary key
      else if (j === 0) keys[0].dir = -keys[0].dir; else keys = [{ k, dir: natural }, ...(keys[0].k === k ? [] : [keys[0]])].slice(0, 2);                 // click: the primary key; the old primary breaks ties
      sort(); paint(); head.replaceChildren(...COLS.map(([k2, t2]) => el("span", { onclick: e => click(k2, e) }, t2 + mark(k2)))); };
    const head = el("div", { class: "thead" }, ...COLS.map(([k, t]) => el("span", { onclick: e => click(k, e) }, t + mark(k))));
    const body = el("div", { class: "tbody" }); const inner = el("div"); body.append(inner);
    const foot = el("div", { class: "tfoot" }, el("label", null, el("input", { type: "checkbox", checked: maskedOnly ? "" : null, onchange: e => { maskedOnly = e.target.checked; build(); paint(); } }), " masked cells only"), el("span", { class: "n" }));
    root.append(head, body, foot);
    const RH = 22;
    const paint = () => { inner.style.height = rows.length * RH + "px"; inner.style.position = "relative"; inner.replaceChildren(); const top = body.scrollTop; const i0 = Math.max(0, Math.floor(top / RH) - 5), i1 = Math.min(rows.length, i0 + Math.ceil(body.clientHeight / RH) + 10);
      for (let i = i0; i < i1; i++) { const r = rows[i]; inner.append(el("div", { class: "trow", style: `position:absolute; top:${i * RH}px; left:0; right:0`, onclick: () => Sequence.open(r.g, r.li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(r.g * D.nloc + r.li); } },
        el("span", null, sw(D.colorOf(r.code)), r.taxon), el("span", null, r.locus), el("span", { class: r.cls ? "cls-" + r.cls : null }, r.cls), el("span", null, r.to), el("span", null, r.bm), el("span", null, r.n_fail), el("span", null, r.n_similar), el("span", null, r.obs === null ? "" : fmt2(r.obs)), el("span", null, r.exp === null ? "" : fmt2(r.exp)), el("span", null, r.excess_m === null ? "" : fmt2(r.excess_m)), el("span", null, r.bm_id === null ? "" : fmt1(r.bm_id)), el("span", null, r.away === null ? "" : fmt1(r.away)), el("span", null, r.lam === null ? "" : fmt1(r.lam)))); }
      $(".n", foot).textContent = `${fmtInt(rows.length)} rows · click a row for the alignment · click a header to sort, shift-click for the tie-breaking key`; };
    body.addEventListener("scroll", paint); paint();
  }
  function download() {
    const D = App.D; const keep = maskedOnly; maskedOnly = false; build(); maskedOnly = keep;
    const lines = ["taxon\tlocus\tstate\tclass\tmoved_towards\taway_frac\tn_fail\tn_similar\tobserved_divergence\texpected_divergence\texcess_m\tbest_match\tbest_match_identity\tpartner\tlrt"];
    for (const r of rows) { const pt = D.num("partner", r.g, r.li); const bm = D.num("best_match", r.g, r.li); lines.push([r.taxon, r.locus, D.labelOf(r.code), r.cls, r.to, r.away === null ? "" : (r.away / 100).toFixed(3), r.n_fail, r.n_similar, r.obs ?? "", r.exp ?? "", r.excess_m ?? "", bm !== null && bm >= 0 ? D.taxa.names[bm] : "", r.bm_id ?? "", pt !== null && pt >= 0 ? D.taxa.names[pt] : "", r.lam ?? ""].join("\t")); }
    const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/tab-separated-values" })); a.download = `${D.man.dataset || "report"}_cells.tsv`; a.click(); build();
  }
  return { render, download };
})();

// ------------------------------------------------------------------ URL hash: view, ordering, viewport, open panels
let hashTimer = null;
function writeHash() {
  clearTimeout(hashTimer); hashTimer = setTimeout(() => {
    const D = App.D; if (!D) return; const V = Overview.V; const p = new URLSearchParams();
    p.set("d", D.man.input.sha256.slice(0, 10)); p.set("view", App.view); p.set("o", `${App.order.taxaMode},${App.order.lociMode}`);
    p.set("v", [V.x0, V.y0, V.cw, V.ch].map(v => +v.toFixed(3)).join(","));
    const dc = Detail.cur(); if (dc && App.panels.detail) { p.set(dc.kind === "taxon" ? "t" : "l", dc.i); const ds = App.panels.detail.P.el.querySelector(".crumb select"); if (ds) p.set("dv", ds.value); }
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
  if (p.get("t")) Detail.open("taxon", +p.get("t")); else if (p.get("l")) Detail.open("locus", +p.get("l"));
  if (p.get("c")) { const [g, li] = p.get("c").split(",").map(Number); Sequence.open(g, li); }
}
const _render = Overview.render;
Overview.render = function () { if (App._viewFromHash) { App._viewFromHash = false; const V = Overview.V; const keep = { ...V }; _render(); Object.assign(V, keep); Overview.labelSpace(); Overview.redraw(); return; } _render(); };
// the alignment panel redraws when its row is resized
new ResizeObserver(debounce(() => { const p = App.panels.sequence; if (p && p.redraw) p.redraw(); }, 120)).observe(document.body);

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
