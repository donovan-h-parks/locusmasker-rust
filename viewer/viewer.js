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
function ramp(cols, n) {                                   // n steps over the 3 anchor colours, as plot_panels.ramp
  if (n <= 0) return []; if (n === 1) return [cols[0]]; if (n === 2) return [cols[0], cols[2]]; if (n === 3) return cols.slice();
  return d3.range(n).map(i => { const u = 2 * i / (n - 1), j = Math.min(Math.floor(u), 1); return mixHex(cols[j], cols[j + 1], u - j); });
}

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
  if (man.format_version > 1) throw new Error(`this report was written in bundle format ${man.format_version}; this viewer reads format 1. Please use the viewer that came with the tool version ${man.version}.`);
  const taxa = await json("taxa.json"), loci = await json("loci.json");
  const D = { zip, man, taxa, loci, n: taxa.names.length, nloc: loci.names.length, K: man.options.K_eff, M: man.options.M_eff, cols: {}, cache: { codes: new Map(), seq: new Map() } };
  const byteCols = ["n_comp", "n_fail", "n_similar", "side", "focal_sites"];
  const numCols = ["med_identity", "med_expected", "excess_m", "closest", "closest_overlap", "closest_identity", "closest_expected", "best_match", "best_match_identity", "partner", "lam", "lam_p"];
  for (const k of byteCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  D.dense = man.carried.numeric_columns === "dense";
  for (const k of numCols) D.cols[k] = (await npy(`cells/${k}.npy`)).data;
  D.comparators = (await npy("cells/comparators.npy")).data;
  D.closest = (await npy("taxa/closest.npy")).data;
  if (!D.dense) {
    const sl = (await npy("cells/sparse_locus.npy")).data, st = (await npy("cells/sparse_taxon.npy")).data;
    D.sparse = new Map(); for (let i = 0; i < sl.length; i++) D.sparse.set(st[i] * D.nloc + sl[i], i);
  }
  // one accessor: a numeric column of cell (g, li), or null when the bundle does not carry it
  D.num = (k, g, li) => {
    const arr = D.cols[k]; if (!arr) return null;
    if (D.dense) { const v = arr[g * D.nloc + li]; return Number.isNaN(v) ? null : v; }
    const r = D.sparse.get(g * D.nloc + li); if (r === undefined) return null; const v = arr[r]; return Number.isNaN(v) ? null : v;
  };
  D.comps = (g, li) => {                                      // the cell's comparators in rank order
    let base; if (D.dense) base = (g * D.nloc + li) * D.K; else { const r = D.sparse.get(g * D.nloc + li); if (r === undefined) return []; base = r * D.K; }
    const out = []; for (let r = 0; r < D.K; r++) { const c = D.comparators[base + r]; if (c >= 0) out.push(c); } return out;
  };
  D.codes = async li => {                                     // {cells: Int32Array, blk: Uint8Array, len}
    if (D.cache.codes.has(li)) return D.cache.codes.get(li);
    const p = (async () => { if (!zip.has(`codes/${li}.npy`)) return null; const b = await npy(`codes/${li}.npy`), c = await npy(`codes/${li}_cells.npy`); return { cells: c.data, blk: b.data, len: b.shape[1], rows: b.shape[0] }; })();
    D.cache.codes.set(li, p); return p;
  };
  D.seq = async li => {                                       // {taxa: Int32Array, blk: Uint8Array, len}
    if (D.cache.seq.has(li)) return D.cache.seq.get(li);
    const p = (async () => { if (!zip.has(`seq/${li}.npy`)) return null; const b = await npy(`seq/${li}.npy`), t = await npy(`seq/${li}_taxa.npy`); return { taxa: t.data, blk: b.data, len: b.shape[1], rows: b.shape[0] }; })();
    D.cache.seq.set(li, p); return p;
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
    else { const f = Math.min(nf[i], M - 1); st[i] = ST.GREY + f; sev[i] = 2 + f; }
  }
  D.state = st; D.sev = sev;
  // palette by state code
  const greys = ramp([cssVar("--grey-1"), cssVar("--grey-2"), cssVar("--grey-3")], M);
  const reds = ramp([cssVar("--red-1"), cssVar("--red-2"), cssVar("--red-3")], K - M + 1);
  const blues = [cssVar("--blue-1"), cssVar("--blue-2"), cssVar("--blue-3")];
  D.colorOf = code => code === ST.ABSENT ? null : code === ST.UNTESTED ? "hatch" : code < ST.RED ? greys[code - ST.GREY] : code < ST.BLUE ? reds[code - ST.RED - M] : blues[code - ST.BLUE - 1];
  D.labelOf = code => code === ST.ABSENT ? "absent" : code === ST.UNTESTED ? "present, not tested" : code < ST.RED ? `not masked · fails ${code - ST.GREY} of ${K}` : code < ST.BLUE ? `masked, too divergent · fails ${code - ST.RED} of ${K}` : `masked, too similar · ${code - ST.BLUE}${code - ST.BLUE === 3 ? "+" : ""} pair${code - ST.BLUE > 1 ? "s" : ""}`;
  D.greys = greys; D.reds = reds; D.blues = blues;
  // RGBA lookup by state code for the pixel buffer
  const surf = hexToRgb(cssVar("--surface")), hatch = hexToRgb(cssVar("--surface-3"));
  const lut = new Uint8ClampedArray(64 * 4);
  for (let c = 0; c < 64; c++) {
    let rgb = surf; const col = D.colorOf(c);
    if (col === "hatch") rgb = hatch; else if (col) rgb = hexToRgb(col);
    lut.set([rgb[0], rgb[1], rgb[2], 255], c * 4);
  }
  D.lut = lut;
  // per-unit counts
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
const App = { D: null, layout: "vertical", panels: {}, hover: null, sel: { rows: new Set(), cols: new Set(), cells: new Set() }, order: { taxa: null, loci: null, taxaMode: "input", lociMode: "input" }, hidePass: false, view: null };

function boot() {
  const app = $("#app");
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
window.LM = { start, App };                                   // the website hands the bundle over through LM.start(bytes)

// ------------------------------------------------------------------ shell: topbar, panels container, hover card
function buildShell() {
  const D = App.D, m = D.man; const app = $("#app"); app.replaceChildren();
  const c = m.counts;
  const theme = el("button", { title: "light / dark", onclick: () => { const r = document.documentElement; r.dataset.theme = r.dataset.theme === "dark" ? "light" : "dark"; deriveStates(D); renderAll(); } }, "◐");
  const layoutBtn = el("button", { title: "stack the panels vertically or side by side", onclick: () => setLayout(App.layout === "vertical" ? "horizontal" : "vertical") }, "⇅ layout");
  const about = el("details", { class: "meta" }, el("summary", null, `${m.tool} ${m.version} · ${m.date} · ${m.platform}`),
    el("div", null, `commit ${m.commit} · python ${m.python} · numpy ${m.numpy} · ${m.system}`),
    el("div", null, `input ${m.input.file} (${(m.input.bytes / 1e6).toFixed(1)} MB, sha256 ${m.input.sha256.slice(0, 12)}…)${m.input.partitions ? ", partitions " + m.input.partitions : ""}`),
    el("div", null, `options: datatype ${m.options.datatype}, n-comp ${m.options.n_comp} (K ${m.options.K_eff}, majority ${m.options.M_eff}), fit-overlap ${m.options.fit_overlap}, alpha ${m.options.alpha}, min-ratio ${m.options.min_ratio}${m.options.mask_similar_only ? ", similar side only" : ""}${m.options.mask_different_only ? ", divergent side only" : ""}`),
    el("div", null, `test: ${fmtInt(m.test.comparisons)} distinct comparisons, Bonferroni k ${m.test.bonferroni_k}, dataset dispersion ${m.test.phi_dataset}; carried: ${m.carried.sequences} sequences, codes for ${m.carried.codes} cells, best matches for ${m.carried.best_match_cells} cells`),
    el("div", null, `command: ${m.command}`));
  app.append(el("div", { class: "topbar" },
    el("h1", null, m.dataset || m.input.file), el("span", { class: "sub" }, `${fmtInt(c.masked)} of ${fmtInt(c.cells)} sequences masked (${pct(c.masked, c.cells)}) · ${fmtInt(c.taxa)} taxa × ${fmtInt(c.loci)} loci`),
    el("span", { class: "spacer" }), about, layoutBtn, theme));
  const panels = el("div", { class: "panels " + App.layout, id: "panels" }); app.append(panels);
  App.hover = el("div", { class: "hover" }); document.body.append(App.hover);
  App.panelsEl = panels;
  App.panels.overview = Overview.create(); panels.append(App.panels.overview.el);
}
function setLayout(mode) { App.layout = mode; App.panelsEl.className = "panels " + mode; writeHash(); Overview.redraw(); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); }
function renderAll() { Overview.render(); for (const k of ["detail", "sequence"]) if (App.panels[k]) App.panels[k].render(); }
function showHover(html, x, y) { const h = App.hover; h.innerHTML = ""; h.append(html); h.style.display = "block"; const r = h.getBoundingClientRect(); const px = x + 14 + r.width > innerWidth ? x - r.width - 10 : x + 14; const py = y + 14 + r.height > innerHeight ? y - r.height - 10 : y + 14; h.style.left = px + "px"; h.style.top = py + "px"; }
function hideHover() { App.hover.style.display = "none"; }
function sw(color) { const s = el("span", { class: "sw" }); if (color === "hatch") s.style.background = "var(--hatch)"; else if (color) s.style.background = color; else { s.style.border = "1px solid var(--spine)"; } return s; }

// a panel frame with breadcrumb and the collapse / maximise / close tools
function makePanel(key, crumbParts, onClose) {
  const body = el("div", { class: "panel-body" });
  const crumb = el("span", { class: "crumb" });
  const setCrumb = parts => { crumb.replaceChildren(); parts.forEach((p, i) => { if (i) crumb.append(el("span", { class: "sep" }, "›")); crumb.append(typeof p === "string" ? p : p); }); };
  setCrumb(crumbParts);
  const frame = el("div", { class: "panel", "data-key": key });
  const tools = el("div", { class: "tools" });
  const collapse = el("button", { title: "collapse", onclick: () => { frame.classList.toggle("collapsed"); collapse.textContent = frame.classList.contains("collapsed") ? "▸" : "▾"; } }, "▾");
  const maxi = el("button", { title: "maximise", onclick: () => { const on = frame.classList.toggle("maximised"); App.panelsEl.classList.toggle("has-max", on); for (const p of Object.values(App.panels)) if (p && p.redraw) p.redraw(); Overview.redraw(); } }, "⤢");
  tools.append(collapse, maxi);
  if (onClose) tools.append(el("button", { title: "close", onclick: onClose }, "✕"));
  frame.append(el("div", { class: "panel-head" }, crumb, tools), body);
  return { el: frame, body, setCrumb, tools };
}
function closePanel(key) { const p = App.panels[key]; if (!p) return; p.el.remove(); App.panels[key] = null; if (key === "detail" && App.panels.sequence) closePanel("sequence"); App.panelsEl.classList.remove("has-max"); writeHash(); Overview.redraw(); }

// ------------------------------------------------------------------ the overview panel
const Overview = (() => {
  let P, matrix, ctx, rowCanvas, colCanvas, strip, stripCtx, wrap, tabs = "matrix", viewEl = {}, summaryUser = null, summaryShown = true;
  function autoSummary() {                                     // the summary folds when the panel has no room for it and the matrix
    const h = P.body.getBoundingClientRect().height; const show = summaryUser === null ? h >= 640 : summaryUser;
    summaryShown = show; viewEl.summary.style.display = show ? "" : "none"; viewEl.sumBtn.textContent = show ? "summary ▾" : "summary ▸";
  }
  // view: first display row/col and cell size in CSS px
  const V = { x0: 0, y0: 0, cw: 4, ch: 4 };
  let hoverCells = null, dragging = null, lastMouse = null;
  function create() {
    P = makePanel("overview", [App.D.man.dataset || App.D.man.input.file, "overview"]);
    P.el.classList.add("overview");
    const D = App.D;
    // control row
    const taxaSel = el("select", { onchange: e => setOrder("taxa", e.target.value) }, ...[["input", "input order"], ["masked", "most masked first"], ["name", "name"], ["grouped", "grouped by shared masks"]].map(([v, t]) => el("option", { value: v }, t)));
    const lociSel = el("select", { onchange: e => setOrder("loci", e.target.value) }, ...[["input", "input order"], ["masked", "most masked first"], ["name", "name"], ["length", "length"], ["grouped", "grouped by shared masks"]].map(([v, t]) => el("option", { value: v }, t)));
    const search = el("input", { type: "search", placeholder: "find a taxon or locus  /", size: 26, oninput: e => doSearch(e.target.value) });
    const hide = el("input", { type: "checkbox", onchange: e => { App.hidePass = e.target.checked; redraw(); } });
    const tabBtns = ["matrix", "histograms", "table"].map(t => el("button", { class: t === tabs ? "primary" : "", "data-tab": t, onclick: () => setTab(t) }, t));
    const legend = el("div", { class: "legend" });
    const legendParts = [];
    legendParts.push(el("span", { class: "fam" }, el("span", { class: "lab" }, "not masked, fails"), ...d3.range(D.M).map(f => el("span", { class: "fam", title: D.labelOf(ST.GREY + f) }, sw(D.greys[f]), String(f)))));
    legendParts.push(el("span", { class: "fam" }, el("span", { class: "lab" }, "masked, too divergent, fails"), ...d3.range(D.M, D.K + 1).map(f => el("span", { class: "fam", title: D.labelOf(ST.RED + f) }, sw(D.reds[f - D.M]), String(f)))));
    legendParts.push(el("span", { class: "fam" }, el("span", { class: "lab" }, "masked, too similar, pairs"), ...[1, 2, 3].map(k => el("span", { class: "fam", title: D.labelOf(ST.BLUE + k) }, sw(D.blues[k - 1]), k === 3 ? "3+" : String(k)))));
    legendParts.push(el("span", { class: "fam", title: "present, not tested" }, sw("hatch"), "not tested"), el("span", { class: "fam", title: "no residue" }, sw(null), "absent"));
    legend.append(...legendParts, el("span", { class: "lab" }, `of ${D.K} comparators`));
    const controls = el("div", { class: "controls" },
      el("span", null, ...tabBtns),
      el("label", null, "taxa ", taxaSel), el("label", null, "loci ", lociSel), search, el("label", null, hide, " hide passing"),
      el("button", { title: "download the cell table as TSV", onclick: Table.download }, "⤓ cells.tsv"),
      el("button", { title: "export this view as PNG", onclick: exportPng }, "⤓ PNG"));
    P.body.append(controls, legend);
    legend.style.cssText = "padding: 2px 10px 6px; border-bottom: 1px solid var(--grid);";
    // stat tiles and worst lists, in a summary block that folds away when the panel is short
    const summary = el("div", { class: "summary" }); P.body.append(summary); viewEl.summary = summary;
    const sumBtn = el("button", { title: "show or hide the summary tiles and the worst-ten lists", onclick: () => { summaryUser = summary.style.display === "none"; autoSummary(); fitAll(); redraw(); } }, "summary ▾");
    controls.insertBefore(sumBtn, controls.children[1]); viewEl.sumBtn = sumBtn;
    const tiles = el("div", { class: "tiles" }); summary.append(tiles); viewEl.tiles = tiles;
    const worst = el("div", { class: "worst" }); summary.append(worst); viewEl.worst = worst;
    // the matrix
    wrap = el("div", { class: "matrix-wrap" });
    rowCanvas = el("canvas", { class: "rowlabels" }); colCanvas = el("canvas", { class: "collabels" }); matrix = el("canvas", { class: "matrix", tabindex: 0 });
    wrap.append(colCanvas, rowCanvas, matrix); P.body.append(wrap); ctx = matrix.getContext("2d");
    const stripWrap = el("div", { class: "strip-wrap" }); strip = el("canvas"); stripCtx = strip.getContext("2d");
    stripWrap.append(strip, el("span", { class: "hint" }, "wheel: zoom · shift+wheel: loci only · alt+wheel: taxa only · drag: pan · shift+drag: zoom to box · double-click: reset · click a label: open its panel · click a cell: open the sequence panel"));
    P.body.append(stripWrap);
    viewEl.hist = el("div", { class: "hist-wrap", style: "display:none; padding: 8px 10px; overflow:auto; flex: 1 1 auto;" }); P.body.append(viewEl.hist);
    viewEl.table = el("div", { class: "table-wrap", style: "display:none; flex: 1 1 auto;" }); P.body.append(viewEl.table);
    viewEl.matrixParts = [wrap, stripWrap];
    P.body.style.display = "flex"; P.body.style.flexDirection = "column";
    bindMatrix(); bindStrip();
    return { el: P.el, redraw, render };
  }
  function setTab(t) { tabs = t; for (const b of P.body.querySelectorAll("[data-tab]")) b.className = b.dataset.tab === t ? "primary" : ""; viewEl.matrixParts.forEach(e => e.style.display = t === "matrix" ? "" : "none"); viewEl.hist.style.display = t === "histograms" ? "" : "none"; viewEl.table.style.display = t === "table" ? "" : "none"; P.body.style.overflowY = t === "histograms" ? "auto" : "hidden"; if (t === "histograms") Hist.render(viewEl.hist); if (t === "table") Table.render(viewEl.table); if (t === "matrix") redraw(); }
  function render() { renderTiles(); renderWorst(); autoSummary(); fitAll(); redraw(); }
  function renderTiles() {
    const D = App.D, c = D.man.counts; const t = viewEl.tiles; t.replaceChildren();
    const tile = (v, k, title) => el("div", { class: "tile", title: title || "" }, el("div", { class: "v" }, v), el("div", { class: "k" }, k));
    const taxaWith = d3.sum(D.taxCount.masked, v => v > 0), lociWith = d3.sum(D.locCount.masked, v => v > 0);
    const near = Object.values(c.near_miss || {}).reduce((a, b) => a + b, 0);
    t.append(tile(fmtInt(c.masked), `sequences masked · ${pct(c.masked, c.cells)} of ${fmtInt(c.cells)} tested`),
      tile(fmtInt(c.masked_divergent), "too divergent" + (c.ambiguous ? ` · ${c.ambiguous} with fewer comparators` : "")),
      tile(fmtInt(c.masked_similar), "too similar" + (c.masked_both ? ` · ${c.masked_both} both sides` : "")),
      tile(fmtInt(near), `near misses · fail 1 to ${D.M - 1} comparators, not masked`),
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
    if (taxa.length) w.append(col("most masked taxa (share of tested loci)", taxa, g => Detail.open("taxon", g)));
    if (loci.length) w.append(col("most masked loci (sequences)", loci, l => Detail.open("locus", l)));
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
    // taxa by the first locus at which they are masked, loci in masked-count order (and the reverse for loci): shared causes become blocks
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
  function labelSpace() { const rl = V.ch >= 10 ? 150 : 0, cl = V.cw >= 10 ? 90 : 0; wrap.style.setProperty("--rowlab", rl + "px"); wrap.style.setProperty("--collab", cl + "px"); return { rl, cl }; }
  function fitAll() { const D = App.D; const { W, H } = size(); const { rl, cl } = labelSpace(); V.cw = Math.max(0.1, (W - rl) / D.nloc); V.ch = Math.max(0.1, (H - cl) / D.n); V.x0 = 0; V.y0 = 0; labelSpace(); }
  function canvasArea() { const { W, H } = size(); const { rl, cl } = labelSpace(); return { w: Math.max(10, W - rl), h: Math.max(10, H - cl), rl, cl }; }
  function redraw() {
    if (!App.D || tabs !== "matrix") return;
    const was = summaryShown; autoSummary(); if (was !== summaryShown) fitAll();
    const D = App.D; const { w, h, rl, cl } = canvasArea(); const dpr = devicePixelRatio || 1;
    matrix.width = Math.round(w * dpr); matrix.height = Math.round(h * dpr); matrix.style.width = w + "px"; matrix.style.height = h + "px";
    // clamp the view
    V.x0 = clamp(V.x0, 0, Math.max(0, D.nloc - w / V.cw)); V.y0 = clamp(V.y0, 0, Math.max(0, D.n - h / V.ch));
    const img = ctx.createImageData(matrix.width, matrix.height); const px = img.data; const lut = D.lut; const sev = D.sev, st = D.state;
    const ot = App.order.taxa, ol = App.order.loci; const nloc = D.nloc;
    const pw = V.cw * dpr, ph = V.ch * dpr; const surf = hexToRgb(cssVar("--surface"));
    const hide = App.hidePass;
    // per pixel column: the display columns it covers
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
        const o = rowBase + x * 4;
        if (hide && bestCode >= ST.GREY && bestCode < ST.RED) { px[o] = surf[0]; px[o + 1] = surf[1]; px[o + 2] = surf[2]; px[o + 3] = 255; }
        else { const l = bestCode * 4; px[o] = lut[l]; px[o + 1] = lut[l + 1]; px[o + 2] = lut[l + 2]; px[o + 3] = 255; }
      }
    }
    ctx.putImageData(img, 0, 0);
    // hatch for untested: draw diagonal lines over untested cells when cells are large enough
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
    if (rl) { const c = rowCanvas.getContext("2d"); c.scale(dpr, dpr); c.fillStyle = cssVar("--ink"); c.font = `${Math.min(12, V.ch - 2)}px ${cssVar("--sans")}`; c.textBaseline = "middle"; c.textAlign = "right";
      const r0 = Math.floor(V.y0), r1 = Math.min(D.n - 1, Math.ceil(V.y0 + h / V.ch));
      for (let r = r0; r <= r1; r++) { const nm = D.taxa.names[ot[r]]; c.fillStyle = D.taxa.status[ot[r]] === "tested" ? cssVar("--ink") : cssVar("--faint"); c.fillText(trunc(c, nm, rl - 8), rl - 4, (r - V.y0) * V.ch + V.ch / 2); } }
    if (cl) { const c = colCanvas.getContext("2d"); c.scale(dpr, dpr); c.fillStyle = cssVar("--ink"); c.font = `${Math.min(12, V.cw - 2)}px ${cssVar("--sans")}`; c.textBaseline = "middle"; c.textAlign = "left";
      const c0 = Math.floor(V.x0), c1 = Math.min(D.nloc - 1, Math.ceil(V.x0 + w / V.cw));
      for (let k = c0; k <= c1; k++) { c.save(); c.translate((k - V.x0) * V.cw + V.cw / 2, cl - 4); c.rotate(-Math.PI / 2); c.fillStyle = D.loci.status[ol[k]] === "tested" ? cssVar("--ink") : cssVar("--faint"); c.fillText(trunc(c, D.loci.names[ol[k]], cl - 8), 0, 0); c.restore(); } }
  }
  function trunc(c, s, maxw) { if (c.measureText(s).width <= maxw) return s; let lo = 0, hi = s.length; while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (c.measureText(s.slice(0, mid) + "…").width <= maxw) lo = mid; else hi = mid - 1; } return s.slice(0, lo) + "…"; }
  function drawStrip() {
    const D = App.D; const dpr = devicePixelRatio || 1; const sw_ = Math.min(700, Math.max(200, size().W - 40)), sh = 48;
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
    matrix.addEventListener("wheel", ev => {
      ev.preventDefault(); const D = App.D; const { x, y, col, row } = cellAt(ev); const f = Math.exp(-ev.deltaY * 0.002);
      const zx = !ev.altKey, zy = !ev.shiftKey;
      if (zx) { const ncw = clamp(V.cw * f, 0.05, 40); V.x0 = V.x0 + x / V.cw - x / ncw; V.cw = ncw; }
      if (zy) { const nch = clamp(V.ch * f, 0.05, 40); V.y0 = V.y0 + y / V.ch - y / nch; V.ch = nch; }
      labelSpace(); redraw(); writeHash();
    }, { passive: false });
    matrix.addEventListener("pointerdown", ev => { matrix.setPointerCapture(ev.pointerId); const c = cellAt(ev); dragging = { x: c.x, y: c.y, x0: V.x0, y0: V.y0, box: ev.shiftKey, moved: false }; });
    matrix.addEventListener("pointermove", ev => {
      const c = cellAt(ev); lastMouse = ev;
      if (dragging) { const dx = c.x - dragging.x, dy = c.y - dragging.y; if (Math.abs(dx) + Math.abs(dy) > 3) dragging.moved = true;
        if (dragging.box) { redraw(); const dpr = devicePixelRatio || 1; ctx.save(); ctx.scale(dpr, dpr); ctx.strokeStyle = cssVar("--focus"); ctx.setLineDash([4, 3]); ctx.strokeRect(dragging.x, dragging.y, dx, dy); ctx.restore(); }
        else if (dragging.moved) { V.x0 = dragging.x0 - dx / V.cw; V.y0 = dragging.y0 - dy / V.ch; redraw(); }
        return; }
      hoverAt(c, ev);
    });
    matrix.addEventListener("pointerup", ev => {
      const c = cellAt(ev); const d = dragging; dragging = null; if (!d) return;
      if (d.box && d.moved) { const { w, h } = canvasArea(); const x0 = Math.min(d.x, c.x), x1 = Math.max(d.x, c.x), y0 = Math.min(d.y, c.y), y1 = Math.max(d.y, c.y);
        const c0 = V.x0 + x0 / V.cw, c1 = V.x0 + x1 / V.cw, r0 = V.y0 + y0 / V.ch, r1 = V.y0 + y1 / V.ch;
        V.cw = clamp(w / Math.max(c1 - c0, 1), 0.05, 40); V.ch = clamp(h / Math.max(r1 - r0, 1), 0.05, 40); V.x0 = c0; V.y0 = r0; labelSpace(); redraw(); writeHash(); return; }
      if (!d.moved) clickCell(c);
      writeHash();
    });
    matrix.addEventListener("dblclick", () => { fitAll(); redraw(); writeHash(); });
    matrix.addEventListener("pointerleave", () => { hideHover(); hoverCells = null; });
    rowCanvas.addEventListener("click", ev => { const r = rowCanvas.getBoundingClientRect(); const row = Math.floor(V.y0 + (ev.clientY - r.top) / V.ch); if (row >= 0 && row < App.D.n) Detail.open("taxon", App.order.taxa[row]); });
    colCanvas.addEventListener("click", ev => { const r = colCanvas.getBoundingClientRect(); const col = Math.floor(V.x0 + (ev.clientX - r.left) / V.cw); if (col >= 0 && col < App.D.nloc) Detail.open("locus", App.order.loci[col]); });
    rowCanvas.style.cursor = "pointer"; colCanvas.style.cursor = "pointer";
    rowCanvas.addEventListener("mousemove", ev => { const r = rowCanvas.getBoundingClientRect(); const row = Math.floor(V.y0 + (ev.clientY - r.top) / V.ch); if (row >= 0 && row < App.D.n) showHover(unitHover("taxon", App.order.taxa[row]), ev.clientX, ev.clientY); });
    colCanvas.addEventListener("mousemove", ev => { const r = colCanvas.getBoundingClientRect(); const col = Math.floor(V.x0 + (ev.clientX - r.left) / V.cw); if (col >= 0 && col < App.D.nloc) showHover(unitHover("locus", App.order.loci[col]), ev.clientX, ev.clientY); });
    rowCanvas.addEventListener("mouseleave", hideHover); colCanvas.addEventListener("mouseleave", hideHover);
    document.addEventListener("keydown", ev => {
      if (ev.target.tagName === "INPUT" || ev.target.tagName === "SELECT" || ev.target.tagName === "TEXTAREA") { if (ev.key === "Escape") ev.target.blur(); return; }
      if (ev.key === "/") { ev.preventDefault(); $("input[type=search]", P.el).focus(); return; }
      if (ev.key === "Escape") { if (App.panels.sequence) closePanel("sequence"); else if (App.panels.detail) closePanel("detail"); return; }
      const step = 40;
      if (ev.key === "ArrowLeft") V.x0 -= step / V.cw; else if (ev.key === "ArrowRight") V.x0 += step / V.cw; else if (ev.key === "ArrowUp") V.y0 -= step / V.ch; else if (ev.key === "ArrowDown") V.y0 += step / V.ch;
      else if (ev.key === "+" || ev.key === "=") { V.cw *= 1.25; V.ch *= 1.25; } else if (ev.key === "-") { V.cw /= 1.25; V.ch /= 1.25; } else return;
      ev.preventDefault(); labelSpace(); redraw(); writeHash();
    });
  }
  function bindStrip() {
    const move = ev => { const D = App.D; const r = strip.getBoundingClientRect(); const { w, h } = canvasArea(); const fx = (ev.clientX - r.left) / r.width, fy = (ev.clientY - r.top) / r.height; V.x0 = fx * D.nloc - w / V.cw / 2; V.y0 = fy * D.n - h / V.ch / 2; redraw(); };
    let down = false; strip.addEventListener("pointerdown", ev => { down = true; strip.setPointerCapture(ev.pointerId); move(ev); }); strip.addEventListener("pointermove", ev => down && move(ev)); strip.addEventListener("pointerup", () => { down = false; writeHash(); });
  }
  function clickCell(c) {
    const D = App.D; if (c.row < 0 || c.row >= D.n || c.col < 0 || c.col >= D.nloc) return;
    if (V.cw < 4 || V.ch < 4) {                              // too small to pick: zoom in on it instead
      const { w, h } = canvasArea(); const fx = c.x / w, fy = c.y / h; V.cw = Math.max(V.cw * 4, 6); V.ch = Math.max(V.ch * 4, 6); V.x0 = c.col - fx * w / V.cw; V.y0 = c.row - fy * h / V.ch; labelSpace(); redraw(); return;
    }
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
  function centerOn(g, li) { const { w, h } = canvasArea(); if (g !== null) V.y0 = App.order.taxa.indexOf(g) - h / V.ch / 2; if (li !== null) V.x0 = App.order.loci.indexOf(li) - w / V.cw / 2; redraw(); }
  return { create, render, redraw, setOrder, V, centerOn, highlight, labelSpace, fitAll };
})();

// hover content for a cell, shared by the matrix, the scatter and the ribbons
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
    if (bm !== null && bm >= 0 && !compact) { const close = Array.from(D.closest.subarray(g * 5, g * 5 + 5)).includes(bm); row("best match", `${D.taxa.names[bm]} (${fmt1(D.num("best_match_identity", g, li))}% identity${close ? ", one of the 5 closest relatives" : ""})`); }
    const pt = D.num("partner", g, li);
    if (pt !== null && pt >= 0) { const lam = D.num("lam", g, li), lp = D.num("lam_p", g, li); row("too similar to", `${D.taxa.names[pt]}${sd & 8 ? " (identified as the recipient)" : " (direction not identified)"}`); if (lam !== null) row("direction LRT", `|2 ln LR| = ${fmt1(lam)}, p = ${lp < 1e-4 ? lp.toExponential(1) : fmt3(lp)} (χ², 1 df)`); }
    if (!compact) row("", el("span", { class: "more" }, "click: sequence panel"));
  }
  box.append(t); return box;
}
function unitHover(kind, i) {
  const D = App.D; const c = kind === "taxon" ? D.taxCount : D.locCount; const nm = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
  const t = el("table"); const row = (k, v) => t.append(el("tr", null, el("td", { class: "k" }, k), el("td", null, v)));
  row(kind, nm); row("status", kind === "taxon" ? D.taxa.status[i] : D.loci.status[i]);
  row("sequences", `${fmtInt(c.present[i])} present, ${fmtInt(c.tested[i])} tested`);
  row("masked", `${fmtInt(c.masked[i])} (${pct(c.masked[i], c.tested[i])}): ${c.div[i]} too divergent, ${c.sim[i]} too similar; ${c.near[i]} near misses`);
  if (kind === "locus") { row("length", `${fmtInt(D.loci.length[i])} sites`); if (D.loci.rate[i] !== null) row("rate", `${fmt3(D.loci.rate[i])} × the typical locus`); row("plateau", `${fmt1(D.loci.plateau[i])}%`); if (D.loci.phi[i] !== null) row("dispersion φ", fmt3(D.loci.phi[i])); }
  else { const cl = Array.from(D.closest.subarray(i * 5, i * 5 + 5)).filter(x => x >= 0).map(x => D.taxa.names[x]); if (cl.length) row("closest by m", cl.join(", ")); }
  return el("div", null, t, el("div", { class: "more" }, "click: open the panel"));
}

// ------------------------------------------------------------------ the taxon / locus panel
const Detail = (() => {
  let cur = null;
  function open(kind, i) {
    const D = App.D; cur = { kind, i };
    if (!App.panels.detail) { const P = makePanel("detail", [], () => closePanel("detail")); App.panels.detail = { el: P.el, P, render, redraw: render }; App.panelsEl.append(P.el); }
    if (App.panels.sequence) closePanel("sequence");
    render(); writeHash(); Overview.highlight(kind === "taxon" ? "row" : "col", i);
  }
  function cellsOf() {
    const D = App.D; const out = [];
    if (cur.kind === "taxon") { const g = cur.i; for (let li = 0; li < D.nloc; li++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    else { const li = cur.i; for (let g = 0; g < D.n; g++) if (D.cols.n_comp[g * D.nloc + li] > 0) out.push([g, li]); }
    return out;
  }
  function render() {
    if (!cur || !App.panels.detail) return;
    const D = App.D; const { kind, i } = cur; const P = App.panels.detail.P; const name = kind === "taxon" ? D.taxa.names[i] : D.loci.names[i];
    const cnt = kind === "taxon" ? D.taxCount : D.locCount;
    P.setCrumb([D.man.dataset || D.man.input.file, kind, el("span", null, name, el("span", { class: "sub" }, `${fmtInt(cnt.masked[i])} of ${fmtInt(cnt.tested[i])} sequences masked (${pct(cnt.masked[i], cnt.tested[i])})`))]);
    P.body.replaceChildren();
    const grid = el("div", { class: "detail" + (P_width() < 660 ? " single" : "") }); P.body.append(grid);
    const cells = cellsOf();
    grid.append(scatterCard(cells), residualCard(cells));
    if (kind === "taxon") grid.append(matchesCard(cells));
    grid.append(ribbonsCard(cells));
  }
  // observed vs expected divergence, log axes with a zero floor, the diagonal and the two effect floors
  function scatterCard(cells) {
    const D = App.D; const card = el("div", { class: "card" }, el("h3", null, "observed vs expected divergence", el("span", { class: "sub" }, "each point is one sequence; medians over its comparators")));
    const wrap = el("div", { class: "scatter-wrap" }); card.append(wrap);
    const pts = cells.map(([g, li]) => { const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li); return mi === null ? null : { g, li, obs: 100 - mi, exp: 100 - me, code: D.state[g * D.nloc + li], sd: D.cols.side[g * D.nloc + li] }; }).filter(Boolean);
    if (!pts.length) { card.append(el("div", { class: "note" }, "the bundle carries no numbers for these cells")); return card; }
    const W = cardWidth(), H = Math.min(W, 420), m = { l: 46, r: 12, t: 8, b: 40 };
    const posX = pts.filter(p => p.exp > 0).map(p => p.exp), posY = pts.filter(p => p.obs > 0).map(p => p.obs);
    const fx = posX.length ? d3.min(posX) / 2 : 0.01, fy = posY.length ? d3.min(posY) / 2 : 0.01;
    const X = p => p.exp > 0 ? p.exp : fx, Y = p => p.obs > 0 ? p.obs : fy;
    const xs = d3.scaleLog().domain([d3.min(pts, X) / 1.3, d3.max(pts, X) * 1.3]).range([m.l, W - m.r]); const ys = d3.scaleLog().domain([d3.min(pts, Y) / 1.3, d3.max(pts, Y) * 1.3]).range([H - m.b, m.t]);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    const ticks = s => s.ticks(5).filter(v => /^[1-9](\.5)?$/.test(d3.format(".1~g")(v).replace(/^0\.0*|^\d+/, mm => mm.replace(/\d/g, "").length ? mm : mm[0])) || true);
    const xAxis = d3.axisBottom(xs).ticks(5, "~g"), yAxis = d3.axisLeft(ys).ticks(5, "~g");
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(xAxis);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(yAxis);
    if (posX.length < pts.length) svg.append("text").attr("x", xs(fx)).attr("y", H - m.b + 14).attr("text-anchor", "middle").text("0");
    if (posY.length < pts.length) svg.append("text").attr("x", m.l - 6).attr("y", ys(fy) + 4).attr("text-anchor", "end").text("0");
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("expected divergence of the sequence (%)");
    svg.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("observed divergence (%)");
    const line = (f, dash, col) => { const [x0, x1] = xs.domain(); const p0 = [Math.max(x0, ys.domain()[0] / f), Math.min(x1, ys.domain()[1] / f)]; if (p0[0] < p0[1]) svg.append("line").attr("x1", xs(p0[0])).attr("y1", ys(p0[0] * f)).attr("x2", xs(p0[1])).attr("y2", ys(p0[1] * f)).attr("stroke", col).attr("stroke-width", f === 1 ? 1.2 : 0.9).attr("stroke-dasharray", dash); };
    line(1, null, cssVar("--model")); line(D.man.options.min_ratio, "4 3", cssVar("--floor")); line(1 / 3, "4 3", cssVar("--floor"));
    wrap.append(svg.node());
    const cv = el("canvas"); const dpr = devicePixelRatio || 1; cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px"; wrap.append(cv);
    const c = cv.getContext("2d"); c.scale(dpr, dpr);
    const drawPt = (p, r, hl) => { const x = xs(X(p)), y = ys(Y(p)); const col = D.colorOf(p.code) || cssVar("--grey-1"); c.fillStyle = col; c.strokeStyle = hl ? cssVar("--focus") : cssVar("--surface"); c.lineWidth = hl ? 2 : 0.8; c.beginPath();
      if (p.code >= ST.BLUE) { c.moveTo(x, y + r); c.lineTo(x - r, y - r * 0.7); c.lineTo(x + r, y - r * 0.7); c.closePath(); } else if (p.code >= ST.RED) { c.moveTo(x, y - r); c.lineTo(x - r, y + r * 0.7); c.lineTo(x + r, y + r * 0.7); c.closePath(); } else c.arc(x, y, r * 0.75, 0, 2 * Math.PI);
      c.fill(); c.stroke(); };
    const order = pts.slice().sort((a, b) => D.sev[a.g * D.nloc + a.li] - D.sev[b.g * D.nloc + b.li]);
    const r = pts.length > 2000 ? 2.4 : pts.length > 500 ? 3.2 : 4.2;
    for (const p of order) drawPt(p, p.code >= ST.RED ? r + 1.5 : r, false);
    // hover and click: nearest point within 9 px
    svg.style("pointer-events", "all");
    const nearest = ev => { const [mx, my] = d3.pointer(ev, svg.node()); let best = null, bd = 81; for (const p of pts) { const dx = xs(X(p)) - mx, dy = ys(Y(p)) - my; const d = dx * dx + dy * dy; if (d < bd) { bd = d; best = p; } } return best; };
    svg.on("mousemove", ev => { const p = nearest(ev); if (p) { showHover(cellHover(p.g, p.li), ev.clientX, ev.clientY); App.sel.cells.clear(); App.sel.cells.add(p.g * D.nloc + p.li); Overview.redraw(); } else hideHover(); });
    svg.on("mouseleave", () => { hideHover(); App.sel.cells.clear(); Overview.redraw(); });
    svg.on("click", ev => { const p = nearest(ev); if (p) Sequence.open(p.g, p.li); });
    card.append(el("div", { class: "note" }, `diagonal: observed = expected; dashed: the effect floors (${D.man.options.min_ratio}× above, ⅓ below). Shapes: circles pass, ▲ masked too divergent, ▼ masked too similar. Hover a point for its numbers; click it for the sequence panel.`));
    return card;
  }
  function P_width() { return App.panels.detail ? App.panels.detail.P.body.getBoundingClientRect().width || 800 : 800; }
  function cardWidth() { const w = P_width() - 20; return w >= 660 ? Math.max(280, Math.min(620, (w - 8) / 2 - 18)) : Math.max(260, w - 18); }   // two cards a row when there is room, else one
  // the residual strip: excess_m on a symlog axis, masked as triangles
  function residualCard(cells) {
    const D = App.D; const card = el("div", { class: "card" }, el("h3", null, "excess divergence beyond expectation", el("span", { class: "sub" }, "in units of a typical locus's divergence (m)")));
    const vals = cells.map(([g, li]) => ({ g, li, x: D.num("excess_m", g, li), code: D.state[g * D.nloc + li] })).filter(p => p.x !== null && Number.isFinite(p.x));
    if (!vals.length) { card.append(el("div", { class: "note" }, "no residuals carried for these cells")); return card; }
    const W = cardWidth(), H = 150, m = { l: 14, r: 14, t: 10, b: 36 };
    const lo = Math.min(d3.min(vals, p => p.x) * 1.3, -30), hi = Math.max(d3.max(vals, p => p.x) * 1.3, 300);
    const xs = d3.scaleSymlog().constant(1).domain([lo, hi]).range([m.l, W - m.r]);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    const tv = [-1000, -100, -10, -1, 0, 1, 10, 100, 1000].filter(v => v >= lo && v <= hi);
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).tickValues(tv).tickFormat(d3.format(",")));
    svg.append("line").attr("x1", xs(0)).attr("x2", xs(0)).attr("y1", m.t).attr("y2", H - m.b).attr("stroke", cssVar("--model"));
    if (cur.kind === "locus") {                                       // where a saturated sequence would land: from the plateau and rate
      const li = cur.i, P = D.loci.plateau[li], rate = D.loci.rate[li];
      if (rate) { const corr = p => -P * Math.log(1 - Math.min(Math.max(p, 0), 0.98 * P) / P); const ce = vals.map(p => (corr(0.98 * P) - corr(100 - (D.num("med_expected", p.g, p.li) ?? 100))) / rate).filter(Number.isFinite).sort(d3.ascending); if (ce.length) { const c10 = d3.quantileSorted(ce, 0.1); svg.append("rect").attr("x", xs(c10)).attr("y", m.t).attr("width", Math.max(0, W - m.r - xs(c10))).attr("height", H - m.b - m.t).attr("fill", cssVar("--surface-3")).attr("opacity", 0.8); svg.append("text").attr("x", W - m.r - 4).attr("y", m.t + 12).attr("text-anchor", "end").text("saturated sequences land here"); } }
    }
    svg.append("text").attr("x", W / 2).attr("y", H - 4).attr("text-anchor", "middle").text("◀ too similar     excess divergence (m)     too divergent ▶");
    const rng = d3.randomUniform.source(d3.randomLcg(42))(-18, 18); const cy = (m.t + H - m.b) / 2;
    const g = svg.append("g");
    for (const p of vals.slice().sort((a, b) => D.sev[a.g * D.nloc + a.li] - D.sev[b.g * D.nloc + b.li])) {
      const col = D.colorOf(p.code) || cssVar("--grey-1"); const x = xs(p.x), y = cy + rng(); const masked = p.code >= ST.RED;
      const node = masked ? g.append("path").attr("d", p.code >= ST.BLUE ? d3.symbol(d3.symbolTriangle, 50)() : d3.symbol(d3.symbolTriangle, 50)()).attr("transform", `translate(${x},${y}) rotate(${p.code >= ST.BLUE ? 180 : 0})`).attr("fill", col).attr("stroke", cssVar("--surface")) : g.append("circle").attr("cx", x).attr("cy", y).attr("r", vals.length > 500 ? 1.6 : 2.4).attr("fill", col).attr("opacity", vals.length > 500 ? 0.5 : 0.8);
      node.style("cursor", "pointer").on("mousemove", ev => showHover(cellHover(p.g, p.li), ev.clientX, ev.clientY)).on("mouseleave", hideHover).on("click", () => Sequence.open(p.g, p.li));
    }
    card.append(svg.node());
    card.append(el("div", { class: "note" }, "cells identical to their comparators sit at the left; a tight cluster is one cause, a scatter is many independent problems"));
    return card;
  }
  // best matches of the masked sequences, with the count expected from the unmasked ones
  function matchesCard(cells) {
    const D = App.D; const g0 = cur.i; const card = el("div", { class: "card wide" }, el("h3", null, "best match of each masked sequence", el("span", { class: "sub" }, "the taxon each masked sequence matches best at its locus")));
    const close = new Set(Array.from(D.closest.subarray(g0 * 5, g0 * 5 + 5)).filter(x => x >= 0));
    const masked = [], unmasked = [];
    for (const [g, li] of cells) { const bm = D.num("best_match", g, li); if (bm === null || bm < 0) continue; ((D.cols.side[g * D.nloc + li] & 3) ? masked : unmasked).push({ g, li, bm, id: D.num("best_match_identity", g, li), code: D.state[g * D.nloc + li] }); }
    if (!masked.length) { card.append(el("div", { class: "note" }, unmasked.length ? "no masked sequence for this taxon" : "best matches were not carried for this taxon (the dataset is large: they are carried for cells with a fail only)")); return card; }
    const cm = d3.rollup(masked, v => v.length, d => d.bm), cu = d3.rollup(unmasked, v => v.length, d => d.bm); const nm = masked.length, nu = unmasked.length;
    let taxaSorted = Array.from(new Set([...cm.keys(), ...Array.from(cu.keys()).filter(u => (cu.get(u) / Math.max(nu, 1)) * nm >= 1)])).sort((a, b) => (cm.get(b) || 0) - (cm.get(a) || 0) || (cu.get(b) || 0) - (cu.get(a) || 0));
    const MAXB = 24; let other = null; if (taxaSorted.length > MAXB) { const rest = taxaSorted.slice(MAXB - 1); taxaSorted = taxaSorted.slice(0, MAXB - 1); other = { n: rest.length, cm: d3.sum(rest, u => cm.get(u) || 0), cu: d3.sum(rest, u => cu.get(u) || 0), set: new Set(rest) }; }
    const W = Math.max(400, P_width() - 40), H = 230, m = { l: 40, r: 10, t: 10, b: 110 };
    const cats = taxaSorted.map(u => ({ key: u, name: D.taxa.names[u], cm: cm.get(u) || 0, cu: cu.get(u) || 0, rows: masked.filter(d => d.bm === u) })); if (other) cats.push({ key: -1, name: `other (${other.n} taxa)`, cm: other.cm, cu: other.cu, rows: masked.filter(d => other.set.has(d.bm)) });
    const xs = d3.scaleBand().domain(cats.map(c => c.name)).range([m.l, W - m.r]).padding(0.18); const top = Math.max(d3.max(cats, c => Math.max(c.cm, c.cu / Math.max(nu, 1) * nm)), 1);
    const ys = d3.scaleLinear().domain([0, top * 1.15]).range([H - m.b, m.t]);
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).ticks(4).tickFormat(d3.format("d")));
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
    card.append(svg.node());
    const pm = 100 * d3.mean(masked, d => close.has(d.bm) ? 1 : 0), pu = nu ? 100 * d3.mean(unmasked, d => close.has(d.bm) ? 1 : 0) : NaN;
    card.append(el("div", { class: "note" }, `bars: masked sequences by best match, stacked by mask kind; the tick is the count expected if they matched like the taxon's ${fmtInt(nu)} unmasked sequences; bold names are the 5 closest relatives genome-wide. Best match is one of the 5 closest for ${fmt1(pm)}% of masked and ${Number.isFinite(pu) ? fmt1(pu) : "–"}% of unmasked sequences; median identity to it ${fmt1(d3.median(masked, d => d.id))}% masked vs ${nu ? fmt1(d3.median(unmasked, d => d.id)) : "–"}% unmasked. Bars on their ticks are damaged copies of the right sequence (a wrong exon); bars shifted onto a distant taxon point to contamination or a shared paralogue.`));
    return card;
  }
  // site ribbons: one pixel-row per sequence, coloured by the site code; masked first, then the near-misses
  function ribbonsCard(cells) {
    const D = App.D; const card = el("div", { class: "card wide" }, el("h3", null, "site ribbons", el("span", { class: "sub" }, "each row is one sequence, sites along x; red where a majority of its comparators differ, blue where a pair shares a state against the others")));
    const wrap = el("div", { class: "ribbons-wrap" }); card.append(wrap);
    const withCodes = cells.filter(([g, li]) => { const i = g * D.nloc + li; return (D.cols.side[i] & 3) || D.cols.n_fail[i] > 0 || D.cols.n_similar[i] > 0; });
    const masked = withCodes.filter(([g, li]) => D.cols.side[g * D.nloc + li] & 3), near = withCodes.filter(([g, li]) => !(D.cols.side[g * D.nloc + li] & 3));
    if (!withCodes.length) { card.append(el("div", { class: "note" }, "no masked sequence or near miss here")); return card; }
    const maxLen = cur.kind === "locus" ? D.loci.length[cur.i] : d3.max(withCodes, ([g, li]) => D.loci.length[li]);
    const addRows = (list, title) => {
      if (!list.length) return; if (title) wrap.append(el("div", { class: "rule" }, title));
      list.sort((a, b) => D.sev[b[0] * D.nloc + b[1]] - D.sev[a[0] * D.nloc + a[1]]);
      for (const [g, li] of list.slice(0, 400)) {
        const cv = el("canvas"); const row = el("div", { class: "rib-row", onclick: () => Sequence.open(g, li), onmouseenter: () => { App.sel.cells.clear(); App.sel.cells.add(g * D.nloc + li); Overview.redraw(); }, onmouseleave: () => { App.sel.cells.clear(); Overview.redraw(); } },
          el("span", { class: "nm", title: cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g] }, sw(D.colorOf(D.state[g * D.nloc + li])), cur.kind === "taxon" ? D.loci.names[li] : D.taxa.names[g]), cv);
        wrap.append(row);
        D.codes(li).then(C => drawRibbon(cv, C, g, li, maxLen, row));
      }
      if (list.length > 400) wrap.append(el("div", { class: "rule" }, `${fmtInt(list.length - 400)} more rows not drawn; use the table view`));
    };
    addRows(masked, null); addRows(near, `near misses: not masked, ${D.M - 1 > 1 ? `1 to ${D.M - 1}` : "1"} of ${D.K} comparators fail`);
    card.append(el("div", { class: "note" }, cur.kind === "locus" ? "rows share the alignment's columns: a vertical band is a column problem (misalignment), a horizontal block in one row is a taxon problem (a wrong exon, a frameshift). Hover a row for the site; click it for the sequence panel." : "each row is a different locus drawn to its own length on a shared axis. Click a row for the sequence panel."));
    return card;
  }
  function drawRibbon(cv, C, g, li, maxLen, row) {
    const D = App.D; if (!C) { cv.title = "codes not carried for this cell"; return; }
    let r = -1; for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) { r = k; break; }
    if (r < 0) { cv.title = "codes not carried for this cell"; return; }
    const len = C.len; const W = Math.max(50, Math.floor(cv.getBoundingClientRect().width || 600)); cv.width = W; cv.height = 1;
    const img = cv.getContext("2d").createImageData(W, 1); const px = img.data; const base = r * len; const surf = hexToRgb(cssVar("--surface"));
    const span = len / maxLen;                                   // the row's share of the shared axis
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
  return { open, render, cur: () => cur };
})();
// site-code helpers shared by the ribbons and the sequence panel
function codeSeverity(D, c) { if (c === 255) return 0; if (c === 254) return 5; if (c >= 128) return 30 + (c - 128); if (c >= D.M) return 20 + c; return 2 + c; }
function codeRgb(D, c) {
  if (c === 255) return hexToRgb(cssVar("--surface")); if (c === 254) return hexToRgb(cssVar("--muted"));
  if (c >= 128) { const j = c - 128; return hexToRgb(D.blues[Math.min(j - 2, 2) < 0 ? 0 : Math.min(j - 2, 2)]); }
  if (c >= D.M) return hexToRgb(D.reds[Math.min(c - D.M, D.reds.length - 1)]);
  return hexToRgb(D.greys[Math.min(c, D.greys.length - 1)]);
}
function codeLabel(D, c) { if (c === 255) return "no residue in the focal sequence"; if (c === 254) return "the pair differs here"; if (c >= 128) return `the pair shares this state while ${c - 128} other comparators differ`; return `${c} of the comparators present differ${c >= D.M ? " (a majority)" : ""}`; }
function alignmentPos(D, li, s) { let left = s; for (const [a, b, step] of D.loci.ranges[li]) { const n = Math.floor((b - a) / step) + 1; if (left < n) return a + left * step; left -= n; } return s; }

// ------------------------------------------------------------------ the sequence panel
const Sequence = (() => {
  let cur = null, pps = null;                                   // pixels per site
  function open(g, li) {
    cur = { g, li }; pps = null;
    if (!App.panels.sequence) { const P = makePanel("sequence", [], () => closePanel("sequence")); App.panels.sequence = { el: P.el, P, render, redraw: render }; App.panelsEl.append(P.el); }
    App.sel.cells.clear(); App.sel.cells.add(g * App.D.nloc + li); Overview.redraw(); render(); writeHash();
  }
  async function render() {
    if (!cur || !App.panels.sequence) return;
    const D = App.D; const { g, li } = cur; const P = App.panels.sequence.P; const i = g * D.nloc + li;
    P.setCrumb([D.man.dataset || D.man.input.file, el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("taxon", g); } }, D.taxa.names[g]), el("a", { href: "#", onclick: ev => { ev.preventDefault(); Detail.open("locus", li); } }, D.loci.names[li])]);
    P.body.replaceChildren(el("div", { class: "loading" }, "loading the sequences…"));
    const [C, S] = await Promise.all([D.codes(li), D.seq(li)]);
    const comps = D.comps(g, li); const partner = D.num("partner", g, li); const sd = D.cols.side[i];
    const rows = [{ g, role: "focal" }]; if (partner !== null && partner >= 0) rows.push({ g: partner, role: "partner" }); comps.forEach((c, r) => { if (c !== partner) rows.push({ g: c, role: `comparator ${r + 1}` }); });
    const seqRow = gg => { if (!S) return null; for (let k = 0; k < S.taxa.length; k++) if (S.taxa[k] === gg) return S.blk.subarray(k * S.len, (k + 1) * S.len); return null; };
    let codeRow = null; if (C) for (let k = 0; k < C.cells.length; k++) if (C.cells[k] === g) codeRow = C.blk.subarray(k * C.len, (k + 1) * C.len);
    const len = D.loci.length[li];
    // header line with the cell's numbers in words
    const head = el("div", { class: "seq-head" });
    const cl = D.num("closest", g, li), co = D.num("closest_overlap", g, li), ci = D.num("closest_identity", g, li), ce = D.num("closest_expected", g, li);
    const parts = [el("span", { class: "state" }, sw(D.colorOf(D.state[i])), D.labelOf(D.state[i]))];
    if (cl !== null && cl >= 0) parts.push(el("span", null, ` · against its closest comparator ${D.taxa.names[cl]}: `, el("span", { class: "num" }, `${fmtInt(Math.round((100 - ci) / 100 * co))} of ${fmtInt(co)} shared sites differ, ${fmtInt(Math.round((100 - ce) / 100 * co))} expected`)));
    if (partner !== null && partner >= 0) { const lam = D.num("lam", g, li), lp = D.num("lam_p", g, li); parts.push(el("span", null, ` · too similar to ${D.taxa.names[partner]}${sd & 8 ? ", identified as the recipient" : ", direction not identified"}${lam !== null ? ` (direction LRT |2 ln LR| = ${fmt1(lam)}, p = ${lp < 1e-4 ? lp.toExponential(1) : fmt3(lp)})` : ""}`)); }
    const xm = D.num("excess_m", g, li); if (xm !== null) parts.push(el("span", null, ` · excess ${fmt2(xm)} m`));
    head.append(...parts);
    const body = el("div", { class: "seq-body" }); const names = el("div", { class: "seq-names" }); const cwrap = el("div", { class: "seq-canvas-wrap" }); const cv = el("canvas"); cwrap.append(cv); body.append(names, cwrap);
    for (const r of rows) names.append(el("div", { class: r.role === "focal" ? "focal" : "", title: `${D.taxa.names[r.g]} (${r.role})` }, `${D.taxa.names[r.g]} `, el("span", { style: "color: var(--faint)" }, r.role === "focal" ? "" : r.role === "partner" ? "partner" : r.role.replace("comparator ", "#"))));
    const foot = el("div", { class: "seq-foot" }, el("span", null, `${fmtInt(len)} sites · wheel: zoom · drag: scroll · letters appear from 7 px per site`), el("span", { style: "flex:1" }),
      el("button", { onclick: () => { pps = cwrap.getBoundingClientRect().width / len; draw(); } }, "whole locus"), el("button", { onclick: () => { pps = 12; draw(); } }, "letters"));
    P.body.replaceChildren(head, body, foot);
    if (!S && !codeRow) { cwrap.append(el("div", { class: "note", style: "padding: 10px" }, "neither sequences nor site codes were carried for this cell in this report")); return; }
    let x0 = 0;                                                   // first visible site
    const dpr = devicePixelRatio || 1; const RH = 16, TOP = 18;
    const seqRows = rows.map(r => ({ ...r, bytes: seqRow(r.g) }));
    const draw = () => {
      const W = Math.max(100, cwrap.getBoundingClientRect().width), H = TOP + RH * rows.length + 4;
      if (pps === null) pps = W / len; pps = clamp(pps, W / len, 16); x0 = clamp(x0, 0, Math.max(0, len - W / pps));
      cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px";
      const c = cv.getContext("2d"); c.scale(dpr, dpr); c.fillStyle = cssVar("--surface"); c.fillRect(0, 0, W, H);
      const s0 = Math.floor(x0), s1 = Math.min(len - 1, Math.ceil(x0 + W / pps));
      // ruler in alignment coordinates
      c.fillStyle = cssVar("--muted"); c.font = `10px ${cssVar("--sans")}`; c.textBaseline = "top"; c.textAlign = "center";
      const stepSites = Math.max(1, Math.pow(10, Math.ceil(Math.log10(60 / pps))));
      for (let s = Math.ceil(s0 / stepSites) * stepSites; s <= s1; s += stepSites) { const x = (s - x0) * pps + pps / 2; c.fillRect(x, TOP - 4, 1, 4); c.fillText(fmtInt(alignmentPos(D, li, s) + 1), x, 1); }
      const letters = pps >= 7; c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; c.textBaseline = "middle";
      const focal = seqRows[0].bytes;
      rows.forEach((r, ri) => {
        const y = TOP + ri * RH; const bytes = seqRows[ri].bytes;
        if (ri === 0) {                                           // the focal row wears the site codes
          if (codeRow) { if (letters) { for (let s = s0; s <= s1; s++) { const rgb = codeRgb(D, codeRow[s]); c.fillStyle = `rgb(${rgb})`; c.fillRect((s - x0) * pps, y + 1, pps, RH - 2); } }
            else { const img = c.createImageData(Math.ceil(W * dpr), 1); const px = img.data; for (let x = 0; x < img.width; x++) { const a = Math.floor(x0 + x / dpr / pps), b = Math.max(a, Math.ceil(x0 + (x + 1) / dpr / pps) - 1); let best = -1, code = 255; for (let s = a; s <= Math.min(b, len - 1); s++) { const sv = codeSeverity(D, codeRow[s]); if (sv > best) { best = sv; code = codeRow[s]; } } const rgb = codeRgb(D, code); px[x * 4] = rgb[0]; px[x * 4 + 1] = rgb[1]; px[x * 4 + 2] = rgb[2]; px[x * 4 + 3] = 255; } for (let yy = 1; yy < RH - 1; yy++) c.putImageData(img, 0, Math.round((y + yy) * dpr)); } }
          else if (bytes) { c.fillStyle = cssVar("--grey-1"); for (let s = s0; s <= s1; s++) if (bytes[s] !== 45 && bytes[s] !== 63) c.fillRect((s - x0) * pps, y + 1, pps, RH - 2); }
          if (letters && bytes) { c.fillStyle = cssVar("--ink"); for (let s = s0; s <= s1; s++) c.fillText(String.fromCharCode(bytes[s]), (s - x0) * pps + pps / 2, y + RH / 2); }
        } else {
          if (!bytes) { c.fillStyle = cssVar("--faint"); c.font = `11px ${cssVar("--sans")}`; c.textAlign = "left"; c.fillText("row not carried in this report", 4, y + RH / 2); c.textAlign = "center"; c.font = `${Math.min(13, pps)}px ${cssVar("--mono")}`; return; }
          if (letters) { for (let s = s0; s <= s1; s++) { const b = bytes[s]; if (b === 45 || b === 63) continue; const diff = focal && focal[s] !== b && focal[s] !== 45 && focal[s] !== 63; c.fillStyle = diff ? cssVar("--ink") : cssVar("--faint"); c.fillText(String.fromCharCode(b), (s - x0) * pps + pps / 2, y + RH / 2); } }
          else { c.fillStyle = cssVar("--grey-1"); for (let s = s0; s <= s1; s++) if (bytes[s] !== 45 && bytes[s] !== 63) c.fillRect((s - x0) * pps, y + 4, Math.max(pps, 1), RH - 8);
            c.fillStyle = cssVar("--ink"); if (focal) for (let s = s0; s <= s1; s++) { const b = bytes[s]; if (b !== 45 && b !== 63 && focal[s] !== b && focal[s] !== 45 && focal[s] !== 63) c.fillRect((s - x0) * pps, y + 2, Math.max(pps, 1), RH - 4); } }
        }
      });
    };
    draw();
    cv.addEventListener("wheel", ev => { ev.preventDefault(); const r = cv.getBoundingClientRect(); const mx = ev.clientX - r.left; const f = Math.exp(-ev.deltaY * 0.002); const np = clamp(pps * f, r.width / len, 16); x0 = x0 + mx / pps - mx / np; pps = np; draw(); }, { passive: false });
    let drag = null; cv.addEventListener("pointerdown", ev => { cv.setPointerCapture(ev.pointerId); drag = { x: ev.clientX, x0 }; }); cv.addEventListener("pointerup", () => drag = null);
    cv.addEventListener("pointermove", ev => {
      if (drag) { x0 = drag.x0 - (ev.clientX - drag.x) / pps; draw(); return; }
      const r = cv.getBoundingClientRect(); const s = Math.floor(x0 + (ev.clientX - r.left) / pps); if (s < 0 || s >= len) { hideHover(); return; }
      const t = el("table"); t.append(el("tr", null, el("td", { class: "k" }, "site"), el("td", null, `${s + 1} of ${fmtInt(len)} · alignment position ${fmtInt(alignmentPos(D, li, s) + 1)}`)));
      if (codeRow) t.append(el("tr", null, el("td", { class: "k" }, "code"), el("td", null, codeLabel(D, codeRow[s]))));
      for (const rr of seqRows) if (rr.bytes) t.append(el("tr", null, el("td", { class: "k" }, D.taxa.names[rr.g]), el("td", { style: "font-family: var(--mono)" }, String.fromCharCode(rr.bytes[s]))));
      showHover(t, ev.clientX, ev.clientY);
    });
    cv.addEventListener("pointerleave", hideHover);
    App.panels.sequence.redraw = draw;
  }
  return { open, render, cur: () => cur };
})();

// ------------------------------------------------------------------ histograms: loci and taxa by masked count
const Hist = (() => {
  function render(root) {
    root.replaceChildren(); const D = App.D;
    root.append(facet("locus"), facet("taxon"));
  }
  function bins(kof) {
    const pos = kof.filter(k => k > 0); const kmax = d3.max(kof) || 0; let w = 1;
    if (pos.length >= 2) { const q = pos.slice().sort(d3.ascending); const iqr = d3.quantileSorted(q, 0.75) - d3.quantileSorted(q, 0.25); w = 2 * iqr / Math.cbrt(pos.length); if (w < 1) w = kmax / Math.max(1, 1 + Math.log2(pos.length)); w = Math.max(1, Math.ceil(w)); }
    const out = [[0, 0]]; for (let lo = 1; lo <= kmax; lo += w) out.push([lo, Math.min(lo + w - 1, kmax)]); return out;
  }
  function facet(by) {
    const D = App.D; const N = by === "locus" ? D.nloc : D.n; const cnt = by === "locus" ? D.locCount : D.taxCount; const names = by === "locus" ? D.loci.names : D.taxa.names;
    const units = d3.range(N).filter(u => cnt.present[u] > 0); const kof = units.map(u => cnt.masked[u]);
    const B = bins(kof); const heights = B.map(([lo, hi]) => units.filter((u, j) => kof[j] >= lo && kof[j] <= hi));
    const W = Math.max(420, Math.min(900, Overview_width() - 40)), H = 220, m = { l: 48, r: 10, t: 26, b: 42 };
    const xs = d3.scaleBand().domain(B.map((b, i) => i)).range([m.l, W - m.r]).padding(0.15);
    const ymax = d3.max(heights, h => h.length) || 1; const ys = d3.scaleSymlog().constant(2).domain([0, ymax * 1.15]).range([H - m.b, m.t]);
    const card = el("div", { class: "card", style: "margin-bottom: 10px" }, el("h3", null, `${by === "locus" ? "loci" : "taxa"} by the number of masked sequences they hold`, el("span", { class: "sub" }, `${fmtInt(units.filter((u, j) => kof[j] > 0).length)} of ${fmtInt(units.length)} have a masked sequence`)));
    const svg = d3.create("svg").attr("width", W).attr("height", H);
    svg.append("g").attr("class", "axis").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(ys).tickValues([0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].filter(v => v <= ymax * 1.15)).tickFormat(d3.format(",")));
    svg.append("g").attr("class", "axis").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(xs).tickFormat(i => { const [lo, hi] = B[i]; return lo === hi ? String(lo) : `${lo}–${hi}`; })).selectAll("text").attr("transform", B.length > 14 ? "rotate(-60)" : null).attr("text-anchor", B.length > 14 ? "end" : "middle");
    svg.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 4).attr("text-anchor", "middle").text(`sequences masked in the ${by}`);
    svg.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text(by === "locus" ? "loci" : "taxa");
    B.forEach(([lo, hi], i) => {
      const us = heights[i]; if (!us.length) return;
      // composition: the zero bar by comparators failed among its sequences; the others by mask kind
      const comp = new Map();
      if (lo === 0) { for (const u of us) for (let k = 0; k < (by === "locus" ? D.n : D.nloc); k++) { const idx = by === "locus" ? k * D.nloc + u : u * D.nloc + k; if (D.cols.n_comp[idx] === 0) continue; const code = D.state[idx]; if (code < ST.RED) comp.set(code, (comp.get(code) || 0) + 1); } }
      else { for (const u of us) for (let k = 0; k < (by === "locus" ? D.n : D.nloc); k++) { const idx = by === "locus" ? k * D.nloc + u : u * D.nloc + k; const code = D.state[idx]; if (code >= ST.RED) comp.set(code, (comp.get(code) || 0) + 1); } }
      const tot = d3.sum(comp.values()) || 1; let y0 = ys(0); const codes = Array.from(comp.keys()).sort(d3.ascending);
      for (const code of codes) { const h = (ys(0) - ys(us.length)) * comp.get(code) / tot; svg.append("rect").attr("x", xs(i)).attr("y", y0 - h).attr("width", xs.bandwidth()).attr("height", h).attr("fill", D.colorOf(code)).attr("stroke", cssVar("--surface")).attr("stroke-width", 0.8).style("cursor", "pointer")
        .on("mousemove", ev => showHover(el("div", null, el("b", null, `${fmtInt(us.length)} ${by === "locus" ? "loci" : "taxa"} with ${lo === hi ? lo : `${lo} to ${hi}`} masked sequence${hi !== 1 ? "s" : ""}`), el("div", null, `${D.labelOf(code)}: ${fmtInt(comp.get(code))} sequences (${pct(comp.get(code), tot)})`), el("div", { class: "more" }, us.slice(0, 8).map(u => names[u]).join(", ") + (us.length > 8 ? ", …" : "") + " · click to highlight them")), ev.clientX, ev.clientY))
        .on("mouseleave", hideHover).on("click", () => { if (by === "locus") { App.sel.cols = new Set(us); App.sel.rows.clear(); } else { App.sel.rows = new Set(us); App.sel.cols.clear(); } $("[data-tab=matrix]").click(); }); y0 -= h; }
      svg.append("text").attr("x", xs(i) + xs.bandwidth() / 2).attr("y", ys(us.length) - 3).attr("text-anchor", "middle").style("fill", cssVar("--ink")).text(fmtInt(us.length));
    });
    card.append(svg.node(), el("div", { class: "note" }, "zero is its own bar, coloured by how many comparators its sequences fail; the other bars are coloured by the kinds of mask they hold. The y axis is linear near zero and logarithmic above. Click a bar to highlight its units in the matrix."));
    return card;
  }
  function Overview_width() { return App.panels.overview ? App.panels.overview.el.getBoundingClientRect().width || 900 : 900; }
  return { render };
})();

// ------------------------------------------------------------------ the table view (virtualised)
const Table = (() => {
  let rows = [], sortKey = "excess_m", sortDir = -1, maskedOnly = true, root = null;
  const COLS = [["taxon", "taxon"], ["locus", "locus"], ["n_fail", "fails"], ["n_similar", "similar"], ["obs", "observed %"], ["exp", "expected %"], ["excess_m", "excess (m)"], ["bm_id", "best match %"], ["lam", "|2 ln LR|"], ["bm", "best match / partner"]];
  function build() {
    const D = App.D; rows = [];
    for (let g = 0; g < D.n; g++) for (let li = 0; li < D.nloc; li++) { const i = g * D.nloc + li; if (D.cols.n_comp[i] === 0) continue; if (maskedOnly && !(D.cols.side[i] & 3) && D.cols.n_fail[i] === 0 && D.cols.n_similar[i] === 0) continue;
      const mi = D.num("med_identity", g, li), me = D.num("med_expected", g, li); const bm = D.num("best_match", g, li), pt = D.num("partner", g, li);
      rows.push({ g, li, taxon: D.taxa.names[g], locus: D.loci.names[li], n_fail: D.cols.n_fail[i], n_similar: D.cols.n_similar[i], obs: mi === null ? null : 100 - mi, exp: me === null ? null : 100 - me, excess_m: D.num("excess_m", g, li), bm_id: D.num("best_match_identity", g, li), lam: D.num("lam", g, li), bm: pt !== null && pt >= 0 ? "partner " + D.taxa.names[pt] : bm !== null && bm >= 0 ? D.taxa.names[bm] : "", code: D.state[i] }); }
    sort();
  }
  function sort() { rows.sort((a, b) => { const x = a[sortKey], y = b[sortKey]; if (x === null || x === undefined) return 1; if (y === null || y === undefined) return -1; return (typeof x === "string" ? x.localeCompare(y) : x - y) * sortDir; }); }
  function render(r) {
    root = r; root.replaceChildren(); build(); const D = App.D;
    const head = el("div", { class: "thead" }, ...COLS.map(([k, t]) => el("span", { onclick: () => { if (sortKey === k) sortDir = -sortDir; else { sortKey = k; sortDir = typeof rows[0]?.[k] === "string" ? 1 : -1; } sort(); paint(); } }, t + (sortKey === k ? (sortDir < 0 ? " ▾" : " ▴") : ""))));
    const body = el("div", { class: "tbody" }); const inner = el("div"); body.append(inner);
    const foot = el("div", { class: "tfoot" }, el("label", null, el("input", { type: "checkbox", checked: maskedOnly ? "" : null, onchange: e => { maskedOnly = e.target.checked; build(); paint(); } }), " masked and near-miss cells only"), el("span", { class: "n" }));
    root.append(head, body, foot);
    const RH = 22;
    const paint = () => { inner.style.height = rows.length * RH + "px"; inner.style.position = "relative"; inner.replaceChildren(); const top = body.scrollTop; const i0 = Math.max(0, Math.floor(top / RH) - 5), i1 = Math.min(rows.length, i0 + Math.ceil(body.clientHeight / RH) + 10);
      for (let i = i0; i < i1; i++) { const r = rows[i]; inner.append(el("div", { class: "trow", style: `position:absolute; top:${i * RH}px; left:0; right:0`, onclick: () => Sequence.open(r.g, r.li), onmouseenter: ev => { App.sel.cells.clear(); App.sel.cells.add(r.g * D.nloc + r.li); } },
        el("span", null, sw(D.colorOf(r.code)), r.taxon), el("span", null, r.locus), el("span", null, r.n_fail), el("span", null, r.n_similar), el("span", null, r.obs === null ? "" : fmt2(r.obs)), el("span", null, r.exp === null ? "" : fmt2(r.exp)), el("span", null, r.excess_m === null ? "" : fmt2(r.excess_m)), el("span", null, r.bm_id === null ? "" : fmt1(r.bm_id)), el("span", null, r.lam === null ? "" : fmt1(r.lam)), el("span", null, r.bm))); }
      $(".n", foot).textContent = `${fmtInt(rows.length)} rows · click a row for the sequence panel`; };
    body.addEventListener("scroll", paint); paint();
  }
  function download() {
    const D = App.D; const keep = maskedOnly; maskedOnly = false; build(); maskedOnly = keep;
    const lines = ["taxon\tlocus\tstate\tn_fail\tn_similar\tobserved_divergence\texpected_divergence\texcess_m\tbest_match\tbest_match_identity\tpartner\tlrt"];
    for (const r of rows) { const pt = D.num("partner", r.g, r.li); const bm = D.num("best_match", r.g, r.li); lines.push([r.taxon, r.locus, D.labelOf(r.code), r.n_fail, r.n_similar, r.obs ?? "", r.exp ?? "", r.excess_m ?? "", bm !== null && bm >= 0 ? D.taxa.names[bm] : "", r.bm_id ?? "", pt !== null && pt >= 0 ? D.taxa.names[pt] : "", r.lam ?? ""].join("\t")); }
    const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([lines.join("\n")], { type: "text/tab-separated-values" })); a.download = `${D.man.dataset || "report"}_cells.tsv`; a.click(); build();
  }
  return { render, download };
})();

// ------------------------------------------------------------------ URL hash: layout, ordering, viewport, open panels
let hashTimer = null;
function writeHash() {
  clearTimeout(hashTimer); hashTimer = setTimeout(() => {
    const D = App.D; if (!D) return; const V = Overview.V; const p = new URLSearchParams();
    p.set("d", D.man.input.sha256.slice(0, 10)); p.set("L", App.layout); p.set("o", `${App.order.taxaMode},${App.order.lociMode}`);
    p.set("v", [V.x0, V.y0, V.cw, V.ch].map(v => +v.toFixed(3)).join(","));
    const dc = Detail.cur(); if (dc && App.panels.detail) p.set(dc.kind === "taxon" ? "t" : "l", dc.i);
    const sc = Sequence.cur(); if (sc && App.panels.sequence) p.set("c", `${sc.g},${sc.li}`);
    history.replaceState(null, "", "#" + p.toString());
  }, 150);
}
function readHash() {
  const D = App.D; const h = location.hash.slice(1); if (!h) return; const p = new URLSearchParams(h);
  if (p.get("d") && p.get("d") !== D.man.input.sha256.slice(0, 10)) return;        // a view of another report
  if (p.get("L")) { App.layout = p.get("L"); App.panelsEl.className = "panels " + App.layout; }
  if (p.get("o")) { const [a, b] = p.get("o").split(","); const sel = App.panels.overview.el.querySelectorAll("select"); if (a) { Overview.setOrder("taxa", a); sel[0].value = a; } if (b) { Overview.setOrder("loci", b); sel[1].value = b; } }
  if (p.get("v")) { const [x0, y0, cw, ch] = p.get("v").split(",").map(Number); if ([x0, y0, cw, ch].every(Number.isFinite)) { Overview.V.x0 = x0; Overview.V.y0 = y0; Overview.V.cw = cw; Overview.V.ch = ch; Overview.labelSpace(); App._viewFromHash = true; } }
  if (p.get("t")) Detail.open("taxon", +p.get("t")); else if (p.get("l")) Detail.open("locus", +p.get("l"));
  if (p.get("c")) { const [g, li] = p.get("c").split(",").map(Number); Sequence.open(g, li); }
}
// the overview's render fits the view unless the hash set one
const _render = Overview.render;
Overview.render = function () { if (App._viewFromHash) { App._viewFromHash = false; Overview.redraw(); const D = App.D; /* tiles and lists still need drawing */ } _render(); };

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot); else boot();
})();
