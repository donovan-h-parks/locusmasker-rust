/* The website shell: the run panel (inputs, options, size estimate, progress, downloads) around the shared viewer.
   The analysis itself is locus_outliers.py running unchanged in a Pyodide worker (worker.js). */
(() => {
"use strict";
const $ = s => document.querySelector(s);
const TOOL_URL = "../locus_outliers.py", VIEWER_DIR = "../viewer/";
const WASM_GFLOPS = 2.4e9, HEAP_LIMIT = 1.5e9, TIME_LIMIT = 900;        // measured: Pyodide pass one at about 2.4 GFLOPS; Wu 115 s and 1.15 GB
const STAGES = [["read", "reading"], ["counts", "counting sites"], ["pairs0", "provisional divergences"], ["rates", "locus rates"], ["pairs", "pair divergences"], ["comparators", "choosing comparators"], ["dispersion", "dispersion"], ["test", "testing loci"], ["write", "writing"]];
const fmtInt = n => n.toLocaleString("en"), fmtTime = s => s < 90 ? `${Math.round(s)} s` : s < 5400 ? `${(s / 60).toFixed(s < 600 ? 1 : 0)} min` : `${(s / 3600).toFixed(1)} h`;
const fmtBytes = b => b < 1e9 ? `${Math.round(b / 1e6)} MB` : `${(b / 1e9).toFixed(1)} GB`;
let header = null, worker = null, outputs = null, bundleBytes = null, live = null, t0 = 0;

// ---- the input's header: taxa, sites, loci, declared datatype
async function readSlice(file, start, end) { return new TextDecoder("latin1").decode(await file.slice(start, end).arrayBuffer()); }
async function scanHeader(file) {
  const head = await readSlice(file, 0, Math.min(file.size, 1 << 18)); const h = { format: "?", taxa: null, sites: null, loci: null, datatype: null };
  if (/^\s*#NEXUS/i.test(head)) {
    h.format = "nexus";
    const m = /dimensions\s+[^;]*ntax\s*=\s*(\d+)[^;]*nchar\s*=\s*(\d+)/i.exec(head) || /dimensions\s+[^;]*nchar\s*=\s*(\d+)[^;]*ntax\s*=\s*(\d+)/i.exec(head);
    if (m) { if (/ntax\s*=\s*\d+[^;]*nchar/i.test(m[0])) { h.taxa = +m[1]; h.sites = +m[2]; } else { h.sites = +m[1]; h.taxa = +m[2]; } }
    const dt = /datatype\s*=\s*(\w+)/i.exec(head); if (dt) h.datatype = /prot/i.test(dt[1]) ? "aa" : /dna|nuc|rna/i.test(dt[1]) ? "dna" : null;
    const tail = await readSlice(file, Math.max(0, file.size - (4 << 20)), file.size); const cs = tail.match(/^\s*charset\s/gim); h.loci = cs ? cs.length : 0;
  } else if (/^\s*>/.test(head)) {
    h.format = "fasta"; let count = 0, first = null, firstLen = 0, inFirst = false; const CH = 8 << 20;
    for (let off = 0; off < file.size; off += CH) {
      const txt = await readSlice(file, off, Math.min(file.size, off + CH));
      for (let i = 0; i < txt.length; i++) { const c = txt[i]; if (c === ">") { count++; inFirst = count === 1; if (count === 2) break; } else if (inFirst && c !== "\n" && c !== "\r") { if (first === null) { if (c === "\n") first = 0; } } }
      if (count >= 2) break;
    }
    // the first record's length: from the first chunk
    const m = /^>[^\n]*\n([\s\S]*?)(?=\n>|$)/.exec(head); if (m) firstLen = m[1].replace(/\s/g, "").length;
    h.taxa = count; h.sites = firstLen || null;
  } else {
    const m = /^\s*(\d+)\s+(\d+)/.exec(head); if (m) { h.format = "phylip"; h.taxa = +m[1]; h.sites = +m[2]; }
  }
  return h;
}
async function countPartitions(file) { const txt = await readSlice(file, 0, Math.min(file.size, 8 << 20)); return (txt.match(/=/g) || []).length; }

// ---- the estimate, from the formulas the tool itself uses
function estimate(h, dtype, inputBytes) {
  if (!h.taxa || !h.sites) return null;
  const n = h.taxa, L = h.sites, nloc = h.loci || 1; const kbar = dtype === "aa" ? 3.3 : 2.0;
  const flops = 2 * n * n * L * (1 + kbar); const cells = 0.9 * n * nloc;
  const tWasm = flops / WASM_GFLOPS + 6e-5 * cells + 2e-9 * nloc * n * n + 3;
  const heap = 1.6 * (2 * n * L + 2 * nloc * n * (n + 1) + 48 * n * n) + 1e8 + inputBytes;
  const t1 = tWasm / 5; return { tWasm, heap, native: { 1: t1, 16: t1 / 8, 32: t1 / 12, 64: t1 / 16 }, ok: heap <= HEAP_LIMIT && tWasm <= TIME_LIMIT };
}
async function updateEstimate() {
  const f = $("#aln").files[0]; const box = $("#estimate"); if (!f) { box.style.display = "none"; $("#go").disabled = true; return; }
  header = await scanHeader(f); const parts = $("#parts").files[0]; if (parts) header.loci = await countPartitions(parts);
  if (header.datatype) document.querySelector(`input[name=dt][value=${header.datatype}]`).checked = true;
  $("#alnInfo").textContent = header.taxa ? `${header.format}: ${fmtInt(header.taxa)} taxa, ${header.sites ? fmtInt(header.sites) + " sites" : "?"}, ${header.loci ? fmtInt(header.loci) + " loci" : "no loci found: a partition file is needed"}` : "format not recognised";
  const dtype = document.querySelector("input[name=dt]:checked").value; const e = estimate(header, dtype, f.size);
  box.style.display = ""; box.className = "estimate" + (e && !e.ok ? " refuse" : "");
  if (!e) { box.textContent = "the header could not be read; the run will tell you if the file is not understood"; $("#go").disabled = !header.taxa; return; }
  if (e.ok) { box.textContent = `estimate for this browser: about ${fmtTime(e.tWasm)} and ${fmtBytes(e.heap)} of memory, single-threaded`; $("#go").disabled = false; }
  else {
    const cmd = `python locus_outliers.py --aln ${f.name}${parts ? " --partitions " + parts.name : ""} --datatype ${dtype} --out ${f.name.replace(/\.[^.]+$/, "")}_masked --threads 16`;
    box.replaceChildren(document.createTextNode(`This dataset is too large to run comfortably in a browser tab (about ${fmtTime(e.tWasm)} and ${fmtBytes(e.heap)} of memory here). Please run the command-line tool, which writes the same report:`),
      Object.assign(document.createElement("pre"), { textContent: cmd }),
      Object.assign(document.createElement("div"), { textContent: `It needs about ${fmtBytes(e.heap / 1.4)} of memory and takes roughly ${fmtTime(e.native[1])} on 1 thread, ${fmtTime(e.native[16])} on 16, ${fmtTime(e.native[32])} on 32 and ${fmtTime(e.native[64])} on 64 threads (estimates from the measured scaling; the first pass dominates). Then open its report.html here or by double-click.` }));
    $("#go").disabled = true;
  }
}
for (const id of ["#aln", "#parts"]) $(id).addEventListener("change", updateEstimate);
for (const r of document.querySelectorAll("input[name=dt]")) r.addEventListener("change", updateEstimate);

// ---- running
function stageRow(key, label) { const row = document.createElement("div"); row.className = "stage"; row.dataset.stage = key; row.innerHTML = `<span>${label}</span><div class="bar"><i></i></div><span class="lab"></span>`; return row; }
async function run() {
  const f = $("#aln").files[0]; const parts = $("#parts").files[0]; if (!f) return;
  const dtype = document.querySelector("input[name=dt]:checked").value;
  const args = ["--datatype", dtype, "--n-comp", $("#ncomp").value, "--fit-overlap", $("#fit").value, "--alpha", $("#alpha").value];
  if ($("#simonly").checked) args.push("--mask-similar-only"); if ($("#divonly").checked) args.push("--mask-different-only");
  args.push("--name", f.name.replace(/\.[^.]+$/, ""));
  $("#go").disabled = true; $("#cancel").style.display = ""; $("#err").textContent = ""; $("#status").textContent = "starting…";
  const prog = $("#progress"); prog.replaceChildren(...STAGES.map(([k, l]) => stageRow(k, l))); prog.style.display = "block";
  $("#log").style.display = ""; $("#log").textContent = ""; logLines = []; $("#live").style.display = "none"; live = null; t0 = performance.now();
  const tool = await (await fetch(TOOL_URL)).text();
  const files = { alignment: { name: f.name, bytes: new Uint8Array(await f.arrayBuffer()) } };
  if (parts) files.partitions = { name: parts.name, bytes: new Uint8Array(await parts.arrayBuffer()) };
  worker = new Worker("worker.js");
  worker.onmessage = ev => {
    const m = ev.data;
    if (m.type === "line") { for (const line of m.text.split("\n")) if (line.length) onLine(line); }
    else if (m.type === "error") { $("#err").textContent = m.text; $("#status").textContent = "stopped"; $("#cancel").style.display = "none"; $("#go").disabled = false; worker.terminate(); worker = null; }
    else if (m.type === "done") { outputs = m.outputs; bundleBytes = m.bundle; $("#status").textContent = `done in ${fmtTime((performance.now() - t0) / 1000)}`; worker.terminate(); worker = null; showResult(); }
  };
  worker.onerror = e => { $("#err").textContent = "the worker failed: " + e.message; $("#cancel").style.display = "none"; $("#go").disabled = false; };
  worker.postMessage({ type: "run", files, args, tool }, [files.alignment.bytes.buffer, ...(files.partitions ? [files.partitions.bytes.buffer] : [])]);
}
$("#go").addEventListener("click", run);
$("#cancel").addEventListener("click", () => { if (worker) worker.terminate(); worker = null; $("#status").textContent = "cancelled"; $("#cancel").style.display = "none"; $("#go").disabled = false; });

function onLine(line) {
  if (line.startsWith("progress\t")) { const [, stage, done, total, elapsed, extra] = line.split("\t"); const row = $(`.stage[data-stage="${stage}"]`); if (!row) return;
    const d = +done, t = +total; row.querySelector("i").style.width = `${100 * d / t}%`; row.querySelector(".lab").textContent = `${fmtInt(d)} / ${fmtInt(t)} · ${fmtTime(+elapsed)}${extra ? " · " + extra : ""}`; if (d >= t) row.classList.add("done");
    $("#status").textContent = `${STAGES.find(s => s[0] === stage)?.[1] || stage}…`; return; }
  if (line.startsWith("event\tlocus\t")) { onLocus(line.split("\t")); return; }
  const m = /^(\d+) taxa, (\d+) sites, (\d+) loci/.exec(line); if (m) initLive(+m[1], +m[3]);
  logLines.push(line); if (logLines.length > 400) logLines.splice(0, logLines.length - 300);
  if (!logTimer) logTimer = setTimeout(() => { logTimer = null; const log = $("#log"); log.textContent = logLines.join("\n"); log.scrollTop = log.scrollHeight; }, 250);
}
let logLines = [], logTimer = null;
// ---- the live matrix: painted as each locus's result arrives
function initLive(n, nloc) {
  const cv = $("#liveCanvas"); const s = Math.max(1, Math.ceil(Math.sqrt(n * nloc / 4e6))); const W = Math.ceil(nloc / s), H = Math.ceil(n / s);
  cv.width = W; cv.height = H; const ctx = cv.getContext("2d"); const img = ctx.createImageData(W, H);
  const surf = getComputedStyle(document.documentElement).getPropertyValue("--surface").trim(); const rgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
  const bg = rgb(surf); for (let i = 0; i < W * H; i++) img.data.set([bg[0], bg[1], bg[2], 255], i * 4);
  live = { n, nloc, s, W, H, ctx, img, grey: rgb(getComputedStyle(document.documentElement).getPropertyValue("--grey-1").trim()), red: rgb(getComputedStyle(document.documentElement).getPropertyValue("--red-2").trim()), blue: rgb(getComputedStyle(document.documentElement).getPropertyValue("--blue-2").trim()), sev: new Uint8Array(W * H), div: 0, sim: 0, loci: 0 };
  $("#live").style.display = ""; cv.style.height = Math.min(300, Math.max(60, H * 2)) + "px"; ctx.putImageData(img, 0, 0);
}
function onLocus(f) {
  if (!live) return; const li = +f[2]; const tested = f[4] || ""; const div = f[5] ? f[5].split(",").filter(Boolean).map(Number) : []; const sim = f[6] ? f[6].split(",").filter(Boolean).map(Number) : [];
  const x = Math.floor(li / live.s); const paint = (g, col, sv) => { const y = Math.floor(g / live.s); const i = y * live.W + x; if (live.sev[i] >= sv) return; live.sev[i] = sv; live.img.data.set([col[0], col[1], col[2], 255], i * 4); };
  for (let byte = 0; byte < tested.length / 2; byte++) { const v = parseInt(tested.slice(byte * 2, byte * 2 + 2), 16); for (let b = 0; b < 8; b++) if (v & (0x80 >> b)) { const g = byte * 8 + b; if (g < live.n) paint(g, live.grey, 1); } }
  for (const g of div) paint(g, live.red, 2); for (const g of sim) paint(g, live.blue, 3);
  live.div += div.length; live.sim += sim.length; live.loci++;
  $("#tDiv").textContent = fmtInt(live.div); $("#tSim").textContent = fmtInt(live.sim); $("#tLoci").textContent = fmtInt(live.loci);
  if (live.loci % 20 === 0 || live.loci === live.nloc) live.ctx.putImageData(live.img, 0, 0);
}

// ---- the result: the viewer, with a downloads bar
async function showResult() {
  $("#run").style.display = "none"; const app = $("#app"); app.style.display = "";
  const bar = document.createElement("div"); bar.className = "downloads"; document.body.insertBefore(bar, app);
  const add = (label, name, bytes, type) => { const b = document.createElement("button"); b.textContent = label; b.title = name; b.onclick = () => { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([bytes], { type })); a.download = name; a.click(); }; bar.append(b); };
  const k = document.createElement("span"); k.className = "k"; k.textContent = "downloads:"; bar.append(k);
  const masked = Object.keys(outputs).find(n => n.includes(".masked.")); if (masked) add("masked alignment", masked, outputs[masked], "text/plain");
  for (const n of ["masked.csv", "cells.tsv", "per_locus.tsv", "similar.tsv", "run.log"]) if (outputs[n]) add(n, n, outputs[n], "text/plain");
  add("report.npz", "report.npz", bundleBytes, "application/zip");
  const rb = document.createElement("button"); rb.textContent = "report.html"; rb.className = "primary"; rb.onclick = async () => { rb.disabled = true; rb.textContent = "assembling…"; const html = await assembleReport(); const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([html], { type: "text/html" })); a.download = "report.html"; a.click(); rb.disabled = false; rb.textContent = "report.html"; }; bar.append(rb);
  const again = document.createElement("button"); again.textContent = "new run"; again.onclick = () => location.reload(); bar.append(Object.assign(document.createElement("span"), { style: "flex:1" }), again);
  window.LM_OUTPUTS = outputs;
  await window.LM.start(bundleBytes.slice());
}
function b64(bytes) { let s = ""; const CH = 1 << 15; for (let i = 0; i < bytes.length; i += CH) s += String.fromCharCode.apply(null, bytes.subarray(i, i + CH)); return btoa(s); }
async function assembleReport() {                          // the same page the tool writes: the template with its scripts inline and the bundle in its tag
  const get = async p => (await fetch(VIEWER_DIR + p)).text();
  let html = await get("viewer.html"); const data = b64(bundleBytes);
  html = html.replace("@@" + "BUNDLE@@", data).replace("@@" + "SIDECAR@@", "");
  for (const m of Array.from(html.matchAll(/<script src="([^"]+)"><\/script>/g))) html = html.replace(m[0], "<script>\n" + await get(m[1]) + "\n</script>");
  for (const m of Array.from(html.matchAll(/<link rel="stylesheet" href="([^"]+)">/g))) html = html.replace(m[0], "<style>\n" + await get(m[1]) + "\n</style>");
  return html;
}
// ---- opening a finished report
$("#openReport").addEventListener("change", async ev => {
  const f = ev.target.files[0]; if (!f) return; const buf = new Uint8Array(await f.arrayBuffer());
  let bytes = null;
  if (f.name.endsWith(".npz")) bytes = buf; else { const text = new TextDecoder().decode(buf); const m = /id="lm-bundle"[^>]*>([^<]*)</.exec(text) || /LM_BUNDLE_B64\s*=\s*"([^"]+)"/.exec(text); if (m && m[1].trim()) { const b = m[1].replace(/\s+/g, ""); const bin = atob(b); bytes = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i); } }
  if (!bytes) { $("#err").textContent = "no report data found in " + f.name; return; }
  $("#run").style.display = "none"; $("#app").style.display = ""; bundleBytes = bytes; await window.LM.start(bytes);
});
})();
