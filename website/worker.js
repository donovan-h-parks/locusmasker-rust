/* The worker that runs locus_outliers.py in Pyodide. Messages in: {type:"run", files:{alignment:{name,bytes}, partitions?:{name,bytes}},
   args:[...], tool:"<python source>"}. Messages out: {type:"line", text} for every stderr line (progress lines included),
   {type:"done", bundle:Uint8Array, outputs:{name:Uint8Array}}, {type:"error", text}. */
const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/";
importScripts(PYODIDE_URL + "pyodide.js");
let pyodideReady = null;
function ready() {
  if (!pyodideReady) pyodideReady = (async () => {
    postMessage({ type: "line", text: "loading the Python runtime…" });
    const py = await loadPyodide({ indexURL: PYODIDE_URL });
    postMessage({ type: "line", text: "loading numpy…" });
    await py.loadPackage("numpy");
    py.setStderr({ batched: text => postMessage({ type: "line", text }) });
    py.setStdout({ batched: text => postMessage({ type: "line", text }) });
    return py;
  })();
  return pyodideReady;
}
onmessage = async ev => {
  const msg = ev.data; if (msg.type !== "run") return;
  try {
    const py = await ready();
    const FS = py.FS;
    try { FS.mkdir("/work"); } catch (e) {}
    try { FS.mkdir("/work/out"); } catch (e) {}
    for (const f of FS.readdir("/work/out")) if (f !== "." && f !== "..") FS.unlink("/work/out/" + f);
    FS.writeFile("/work/locus_outliers.py", msg.tool);
    FS.writeFile("/work/" + msg.files.alignment.name, msg.files.alignment.bytes);
    if (msg.files.partitions) FS.writeFile("/work/" + msg.files.partitions.name, msg.files.partitions.bytes);
    const argv = ["locus_outliers.py", "--aln", "/work/" + msg.files.alignment.name, "--out", "/work/out", "--threads", "1", "--progress", "lines", ...msg.args];
    if (msg.files.partitions) argv.push("--partitions", "/work/" + msg.files.partitions.name);
    py.globals.set("lm_argv", py.toPy(argv));
    const code = `
import sys, runpy
sys.argv = list(lm_argv)
_err = None
try:
    runpy.run_path("/work/locus_outliers.py", run_name="__main__")
except SystemExit as e:
    _err = ("exit", str(e))
except MemoryError as e:
    _err = ("memory", "the browser ran out of memory for this dataset")
_err
`;
    const err = await py.runPythonAsync(code);
    if (err) { const [kind, text] = err.toJs(); postMessage({ type: "error", text: kind === "memory" ? text : text || "the tool stopped" }); return; }
    const outputs = {};
    for (const f of FS.readdir("/work/out")) { if (f === "." || f === ".." || f === "report.html" || f === "pairs.npy" || f === "trends.tsv" || f.endsWith(".npy") && f !== "report.npz") continue; outputs[f] = FS.readFile("/work/out/" + f); }
    const bundle = outputs["report.npz"]; delete outputs["report.npz"];
    // free the inputs and outputs held in the worker's file system
    for (const f of FS.readdir("/work/out")) if (f !== "." && f !== "..") FS.unlink("/work/out/" + f);
    FS.unlink("/work/" + msg.files.alignment.name);
    postMessage({ type: "done", bundle, outputs }, [bundle.buffer, ...Object.values(outputs).map(u => u.buffer)]);
  } catch (e) {
    const text = String(e && e.message || e);
    postMessage({ type: "error", text: /memory|OOM|allocat/i.test(text) ? "the browser ran out of memory for this dataset: " + text.split("\n")[0] : text });
  }
};
