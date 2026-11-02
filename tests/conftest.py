"""Shared fixtures: the tool as a module, the fixture datasets, an in-process runner with a per-session cache, and a
report loader. The tool is imported through TOOL_MODULE in one place, so a rename is a one-line change here."""
import csv, fcntl, hashlib, io, json, os, shutil, subprocess, sys, zipfile
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ.setdefault(_v, "1")   # before numpy loads, as the tool does: one BLAS thread per process
os.environ.setdefault("NUMBA_NUM_THREADS", "8")      # the compiled pair stage may use up to 8 threads in the suite, whichever test builds it first
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL_MODULE = "locus_outliers"
TOOL_PATH = os.path.join(ROOT, TOOL_MODULE + ".py")
DATA = os.path.join(ROOT, "tests", "data")
PY = sys.executable
if ROOT not in sys.path: sys.path.insert(0, ROOT)
lo = __import__(TOOL_MODULE)
CONSTANTS = ("DENSE_CELLS_MAX", "BEST_MATCH_ALL_MAX_TAXA", "BEST_MATCH_MIN_SITES", "REPORT_SEQ_ALL_MAX", "PASS1_FLOPS_PER_THREAD", "DISPERSION_PAIRS_PER_WORKER", "COMPARATOR_SECONDS_PER_WORKER",
             "TEST_BATCH_COMPARISONS", "TEST_BATCH_MIN", "PAIRS_PER_THREAD", "LSET_ELEMENTS", "NUMBA_PAIR_LOCI_MIN", "NUMBA_PAIR_SLAB", "PAIR_BLOCK", "SWEEP_BLOCK_MIN", "SWEEP_BLOCKS", "SWEEP_BLOCK_FLOOR", "CANDIDATES", "FULL_M_MAX_TAXA", "BAND_SLOTS", "LSET_CHUNK", "MIN_LOCI_PER_TAXON", "MIN_LOCI_FOR_M", "MIN_FIT_PAIRS", "PHI_MIN_E", "PHI_N0", "HALF", "_mem_available", "_cpus")
PRISTINE = {k: getattr(lo, k) for k in CONSTANTS}

try:
    from hypothesis import settings, HealthCheck
    settings.register_profile("suite", max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large])
    settings.register_profile("quick", max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large])
    settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "suite"))
except ImportError:
    pass


# ------------------------------------------------------------------ the fixture datasets
SYNTH = ["dna_small", "aa_small", "dna_clean", "aa_clean", "dna_dense", "dna_rc", "dna_blocks"]
REAL = ["cannon_aa", "near_dna"]
SIM = ["aa_random_trees"]        # simulated by AliSim (tests/data/make_random_trees.py); nothing planted, so every mask is a false positive


class Dataset:
    def __init__(self, name):
        self.name = name
        if name in SYNTH or name in SIM:
            self.dir = os.path.join(DATA, "synth" if name in SYNTH else "sim", name); self.truth = json.load(open(os.path.join(self.dir, "truth.json")))
            self.dtype = self.truth["datatype"]; self.aln = os.path.join(self.dir, "alignment.fasta"); self.parts = os.path.join(self.dir, "partitions.txt")
            self.events = self.truth["events"]; self.ntaxa = self.truth["taxa"]; self.nloci = self.truth["loci"]
        else:
            self.dir = os.path.join(DATA, "real", name); meta = json.load(open(os.path.join(self.dir, "meta.json")))
            self.dtype = meta["datatype"]; self.aln = os.path.join(self.dir, "alignment.nex"); self.parts = None; self.events = []
            self.ntaxa = meta["taxa"]; self.nloci = meta["loci"]; self.truth = meta

    def file(self, fmt):
        """an input file of the dataset in a given format (dna_small is written in every format)"""
        p = {"fasta": "alignment.fasta", "nexus": "alignment.nex", "nexus_quirky": "alignment_quirky.nex", "phylip": "alignment.phy", "phylip_interleaved": "alignment_interleaved.phy"}[fmt]
        path = os.path.join(self.dir, p)
        if not os.path.exists(path): pytest.skip(f"{self.name} is not written in {fmt}")
        return path

    def ev(self, kind): return [e for e in self.events if e["kind"] == kind]


def dataset(name): return Dataset(name)


# ------------------------------------------------------------------ running the tool
class Report:
    """report.npz opened for the tests: members as arrays, the manifest and the two JSON tables, and cell lookups"""
    def __init__(self, path):
        self.path = path; self.z = zipfile.ZipFile(path); self.names = self.z.namelist()
        self.man = json.loads(self.z.read("manifest.json")); self.taxa = json.loads(self.z.read("taxa.json")); self.loci = json.loads(self.z.read("loci.json"))
        self.n = len(self.taxa["names"]); self.nloc = len(self.loci["names"]); self.K = self.man["options"]["K_eff"]; self.M = self.man["options"]["M_eff"]
        self._cache = {}
        self.dense = self.man["carried"]["numeric_columns"] == "dense"
        if not self.dense:
            sl, st = self.arr("cells/sparse_locus.npy"), self.arr("cells/sparse_taxon.npy"); self.sparse = {(int(t), int(l)): i for i, (t, l) in enumerate(zip(st, sl))}

    def arr(self, member):
        if member not in self._cache: self._cache[member] = np.load(io.BytesIO(self.z.read(member)))
        return self._cache[member]

    def has(self, member): return member in self.names
    def js(self, member): return json.loads(self.z.read(member))
    def tx(self, name): return self.taxa["names"].index(name)
    def lx(self, name): return self.loci["names"].index(name)

    def cell(self, col, g, li):
        """a per-cell column at (taxon, locus) indices; None where a sparse report carries nothing"""
        a = self.arr(f"cells/{col}.npy")
        if self.dense or col in ("n_comp", "n_fail", "n_similar", "side", "focal_sites"): return a[g, li]
        i = self.sparse.get((g, li)); return None if i is None else a[i]

    def side(self, g, li): return int(self.arr("cells/side.npy")[g, li])
    def comps(self, g, li):
        c = self.cell("comparators", g, li); return [] if c is None else [int(x) for x in c if x >= 0]

    def similar(self, li): return self.js(f"similar/{li}.json") if self.has(f"similar/{li}.json") else []


class Run:
    def __init__(self, out, log, rc=0):
        self.out = out; self.log = log; self.rc = rc; self._rep = None
    @property
    def masked(self):
        p = os.path.join(self.out, "masked.csv"); return list(csv.DictReader(open(p, newline=""))) if os.path.exists(p) else []
    def masked_at(self, taxon, locus):
        for r in self.masked:
            if r["taxon"] == taxon and r["locus"] == locus: return r
        return None
    @property
    def report(self):
        if self._rep is None: self._rep = Report(os.path.join(self.out, "report.npz"))
        return self._rep
    def file(self, name): return os.path.join(self.out, name)
    def exists(self, name): return os.path.exists(os.path.join(self.out, name))


def run_tool(aln, out, dtype, parts=None, extra=(), progress="none"):
    """the tool in this process (so coverage sees it): main() with sys.argv set, stderr captured, the module's progress
    and log state restored afterwards. SystemExit is caught and its code returned."""
    argv = [TOOL_PATH, "--aln", aln, "--out", out, "--datatype", dtype, "--progress", progress] + (["--partitions", parts] if parts else []) + list(extra)
    old_argv, old_err = sys.argv, sys.stderr; sys.argv = argv; err = io.StringIO(); sys.stderr = err
    old_prog, old_bar, old_fh = lo.PROG, lo.LOG.bar, lo.LOG.fh; lo.LOG.real = err; rc = 0
    try:
        try: lo.main()
        except SystemExit as e: rc = e.code if isinstance(e.code, int) else 1; err.write(str(e.code) + "\n") if e.code and not isinstance(e.code, int) else None
    finally:
        sys.argv, sys.stderr = old_argv, old_err; lo.PROG = old_prog; lo.LOG.bar = old_bar; lo.LOG.real = sys.stderr
        if lo.LOG.fh is not None: lo.LOG.close()
        lo.LOG.fh = old_fh; lo.MEM_BUDGET = None; lo.MEM_RESIDENT_EST = 0; lo.USE_NUMBA = True
    return Run(out, err.getvalue(), rc)


def run_cli(args, cwd=None, timeout=600):
    """the tool as a separate process: the command line a user types"""
    return subprocess.run([PY, TOOL_PATH] + list(args), capture_output=True, text=True, cwd=cwd, timeout=timeout)


@pytest.fixture(scope="session")
def shared_tmp(tmp_path_factory):
    """one directory shared by every xdist worker of the session (the parent of the per-worker basetemp)"""
    base = tmp_path_factory.getbasetemp()
    root = base.parent if base.name.startswith("popen-") else base
    d = root / "lm_shared"; d.mkdir(exist_ok=True); return d


@pytest.fixture(scope="session")
def runs(shared_tmp):
    """the standard run of each fixture dataset, made once per session (under a file lock, so parallel workers share it)"""
    cache = {}
    def get(name, extra=(), key=None):
        key = key or name
        if key in cache: return cache[key]
        out = shared_tmp / "runs" / key; out.mkdir(parents=True, exist_ok=True); lock = shared_tmp / "runs" / (key + ".lock")
        with open(lock, "w") as lk:
            fcntl.flock(lk, fcntl.LOCK_EX)
            try:
                if not (out / "done").exists():
                    patched = {k: getattr(lo, k) for k in CONSTANTS}
                    for k, v in PRISTINE.items(): setattr(lo, k, v)                 # never build the shared run under a test's monkeypatch
                    try: ds = dataset(name); r = run_tool(ds.aln, str(out), ds.dtype, ds.parts, extra)
                    finally:
                        for k, v in patched.items(): setattr(lo, k, v)
                    if r.rc != 0: raise RuntimeError(f"the standard run of {name} failed:\n{r.log}")
                    (out / "log.txt").write_text(r.log); (out / "done").write_text("ok")
            finally: fcntl.flock(lk, fcntl.LOCK_UN)
        cache[key] = Run(str(out), (out / "log.txt").read_text()); return cache[key]
    return get


@pytest.fixture
def tmp_out(tmp_path):
    return str(tmp_path / "out")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def npz_hashes(path):
    z = zipfile.ZipFile(path); return {k: hashlib.sha256(z.read(k)).hexdigest() for k in sorted(z.namelist()) if k != "manifest.json"}
