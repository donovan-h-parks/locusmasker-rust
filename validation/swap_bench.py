"""Empirical validation of locus_outliers --test z by random locus overriding.

For each benchmark dataset: convert to FASTA+partitions once, run a baseline, then in each
replicate plant swap events (pick two random taxa, exchange their sequence in one locus) with
partners drawn log-uniformly across the divergence spectrum, rerun the tool, and score:
  - detection per planted cell (both sides of each swap), with the pair divergence m(a,b)
  - false positives: cells flagged in the run but neither planted nor flagged in the baseline
Results: one JSON per job in bench/results/. Baselines keep their masked.csv (natural flags).

Usage: python swap_bench.py pilot|full
"""
import json, os, subprocess, sys, time, zlib
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = f"{REPO}/test_data"
HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = f"{HERE}/bench"
TOOL = f"{REPO}/locus_outliers.py"
HEAVY = {"Wu_2018_dna", "Wu_2018_aa", "Shen_2018", "Jarvis_2015", "Irisarri_2017", "Misof_2014"}
sys.path.insert(0, REPO)
import locus_outliers as lo


def compress_ranges(sites):
    """0-based sorted site array -> RAxML range string (1-based, with \\step runs)"""
    runs, i = [], 0
    while i < len(sites):
        j = i + 1
        step = sites[j] - sites[i] if j < len(sites) else 1
        while j + 1 < len(sites) and sites[j + 1] - sites[j] == step: j += 1
        if j == i + 1 and step != 1:                       # lone pair with odd step: emit singly
            runs.append(f"{sites[i] + 1}"); i += 1; continue
        a, b = sites[i], sites[j if j < len(sites) else i]
        runs.append(f"{a + 1}" if a == b else f"{a + 1}-{b + 1}" + (f"\\{step}" if step != 1 else ""))
        i = j + 1
    return ", ".join(runs)


def prep(ds):
    out = f"{BENCH}/{ds}"; os.makedirs(out, exist_ok=True)
    if os.path.exists(f"{out}/meta.json"): return ds, "cached"
    try:
        names, seqs, loci, fmt, dtype, _ = lo.read_input(f"{DATA}/{ds}/alignment.nex", None)
    except SystemExit as e:
        return ds, f"unusable: {e}"
    if len(loci) < 3: return ds, f"unusable: {len(loci)} loci"
    with open(f"{out}/base.fasta", "w") as fh:
        for t, s in zip(names, seqs): fh.write(f">{t}\n{s}\n")
    with open(f"{out}/base.parts.txt", "w") as fh:
        for name, sites in loci: fh.write(f"DNA, {name} = {compress_ranges(sites)}\n")
    json.dump({"ntax": len(names), "nchar": len(seqs[0]), "nloci": len(loci), "dtype": dtype},
              open(f"{out}/meta.json", "w"))
    return ds, "ok"


def run_tool(aln, parts, out, jobdir, dtype):
    r = subprocess.run([sys.executable, TOOL, "--aln", aln, "--partitions", parts, "--out", out,
                        "--datatype", "aa" if dtype == "protein" else "dna"],
                       capture_output=True, text=True, timeout=10800, cwd=jobdir)
    if r.returncode != 0: raise RuntimeError(r.stderr.strip().splitlines()[-1] if r.stderr else "failed")


def read_flags(out):
    """flags = cells masked on either side (from masked.csv); cells = cells.tsv rows by header"""
    import csv as _csv
    cells = {(r["locus"], r["taxon"]): r for r in _csv.DictReader(open(f"{out}/cells.tsv"), delimiter="\t")}
    flags = {(r["locus"], r["taxon"]) for r in _csv.DictReader(open(f"{out}/masked.csv"))}
    return flags, cells


def baseline(ds):
    d = f"{BENCH}/{ds}"; out = f"{d}/baseline"
    if os.path.exists(f"{out}/.done"): return ds, "cached"
    t0 = time.time()
    run_tool(f"{d}/base.fasta", f"{d}/base.parts.txt", out, d, json.load(open(f"{d}/meta.json"))["dtype"])
    flags, _ = read_flags(out)
    tested = [f.split("\t")[0] for f in open(f"{out}/per_locus.tsv")
              if f.rstrip().endswith("\ttested")]
    json.dump({"flags": sorted(map(list, flags)), "tested_loci": tested,
               "runtime": time.time() - t0}, open(f"{d}/baseline.json", "w"))
    for f in ("pairs.npy", "cells.tsv", "trends.tsv"):                 # keep the small files only
        p = f"{out}/{f}"
        if os.path.getsize(p) > 50e6: os.remove(p)
    open(f"{out}/.done", "w").close()
    return ds, f"ok {time.time() - t0:.0f}s, {len(flags)} natural flags, {len(tested)} loci tested"


def swap_job(ds, rep):
    d = f"{BENCH}/{ds}"; res = f"{BENCH}/results/{ds}_rep{rep}.json"
    if os.path.exists(res): return ds, rep, "cached"
    t0 = time.time()
    meta = json.load(open(f"{d}/meta.json")); base = json.load(open(f"{d}/baseline.json"))
    names, seqs = lo.read_fasta(open(f"{d}/base.fasta").read())
    n, L = len(names), len(seqs[0])
    loci = lo.parse_raxml_partitions(f"{d}/base.parts.txt", L)
    tidx = {t: i for i, t in enumerate(names)}
    m = np.load(f"{d}/baseline/m.npy")
    border = [t.strip() for t in open(f"{d}/baseline/taxa.txt")]
    bidx = {t: i for i, t in enumerate(border)}
    alphabet = lo.DNA if meta["dtype"] == "dna" else lo.PROTEIN
    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    arr = np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L).copy()
    pres = lut[arr] >= 0
    tested = [(nm, s) for nm, s in loci if nm in set(base["tested_loci"])]
    rng = np.random.default_rng(zlib.crc32(f"{ds}_rep{rep}".encode()))     # reproducible across processes
    used, events, percap = set(), [], {}
    target = min(24, max(2, n // 6))
    for _ in range(40 * target):
        if len(events) >= target: break
        nm, sites = tested[rng.integers(len(tested))]
        npl = pres[:, sites].sum(1)
        cand = [i for i in range(n) if npl[i] >= 20 and i not in used]
        cap = max(1, int(0.02 * (npl >= 20).sum()))
        if len(cand) < 2 or percap.get(nm, 0) >= cap: continue
        a = int(rng.choice(cand))
        parts = [(m[bidx[names[a]], bidx[names[b]]], b) for b in cand if b != a]
        parts = sorted((p for p in parts if np.isfinite(p[0])))
        if not parts: continue
        j = min(int(len(parts) ** rng.random()) - 1, len(parts) - 1)
        m_ab, b = parts[max(j, 0)]
        if (pres[a, sites] & pres[b, sites]).sum() < 20: continue
        tmp = arr[a, sites].copy(); arr[a, sites] = arr[b, sites]; arr[b, sites] = tmp
        used.update((a, b)); percap[nm] = percap.get(nm, 0) + 1
        events.append({"locus": nm, "a": names[a], "b": names[b], "m_ab": float(m_ab),
                       "lsites": len(sites)})
    jobdir = f"{d}/rep{rep}"; os.makedirs(jobdir, exist_ok=True)
    mod = f"{jobdir}/mod.fasta"
    with open(mod, "w") as fh:
        for i, t in enumerate(names): fh.write(f">{t}\n{arr[i].tobytes().decode()}\n")
    out = f"{jobdir}/scan"
    err = None
    try:
        run_tool(mod, f"{d}/base.parts.txt", out, jobdir, meta["dtype"])
        flags, cells = read_flags(out)
    except Exception as e:
        err = str(e); flags, cells = set(), {}
    planted = set()
    if not err:
        for ev in events:
            for side in ("a", "b"):
                key = (ev["locus"], ev[side]); planted.add(key)
                row = cells.get(key)
                if row is None:
                    ev[f"det_{side}"] = "untested"; continue
                ev[f"det_{side}"] = "flag" if row["flag"] == "1" else ("amb" if row["ambiguous"] == "1" else ("sim" if key in flags else "miss"))
                ev[f"exc_{side}"], ev[f"nf_{side}"] = row["med_excess"], row["n_fail"]
        basefl = set(map(tuple, base["flags"]))
        fp = sorted(flags - planted - basefl)
        lost = sorted(basefl - flags)
        # collateral vs clean FPs: did the FP cell's comparator set contain a swapped taxon?
        evloc = {}
        for ev in events: evloc.setdefault(ev["locus"], set()).update((ev["a"], ev["b"]))
        comp = {tuple(k): set() for k in fp}
        if comp:
            P = np.load(f"{out}/pairs.npy"); taxa = open(f"{out}/taxa.txt").read().split()
            lnames = [r.split("\t", 1)[0] for r in open(f"{out}/per_locus.tsv").read().splitlines()[1:]]
            for q in P:
                k = (lnames[q["locus"]], taxa[q["taxon"]])
                if k in comp: comp[k].add(taxa[q["comparator"]])
            for ln in open(f"{out}/similar.tsv"):
                f = ln.split("\t", 3)
                for k in ((f[0], f[1]), (f[0], f[2])):
                    if k in comp: comp[k].add(f[2] if k[1] == f[1] else f[1])
        fp_adj = [list(k) for k in fp if comp[tuple(k)] & evloc.get(k[0], set())]
        fp_clean = [list(k) for k in fp if not (comp[tuple(k)] & evloc.get(k[0], set()))]
    else:
        fp, lost, fp_adj, fp_clean = [], [], [], []
    json.dump({"dataset": ds, "rep": rep, "error": err, "events": events,
               "fp_adjacent": fp_adj, "fp_clean": fp_clean, "lost": lost, "n_cells": len(cells),
               "runtime": time.time() - t0}, open(res, "w"))
    subprocess.run(["rm", "-rf", jobdir])
    return ds, rep, err or f"ok {time.time() - t0:.0f}s, {len(events)} events, {len(fp)} FP"


def main():
    os.makedirs(f"{BENCH}/results", exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = "4"
    pilot = sys.argv[1:] and sys.argv[1] == "pilot"
    datasets = (["Dornburg_2012", "Oaks_2011", "Wainwright_2012"] if pilot
                else sorted(d for d in os.listdir(DATA) if os.path.isdir(f"{DATA}/{d}")))
    with ProcessPoolExecutor(max_workers=16) as ex:
        usable = []
        for ds, st in ex.map(prep, datasets):
            print(f"prep {ds}: {st}", flush=True)
            if not st.startswith("unusable"): usable.append(ds)
        futs = {ex.submit(baseline, ds): ds for ds in usable}
        okbase = []
        for fu in as_completed(futs):
            try: ds, st = fu.result(); print(f"baseline {ds}: {st}", flush=True); okbase.append(ds)
            except Exception as e: print(f"baseline {futs[fu]}: FAILED {e}", flush=True)
        jobs = [(ds, r) for ds in sorted(okbase) for r in range((1 if pilot else 3) if (ds in HEAVY or pilot) else 8)]
        futs = {ex.submit(swap_job, ds, r): (ds, r) for ds, r in jobs}
        for fu in as_completed(futs):
            try: ds, r, st = fu.result(); print(f"swap {ds} rep{r}: {st}", flush=True)
            except Exception as e: print(f"swap {futs[fu]}: FAILED {e}", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
