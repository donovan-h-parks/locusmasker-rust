"""Score the planted-swap benchmark (bench/results) : detection vs relative divergence, FPs."""
import json, glob, os, numpy as np
res = [json.load(open(f)) for f in glob.glob("bench/results/*.json")]
meta = {d: json.load(open(f"bench/{d}/meta.json")) for d in os.listdir("bench") if os.path.exists(f"bench/{d}/meta.json")}
nbh = {}
for d in meta:
    try:
        M = np.load(f"bench/{d}/baseline/m.npy"); taxa = [t.strip() for t in open(f"bench/{d}/baseline/taxa.txt")]
        Ms = np.sort(np.where(np.isfinite(M), M, np.inf), axis=1); nbh[d] = dict(zip(taxa, np.median(Ms[:, :5], axis=1)))
    except FileNotFoundError: pass
rel, det, side, cells, fpc, fpa, lost, errs = [], [], [], 0, 0, 0, 0, 0
for r in res:
    if r["error"]: errs += 1; continue
    cells += r["n_cells"]; fpc += len(r["fp_clean"]); fpa += len(r["fp_adjacent"]); lost += len(r["lost"])
    for ev in r["events"]:
        for s in ("a", "b"):
            st = ev.get(f"det_{s}", "untested")
            if st == "untested": continue
            nb = nbh.get(r["dataset"], {}).get(ev[s], np.nan)
            if not np.isfinite(nb): continue
            rel.append(ev["m_ab"] / max(nb, 0.05)); det.append(st in ("flag", "amb", "sim")); side.append(st)
rel, det = np.array(rel), np.array(det)
print(f"{len(res)} runs ({errs} errors), {cells:,} cells tested, {len(rel)} planted cells scored")
print(f"FP clean {fpc} ({100*fpc/max(cells,1):.4f}%), adjacent {fpa}, baseline masks lost {lost} ({lost/max(len(res)-errs,1):.1f}/run)")
for lo, hi in [(0, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 10), (10, 30), (30, 1e9)]:
    s = (rel >= lo) & (rel < hi)
    if s.sum(): print(f"  rel {lo:>4}-{hi if hi < 1e9 else 'inf':<4} n={int(s.sum()):5d} detected {100*det[s].mean():5.1f}%")
s = rel >= 3
print(f"rel >= 3: {100*det[s].mean():.1f}% of {int(s.sum())}; by side among detected: {dict(zip(*np.unique(np.array(side)[s & det], return_counts=True)))}")
