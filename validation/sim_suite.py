"""Simulation suite for locus_outliers.py: random trees, planted truths, scoring.

Each config: a random bifurcating tree, 14 loci (lengths 200-3000, rates 0.3-3x, gamma site
rates, lognormal branch noise, optional fragment missingness), then planted events on disjoint
taxa: contaminants (random residues), swaps (random partner), near-copies (partner among the 5
patristically nearest; recipient recorded), a chimera (5 loci from one distant donor), and a
misalignment window (30% of a cell's sites shuffled). Runs locus_outliers.py (LO_TOOL overrides the
path; LO_BASE the working directory), scores masked.csv / similar.tsv against truth. Nothing is kept.

usage: python sim_suite.py [ncpu]
"""
import json, os, subprocess, sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.environ.get("LO_TOOL", os.path.join(os.path.dirname(HERE), "locus_outliers.py"))
BASE = os.environ.get("LO_BASE", f"{HERE}/simsuite")
PROTEIN = "ACDEFGHIKLMNPQRSTVWY"


def sim_dataset(cfg, rng):
    n, depth, A, pmiss = cfg["n"], cfg["depth"], cfg["alpha_size"], cfg["pmiss"]
    nloc = cfg.get("nloc", 14)
    lens = rng.integers(200, 3000, nloc)
    locrate = np.exp(rng.uniform(np.log(0.3), np.log(3), nloc))
    # random bifurcating tree: split a random leaf until n, exponential branch lengths
    parent, blen = {0: None}, {0: 0.0}
    leaves, nxt = [0], 1
    while len(leaves) < n:
        i = leaves.pop(int(rng.integers(len(leaves))))
        for _ in range(2):
            parent[nxt] = i; blen[nxt] = float(rng.exponential(depth / np.log2(n)))
            leaves.append(nxt); nxt += 1
    # patristic distances between leaves (for partner choice/reporting)
    def path_to_root(i):
        p = []
        while i is not None: p.append(i); i = parent[i]
        return p
    paths = {l: path_to_root(l) for l in leaves}
    depth_of = {}
    for l in leaves:
        d = 0.0
        for v in paths[l]: d += blen[v]
        depth_of[l] = d
    def patristic(a, b):
        pa, pb = set(paths[a]), paths[b]
        lca = next(v for v in pb if v in pa)
        dl = sum(blen[v] for v in paths[a][:paths[a].index(lca)])
        dr = sum(blen[v] for v in pb[:pb.index(lca)])
        return dl + dr
    # evolve each locus down the tree
    seqs = {l: np.empty(0, np.int8) for l in leaves}
    pieces = {l: [] for l in leaves}
    for li in range(nloc):
        Ll = int(lens[li])
        srate = rng.gamma(0.5, 2.0, Ll)
        node_seq = {0: rng.integers(0, A, Ll).astype(np.int8)}
        order = sorted(parent, key=lambda v: 0 if parent[v] is None else 1)
        for v in sorted(parent):
            if parent[v] is None: continue
            base = node_seq[parent[v]]
            d = blen[v] * locrate[li] * float(np.exp(rng.normal(0, 0.4)))
            p = 1 - np.exp(-d * srate)
            s = base.copy(); mut = rng.random(Ll) < p
            s[mut] = (s[mut] + rng.integers(1, A, int(mut.sum()))) % A
            node_seq[v] = s
        for l in leaves: pieces[l].append(node_seq[l])
    seqs = {l: np.concatenate(pieces[l]) for l in leaves}
    names = [f"t{j:03d}" for j in range(len(leaves))]
    lmap = dict(zip(leaves, names))
    bounds = np.r_[0, np.cumsum(lens)]
    # missingness: contiguous fragment dropped per (leaf, locus)
    present = {}
    chars = {}
    for l in leaves:
        arr = seqs[l].copy()
        for li in range(nloc):
            if rng.random() < pmiss:
                a, b = bounds[li], bounds[li + 1]
                frac = rng.uniform(0.3, 1.0); w = int((b - a) * frac)
                st = a + int(rng.integers(0, (b - a) - w + 1)) if w < b - a else a
                arr[st:st + w] = -1
        chars[lmap[l]] = arr
    # plant events on disjoint taxa
    pool = list(names); rng.shuffle(pool); poolset = set(pool)
    events = []
    def take():
        if not pool: return None
        x = pool.pop(); poolset.discard(x); return x
    leafof = {v: k for k, v in lmap.items()}
    for _ in range(cfg["contam"]):
        t = take(); li = int(rng.integers(nloc))
        if t is None: break
        a, b = bounds[li], bounds[li + 1]
        seg = chars[t][a:b]; pres = seg >= 0
        seg[pres] = rng.integers(0, A, int(pres.sum()))
        events.append({"kind": "contam", "taxon": t, "locus": f"locus{li + 1}"})
    for _ in range(cfg["swap"]):
        t1, t2 = take(), take()
        if t2 is None: break
        li = int(rng.integers(nloc)); a, b = bounds[li], bounds[li + 1]
        tmp = chars[t1][a:b].copy(); chars[t1][a:b] = chars[t2][a:b]; chars[t2][a:b] = tmp
        events.append({"kind": "swap", "taxon": t1, "taxon2": t2, "locus": f"locus{li + 1}",
                       "pat": patristic(leafof[t1], leafof[t2])})
    for _ in range(cfg["nearcopy"]):
        t1 = take()
        if t1 is None: break
        near = sorted((patristic(leafof[t1], leafof[x]), x) for x in names if x != t1 and x in poolset)[:5]
        if not near: break
        t2 = near[int(rng.integers(len(near)))][1]; pool.remove(t2); poolset.discard(t2)
        li = int(rng.integers(nloc)); a, b = bounds[li], bounds[li + 1]
        chars[t2][a:b] = chars[t1][a:b].copy()                      # t1's copy INTO t2
        events.append({"kind": "nearcopy", "donor": t1, "recipient": t2,
                       "locus": f"locus{li + 1}", "pat": near[0][0]})
    for _ in range(cfg["chimera"]):
        t = take()
        if t is None: break
        far = max((patristic(leafof[t], leafof[x]), x) for x in names if x != t)[1]
        lis = rng.choice(nloc, 5, replace=False)
        for li in lis:
            a, b = bounds[li], bounds[li + 1]
            chars[t][a:b] = chars[far][a:b].copy()
        events.append({"kind": "chimera", "taxon": t, "donor": far,
                       "loci": [f"locus{li + 1}" for li in lis]})
    for _ in range(cfg["misalign"]):
        t = take()
        if t is None: break
        li = int(rng.integers(nloc)); a, b = bounds[li], bounds[li + 1]
        w = (b - a) * 3 // 10; st = a + int(rng.integers(0, b - a - w))
        seg = chars[t][st:st + w]
        chars[t][st:st + w] = seg[rng.permutation(len(seg))]
        events.append({"kind": "misalign", "taxon": t, "locus": f"locus{li + 1}"})
    return names, chars, bounds, events


def run_config(cfg):
    rng = np.random.default_rng(cfg["seed"] * 1000 + cfg["id"])
    names, chars, bounds, events = sim_dataset(cfg, rng)
    tag = cfg["tag"]
    d = f"{BASE}/{tag}"; os.makedirs(d, exist_ok=True)
    AB = PROTEIN if cfg["alpha_size"] == 20 else "ACGT"
    with open(f"{d}/a.fasta", "w") as fh:
        for t in names:
            s = "".join(AB[c] if c >= 0 else "-" for c in chars[t])
            fh.write(f">{t}\n{s}\n")
    with open(f"{d}/a.parts.txt", "w") as fh:
        for i in range(len(bounds) - 1):
            fh.write(f"DNA, locus{i + 1} = {bounds[i] + 1}-{bounds[i + 1]}\n")
    r = subprocess.run([sys.executable, TOOL, "--aln", f"{d}/a.fasta", "--partitions",
                        f"{d}/a.parts.txt", "--out", f"{d}/scan",
                        "--datatype", "aa" if cfg["alpha_size"] == 20 else "dna"],
                       capture_output=True, text=True, timeout=cfg.get("timeout", 1800))
    if r.returncode != 0:
        return {"tag": tag, "error": r.stderr.strip().splitlines()[-1] if r.stderr else "?"}
    div, sim = set(), {}
    for ln in open(f"{d}/scan/masked.csv").read().splitlines()[1:]:
        f = ln.split(",")
        if f[3] == "divergent": div.add((f[0], f[1]))
    for ln in open(f"{d}/scan/similar.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        sim[(f[0], frozenset((f[1], f[2])))] = f
    planted = set()
    res = dict(tag=tag, error=None, contam=[0, 0], swap=[0, 0], nearcopy=[0, 0],
               chimera=[0, 0], misalign=[0, 0], dir_ok=0, dir_called=0, fp_div=0, fp_sim=0)
    for ev in events:
        k = ev["kind"]
        if k == "chimera":
            cells = [(l, ev["taxon"]) for l in ev["loci"]]
            res[k][1] += len(cells); res[k][0] += sum(c in div for c in cells)
            planted |= set(cells) | {(l, ev["donor"]) for l in ev["loci"]}
        elif k == "swap":
            cells = [(ev["locus"], ev["taxon"]), (ev["locus"], ev["taxon2"])]
            res[k][1] += 2; res[k][0] += sum(c in div for c in cells)
            planted |= set(cells)
        elif k == "nearcopy":
            key = (ev["locus"], frozenset((ev["donor"], ev["recipient"])))
            res[k][1] += 1
            hit = key in sim
            res[k][0] += hit
            planted |= {(ev["locus"], ev["donor"]), (ev["locus"], ev["recipient"])}
            if hit:
                f = sim[key]
                if f[10] != "NA":
                    res["dir_called"] += 1
                    res["dir_ok"] += f[10] == ev["recipient"]
        else:
            c = (ev["locus"], ev["taxon"])
            res[k][1] += 1; res[k][0] += c in div
            planted.add(c)
    res["fp_div"] = len(div - planted)
    res["fp_sim"] = sum(1 for (l, pr) in sim
                        if not any((l, t) in planted for t in pr))
    res["cells"] = sum(1 for _ in open(f"{d}/scan/cells.tsv")) - 1
    return res


# the large datasets (10,000 taxa x 100 loci), run one at a time with all threads: mid depth like the
# 80-taxon "mid" sets, very deep (saturation between distant taxa), and very shallow (dense sampling,
# nearest neighbours a fraction of a percent apart, where parsimony would do)
_huge = dict(seed=0, n=10000, alpha_size=4, pmiss=0.15, nloc=100, contam=200, swap=200, nearcopy=100, chimera=20, misalign=100, timeout=6 * 3600)
HUGE = {"mid": dict(_huge, id=99, tag="dna_huge_s0", depth=0.2),
        "deep": dict(_huge, id=98, tag="dna_huge_deep_s0", depth=1.5),
        "shallow": dict(_huge, id=97, tag="dna_huge_shallow_s0", depth=0.01)}


def main():
    ncpu = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    os.makedirs(BASE, exist_ok=True)
    cfgs = []
    i = 0
    if len(sys.argv) > 2 and sys.argv[2] == "huge":                 # one large dataset, run alone with all threads
        os.environ["OMP_NUM_THREADS"] = str(ncpu)
        r = run_config(HUGE[sys.argv[3] if len(sys.argv) > 3 else "mid"])
        if r.get("error"): print(f"{r['tag']}: ERROR {r['error']}"); return
        print(f"{r['tag']:16s} contam {r['contam'][0]}/{r['contam'][1]}  swap {r['swap'][0]}/{r['swap'][1]}"
              f"  nearcopy {r['nearcopy'][0]}/{r['nearcopy'][1]} (dir {r['dir_ok']}/{r['dir_called']})"
              f"  chimera {r['chimera'][0]}/{r['chimera'][1]}  misalign {r['misalign'][0]}/{r['misalign'][1]}"
              f"  FPdiv {r['fp_div']}  FPsim {r['fp_sim']}  cells {r['cells']:,}")
        return
    for seed in range(4):
        for depth, dn in ((0.02, "shallow"), (0.2, "mid"), (0.8, "deep")):
            cfgs.append(dict(id=i, seed=seed, tag=f"dna_{dn}_s{seed}", n=80, depth=depth,
                             alpha_size=4, pmiss=0.15, contam=3, swap=4, nearcopy=3,
                             chimera=1, misalign=2)); i += 1
    for seed in range(2):
        cfgs.append(dict(id=i, seed=seed, tag=f"prot_mid_s{seed}", n=60, depth=0.5,
                         alpha_size=20, pmiss=0.15, contam=3, swap=4, nearcopy=3,
                         chimera=1, misalign=2)); i += 1
        cfgs.append(dict(id=i, seed=seed, tag=f"dna_patchy_s{seed}", n=80, depth=0.2,
                         alpha_size=4, pmiss=0.55, contam=3, swap=4, nearcopy=3,
                         chimera=1, misalign=2)); i += 1
        cfgs.append(dict(id=i, seed=seed, tag=f"dna_null_s{seed}", n=80, depth=0.2,
                         alpha_size=4, pmiss=0.15, contam=0, swap=0, nearcopy=0,
                         chimera=0, misalign=0)); i += 1
        cfgs.append(dict(id=i, seed=seed, tag=f"dna_small_s{seed}", n=14, depth=0.1,
                         alpha_size=4, pmiss=0.1, contam=2, swap=2, nearcopy=1,
                         chimera=0, misalign=1)); i += 1
    os.environ["OMP_NUM_THREADS"] = "2"
    with ProcessPoolExecutor(max_workers=min(ncpu // 2, 16)) as ex:
        results = list(ex.map(run_config, cfgs))
    tot = dict(contam=[0, 0], swap=[0, 0], nearcopy=[0, 0], chimera=[0, 0], misalign=[0, 0])
    fp_div = fp_sim = cells = dir_ok = dir_called = 0
    for r in results:
        if r.get("error"):
            print(f"{r['tag']}: ERROR {r['error']}"); continue
        for k in tot:
            tot[k][0] += r[k][0]; tot[k][1] += r[k][1]
        fp_div += r["fp_div"]; fp_sim += r["fp_sim"]; cells += r["cells"]
        dir_ok += r["dir_ok"]; dir_called += r["dir_called"]
        print(f"{r['tag']:16s} contam {r['contam'][0]}/{r['contam'][1]}  swap {r['swap'][0]}/{r['swap'][1]}"
              f"  nearcopy {r['nearcopy'][0]}/{r['nearcopy'][1]} (dir {r['dir_ok']}/{r['dir_called']})"
              f"  chimera {r['chimera'][0]}/{r['chimera'][1]}  misalign {r['misalign'][0]}/{r['misalign'][1]}"
              f"  FPdiv {r['fp_div']}  FPsim {r['fp_sim']}")
    print(f"\nTOTALS over {cells:,} cells: " + "  ".join(f"{k} {v[0]}/{v[1]}" for k, v in tot.items())
          + f"  direction {dir_ok}/{dir_called}  FPdiv {fp_div}  FPsim {fp_sim}")


if __name__ == "__main__":
    main()
