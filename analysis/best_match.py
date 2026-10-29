#!/usr/bin/env python3
"""For one taxon: which other taxon each of its sequences matches most closely, masked cells vs the rest.
    python best_match.py DIR ALIGNMENT TAXON [--out FILE] [--threads 8]
DIR is a LocusMasker run (0.11 or later: masked.csv and report.npz); ALIGNMENT the NEXUS file it was run on. A run writes no
divergence matrix, so the genome-wide divergences m are recomputed from ALIGNMENT with the run's options
(locus_outliers.divergence_matrix: the same m the run used).
"""
import argparse, csv, hashlib, io, json, os, re, sys, zipfile, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); import locus_outliers as lo
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("aln"); ap.add_argument("taxon"); ap.add_argument("--out"); ap.add_argument("--min-overlap", type=int, default=30)
ap.add_argument("--threads", type=int, default=8, help="threads for recomputing m")
a = ap.parse_args()

def read_nexus_simple(path):
    names, seqs, charsets, inmat = [], [], {}, False
    for line in open(path):
        s = line.strip()
        if not inmat:
            if s.lower() == "matrix": inmat = True; continue
            m = re.match(r"charset\s+(\S+)\s*=\s*(.+);", s, re.I)
            if m:
                sites = []
                for part in m.group(2).split():
                    r = re.match(r"^(\d+)-(\d+)$", part)
                    if r: sites.extend(range(int(r.group(1)) - 1, int(r.group(2))))
                    elif part.isdigit(): sites.append(int(part) - 1)
                charsets[m.group(1)] = np.array(sites, dtype=int)
            continue
        if s == ";": inmat = False; continue
        if not s: continue
        p = s.split(None, 1); names.append(p[0]); seqs.append(p[1].replace(" ", "") if len(p) > 1 else "")
    return names, seqs, charsets

names, seqs, charsets = read_nexus_simple(a.aln)
L = len(seqs[0]); assert all(len(s) == L for s in seqs), "interleaved nexus: not handled by this prototype"
X = np.frombuffer("".join(s.upper() for s in seqs).encode(), dtype=np.uint8).reshape(len(names), L)
missing = np.isin(X, np.frombuffer(b"-?XN*", dtype=np.uint8))
ti = names.index(a.taxon)
z = zipfile.ZipFile(f"{a.dir}/report.npz"); npy = lambda k: np.load(io.BytesIO(z.read(k))); man = json.loads(z.read("manifest.json")); opt = man["options"]
with open(a.aln, "rb") as fh:
    if hashlib.sha256(fh.read()).hexdigest() != man["input"]["sha256"]: print(f"warning: {a.aln} is not the file the run read ({man['input']['file']}); m is recomputed from it", file=sys.stderr)
mt, m, _ = lo.divergence_matrix(a.aln, None, opt["datatype"], fit_overlap=opt["fit_overlap"], n_comp=opt["n_comp"], threads=a.threads)
if list(mt) != json.loads(z.read("taxa.json"))["names"]: sys.exit(f"{a.aln} does not hold the taxa of the run in {a.dir}")
mi = {t: i for i, t in enumerate(mt)}; lnames = json.loads(z.read("loci.json"))["names"]
# genome-wide divergence of every taxon from the focal one, in alignment order (NaN: excluded taxon)
mrow = np.array([m[mi[a.taxon], mi[t]] if t in mi and t != a.taxon else np.nan for t in names])
# every taxon in order of m from the focal one: the run's own 5 closest first (so its tie-break), then the rest
head = [names.index(mt[j]) for j in npy("taxa/closest.npy")[mi[a.taxon]] if j >= 0]
order = np.array(head + [i for i in np.argsort(np.where(np.isnan(mrow), np.inf, mrow), kind="stable").tolist() if i not in head])
rank = np.full(len(names), 10 ** 6); rank[order] = np.arange(len(names))
nn = [names[i] for i in order[:5]]
masked = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/masked.csv")) if r["taxon"] == a.taxon}
tested = [lnames[li] for li in np.flatnonzero(npy("cells/n_comp.npy")[mi[a.taxon]] > 0)]   # the loci where the run tested it (it had comparators there)
rows = []
for loc in tested:
    sites = charsets[loc]; Y = X[:, sites]; ok = ~missing[:, sites]
    both = ok & ok[ti][None, :]; ov = both.sum(1); nf = int(ok[ti].sum())
    mism = nf - ((Y == Y[ti][None, :]) & both).sum(1)                    # over the focal sequence's sites: a site the other taxon lacks counts as a mismatch
    mism = mism.astype(float); mism[ti] = np.nan; mism[ov < a.min_overlap] = np.nan; ov = np.full_like(ov, nf)
    if np.all(np.isnan(mism)): continue
    best = np.nanmin(mism); tied = np.flatnonzero(mism == best)
    j = int(tied[np.argmin(rank[tied])])                                   # ties go to the closest relative genome-wide
    j0 = int(order[0]) if not np.isnan(mism[order[0]]) else j             # the taxon's closest relative genome-wide
    r = masked.get(loc)
    rows.append(dict(locus=loc, best=names[j], ident=100 * (1 - best / ov[j]), overlap=int(ov[j]), ntied=len(tied),
                     margin=int(mism[j0] - best) if not np.isnan(mism[j0]) else -1, ident_nn=100 * (1 - mism[j0] / ov[j0]) if not np.isnan(mism[j0]) else np.nan,
                     best_rank=int(rank[j]), masked=r is not None,
                     reason=r["reason"] if r else "", divergent_fails=int(r["divergent_fails"] or 0) if r else 0,
                     similar_partners=int(r["similar_partners"] or 0) if r else 0, similar_partner=r["similar_partner"] if r else ""))
print("closest relatives genome-wide:", nn)
with open(f"{a.dir}/bestmatch_{a.taxon}.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
import collections
cm = collections.Counter(r["best"] for r in rows if r["masked"]); cu = collections.Counter(r["best"] for r in rows if not r["masked"])
print(f"{a.taxon}: {sum(r['masked'] for r in rows)} masked of {len(rows)} tested cells")
print("best match of MASKED cells:"); [print(f"  {t:14s} {n:4d}") for t, n in cm.most_common(15)]
print("best match of UNMASKED cells:"); [print(f"  {t:14s} {n:4d}") for t, n in cu.most_common(15)]
if not cm: sys.exit(f"{a.taxon}: no masked sequences, so no comparison or figure (wrote {a.dir}/bestmatch_{a.taxon}.tsv)")
im = np.array([r["ident"] for r in rows if r["masked"]]); iu = np.array([r["ident"] for r in rows if not r["masked"]])
mm = np.array([r["margin"] for r in rows if r["masked"]]); mu = np.array([r["margin"] for r in rows if not r["masked"]])
print(f"margin (mismatches fewer than to the closest relative): masked median {np.median(mm):.0f}, share > 0: {np.mean(mm > 0):.2f}, > 5: {np.mean(mm > 5):.2f}; unmasked share > 0: {np.mean(mu > 0):.2f}, > 5: {np.mean(mu > 5):.2f}")
rm = np.array([r["best_rank"] for r in rows if r["masked"]]); ru = np.array([r["best_rank"] for r in rows if not r["masked"]])
print(f"best match among the 5 closest relatives: masked {np.mean(rm < 5):.2f}, unmasked {np.mean(ru < 5):.2f}")
print(f"identity to best match: masked median {np.median(im):.1f}% (10th {np.percentile(im,10):.1f}, 90th {np.percentile(im,90):.1f}); unmasked median {np.median(iu):.1f}% (10th {np.percentile(iu,10):.1f}, 90th {np.percentile(iu,90):.1f})")

# ------------------------------------------------------------------ figure
import os
from matplotlib.patches import Patch
INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
RED3 = ["#d9433f", "#a8201f", "#6f0d0f"]; BLUE3 = ["#86b6ef", "#3987e5", "#184f95"]
NC = 5; FAILS = (3, 4, 5); REDS = dict(zip(FAILS, RED3)); BLUES = {1: BLUE3[0], 2: BLUE3[1], 3: BLUE3[2]}
cats = [(f"masked, too divergent · fails {j} of {NC}", REDS[j], lambda r, j=j: r["reason"] != "similar" and min(max(r["divergent_fails"], 3), NC) == j) for j in FAILS] + \
       [(f"masked, copied · {j}{'+' if j == 3 else ''} too-similar partner{'s' if j > 1 else ''}", BLUES[j], lambda r, j=j: r["reason"] == "similar" and min(r["similar_partners"], 3) == j) for j in (1, 2, 3)]
nm = sum(r["masked"] for r in rows)
taxa_sorted = sorted(cm, key=lambda t: (-cm[t], t))
comp = np.array([[sum(1 for r in rows if r["masked"] and r["best"] == t and f(r)) for _, _, f in cats] for t in taxa_sorted], float)
obs = comp.sum(1); close = np.array([t in nn for t in taxa_sorted])
W = max(10.0, 0.22 * len(taxa_sorted) + 1.5)
fig, ax = plt.subplots(figsize=(W, 4.8), dpi=150); fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
fig.subplots_adjust(left=0.6 / W, right=1 - 0.15 / W, top=0.71, bottom=0.25)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
ax.spines["left"].set_color(SPINE); ax.spines["bottom"].set_color(SPINE); ax.tick_params(colors=MUTED, labelsize=7, length=0)
ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
xs = np.arange(len(taxa_sorted)); bottom = np.zeros(len(xs))
for k, (lab, c, _) in enumerate(cats):
    h = comp[:, k]
    if h.any(): ax.bar(xs, h, bottom=bottom, width=0.82, color=c, edgecolor=SURF, lw=0.4)
    bottom += h
for x_, o_ in zip(xs, obs): ax.text(x_, o_ + obs.max() * 0.012, f"{int(o_)}", ha="center", va="bottom", fontsize=6.2, color=INK)
ax.set_xticks(xs); ax.set_xticklabels(taxa_sorted, rotation=90, fontsize=6.4, color=INK)
for lab, c_ in zip(ax.get_xticklabels(), close):
    if c_: lab.set_fontweight("bold")
ax.set_xlim(-0.7, len(xs) - 0.3); ax.set_ylabel("masked sequences", color=MUTED, fontsize=8); ax.set_xlabel("the taxon the masked sequence matches best  (identity over the masked sequence's sites; a site the other taxon lacks counts as a mismatch)", color=MUTED, fontsize=8)
handles = [Patch(facecolor=c, edgecolor=SURF, label=lab) for lab, c, _ in cats if comp[:, [l for l, _ in enumerate(cats)][cats.index((lab, c, _))]].any()]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.6 / W, 0.83), ncol=3, frameon=False, fontsize=6.5, labelcolor=INK, columnspacing=2.0, handletextpad=0.5)
ds = os.path.basename(os.path.normpath(a.dir))
fig.text(0.6 / W, 0.97, f"{ds}  ·  {a.taxon}  ·  best match of each masked sequence", fontsize=11, fontweight="bold", color=INK, ha="left", va="top")
fig.text(0.6 / W, 0.915, f"{nm} of {len(rows):,} sequences masked ({100 * nm / len(rows):.1f}%)  ·  median identity to the best match: unmasked {np.median(iu):.1f}%, masked {np.median(im):.1f}%", fontsize=7.4, color=INK, ha="left", va="top")
fig.text(0.6 / W, 0.885, f"bold: the 5 closest relatives genome-wide, which are the best match of {100 * np.mean(ru < 5):.0f}% of unmasked sequences and {100 * np.mean(rm < 5):.0f}% of masked ones", fontsize=7.4, color=INK, ha="left", va="top")
out = a.out or f"{a.dir}/bestmatch_{a.taxon}.png"; fig.savefig(out, facecolor=SURF); print("wrote", out)
