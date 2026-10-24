#!/usr/bin/env python3
"""For one taxon: which other taxon each of its sequences matches most closely, masked cells vs the rest.
    python best_match.py DIR ALIGNMENT TAXON [--out FILE]
"""
import argparse, csv, re, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("aln"); ap.add_argument("taxon"); ap.add_argument("--out"); ap.add_argument("--min-overlap", type=int, default=30)
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
m = np.load(f"{a.dir}/m.npy"); mt = open(f"{a.dir}/taxa.txt").read().split(); mi = {t: i for i, t in enumerate(mt)}
# genome-wide divergence of every taxon from the focal one, in alignment order (NaN: excluded taxon)
mrow = np.array([m[mi[a.taxon], mi[t]] if t in mi and t != a.taxon else np.nan for t in names])
rank = np.full(len(names), 10 ** 6); order = np.argsort(np.where(np.isnan(mrow), np.inf, mrow)); rank[order] = np.arange(len(names))
nn = [names[i] for i in order[:5]]
masked = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/masked.csv")) if r["taxon"] == a.taxon}
tested = [r["locus"] for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t") if r["taxon"] == a.taxon]
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
                     side=r["side"] if r else "", n_fail=int(r["n_fail"]) if r else 0, n_sim=int(r["n_similar"]) if r else 0,
                     partner=r["partner"] if r else ""))
print("closest relatives genome-wide:", nn)
import collections
cm = collections.Counter(r["best"] for r in rows if r["masked"]); cu = collections.Counter(r["best"] for r in rows if not r["masked"])
print(f"{a.taxon}: {sum(r['masked'] for r in rows)} masked of {len(rows)} tested cells")
print("best match of MASKED cells:"); [print(f"  {t:14s} {n:4d}") for t, n in cm.most_common(15)]
print("best match of UNMASKED cells:"); [print(f"  {t:14s} {n:4d}") for t, n in cu.most_common(15)]
im = np.array([r["ident"] for r in rows if r["masked"]]); iu = np.array([r["ident"] for r in rows if not r["masked"]])
mm = np.array([r["margin"] for r in rows if r["masked"]]); mu = np.array([r["margin"] for r in rows if not r["masked"]])
print(f"margin (mismatches fewer than to the closest relative): masked median {np.median(mm):.0f}, share > 0: {np.mean(mm > 0):.2f}, > 5: {np.mean(mm > 5):.2f}; unmasked share > 0: {np.mean(mu > 0):.2f}, > 5: {np.mean(mu > 5):.2f}")
rm = np.array([r["best_rank"] for r in rows if r["masked"]]); ru = np.array([r["best_rank"] for r in rows if not r["masked"]])
print(f"best match among the 5 closest relatives: masked {np.mean(rm < 5):.2f}, unmasked {np.mean(ru < 5):.2f}")
print(f"identity to best match: masked median {np.median(im):.1f}% (10th {np.percentile(im,10):.1f}, 90th {np.percentile(im,90):.1f}); unmasked median {np.median(iu):.1f}% (10th {np.percentile(iu,10):.1f}, 90th {np.percentile(iu,90):.1f})")
with open(f"{a.dir}/bestmatch_{a.taxon}.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ figure
import os
from matplotlib.patches import Patch
INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
RED3 = ["#d9433f", "#a8201f", "#6f0d0f"]; BLUE3 = ["#86b6ef", "#3987e5", "#184f95"]
NC = 5; FAILS = (3, 4, 5); REDS = dict(zip(FAILS, RED3)); BLUES = {1: BLUE3[0], 2: BLUE3[1], 3: BLUE3[2]}
cats = [(f"masked, too divergent · fails {j} of {NC}", REDS[j], lambda r, j=j: r["side"] != "similar" and min(max(r["n_fail"], 3), NC) == j) for j in FAILS] + \
       [(f"masked, too similar · {j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}", BLUES[j], lambda r, j=j: r["side"] == "similar" and min(r["n_sim"], 3) == j) for j in (1, 2, 3)]
nm = sum(r["masked"] for r in rows); nu = len(rows) - nm
taxa_sorted = sorted({t for t in cm} | {t for t in cu if cu[t] / max(nu, 1) * nm >= 1}, key=lambda t: (-cm[t], -cu[t]))
comp = np.array([[sum(1 for r in rows if r["masked"] and r["best"] == t and f(r)) for _, _, f in cats] for t in taxa_sorted], float)
obs = comp.sum(1); expd = np.array([cu[t] / max(nu, 1) * nm for t in taxa_sorted]); close = np.array([t in nn for t in taxa_sorted])
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
ax.hlines(expd, xs - 0.41, xs + 0.41, color=INK, lw=1.1, zorder=4)
for x_, o_ in zip(xs, obs): ax.text(x_, o_ + obs.max() * 0.012, f"{int(o_)}", ha="center", va="bottom", fontsize=6.2, color=INK)
ax.set_xticks(xs); ax.set_xticklabels(taxa_sorted, rotation=90, fontsize=6.4, color=INK)
for lab, c_ in zip(ax.get_xticklabels(), close):
    if c_: lab.set_fontweight("bold")
ax.set_xlim(-0.7, len(xs) - 0.3); ax.set_ylabel("masked sequences", color=MUTED, fontsize=8); ax.set_xlabel("the taxon the masked sequence matches best  (identity over the masked sequence's sites; a site the other taxon lacks counts as a mismatch)", color=MUTED, fontsize=8)
handles = [Patch(facecolor=c, edgecolor=SURF, label=lab) for lab, c, _ in cats if comp[:, [l for l, _ in enumerate(cats)][cats.index((lab, c, _))]].any()] + \
          [plt.Line2D([], [], color=INK, lw=1.1, label="count expected if masked sequences matched like unmasked ones")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.6 / W, 0.83), ncol=3, frameon=False, fontsize=6.5, labelcolor=INK, columnspacing=2.0, handletextpad=0.5)
ds = os.path.basename(os.path.normpath(a.dir))
fig.text(0.6 / W, 0.97, f"{ds}  ·  {a.taxon}  ·  best match of each masked sequence", fontsize=11, fontweight="bold", color=INK, ha="left", va="top")
fig.text(0.6 / W, 0.915, f"{nm} of {len(rows):,} sequences masked ({100 * nm / len(rows):.1f}%)  ·  median identity to the best match: unmasked {np.median(iu):.1f}%, masked {np.median(im):.1f}%", fontsize=7.4, color=INK, ha="left", va="top")
fig.text(0.6 / W, 0.885, f"bold: the 5 closest relatives genome-wide, which are the best match of {100 * np.mean(ru < 5):.0f}% of unmasked sequences and {100 * np.mean(rm < 5):.0f}% of masked ones", fontsize=7.4, color=INK, ha="left", va="top")
out = a.out or f"{a.dir}/bestmatch_{a.taxon}.png"; fig.savefig(out, facecolor=SURF); print("wrote", out)
