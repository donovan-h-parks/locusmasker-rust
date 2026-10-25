#!/usr/bin/env python3
"""One taxon pair, every locus: observed against expected divergence, outliers highlighted.

    python plot_pair.py DIR TAXON_A TAXON_B [--out FILE] [--label N]

Reads pairs.npy (so the pair must have been compared: one must be among the other's comparators)
and masked.csv. Left: one point per locus, x = the divergence the locus expects for this pair
(its rate at the pair's m), y = what it shows; the diagonal is the expectation, the dashed lines
the effect floors (2x above, 1/3 below). Right: the same ratios as a histogram. Colours are the
pair's own verdict at each locus: red = too divergent, blue = too similar, grey = within
expectation. A dark ring marks loci where TAXON_A's cell was masked (majority over its
comparators, or an attributed similar event). The --label N most extreme loci are named."""
import argparse, csv, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--out"); ap.add_argument("--label", type=int, default=6)
a = ap.parse_args()
if not os.path.exists(f"{a.dir}/pairs.npy"): raise SystemExit(f"{a.dir}/pairs.npy not found: it is written by every run of locus_outliers.py")
P = np.load(f"{a.dir}/pairs.npy"); taxa = open(f"{a.dir}/taxa.txt").read().split()
lnames = [r.split("\t", 1)[0] for r in open(f"{a.dir}/per_locus.tsv").read().splitlines()[1:]]
ia, ib = taxa.index(a.a), taxa.index(a.b)
rows = {}
for q in P[((P["taxon"] == ia) & (P["comparator"] == ib)) | ((P["taxon"] == ib) & (P["comparator"] == ia))]:
    l = lnames[q["locus"]]
    if q["taxon"] == ia or l not in rows:
        rows[l] = dict(locus=l, taxon=taxa[q["taxon"]], comparator=taxa[q["comparator"]], pair_div=str(q["pair_div"]), overlap=str(int(q["overlap"])),
                       identity=str(q["identity"]), expected=str(q["expected"]), excess=str(q["excess"]), fail=str(int(q["fail"])), fail_similar=str(int(q["fail_similar"])))
if not rows: raise SystemExit(f"{a.a} and {a.b} were never compared (neither is among the other's comparators)")
masked = {}
if os.path.exists(f"{a.dir}/masked.csv"):
    for r in csv.DictReader(open(f"{a.dir}/masked.csv")):
        if r["taxon"] == a.a: masked[r["locus"]] = r["side"]
loci = sorted(rows); R = [rows[l] for l in loci]
obs = np.array([100.0 - float(r["identity"]) for r in R]); exp = np.array([100.0 - float(r["expected"]) for r in R])
o = np.array([int(r["overlap"]) for r in R]); m = float(R[0]["pair_div"])
hi = np.array([r["fail"] == "1" for r in R]); lo = np.array([r["fail_similar"] == "1" for r in R])
ring = np.array([l in masked for l in loci])
eps = 0.5 / o * 100                                          # half a mismatch: keeps zero divergences on a log axis
y = np.maximum(obs, eps); x = np.maximum(exp, eps); ratio = np.log2(y / x)
inf = (obs * o / 100 >= 1) | (exp * o / 100 >= 1)              # loci where at least one side expects or shows a difference

INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
GREY, RED, BLUE, MODEL, FLOOR = "#9a9893", "#e34948", "#3987e5", "#3b3a37", "#a5a39d"
col = np.where(hi, RED, np.where(lo, BLUE, GREY))
fig, (ax, hx) = plt.subplots(1, 2, figsize=(10.5, 4.6), dpi=150, gridspec_kw=dict(width_ratios=(3, 2)))
order = [i for i in np.argsort(hi | lo) if inf[i]]              # outliers drawn last; uninformative loci not drawn
ax.scatter(x[order], y[order], s=14, c=col[order], alpha=0.85, linewidths=0, zorder=3)
if (ring & inf).any(): ax.scatter(x[ring & inf], y[ring & inf], s=40, facecolors="none", edgecolors=INK, linewidths=0.8, zorder=4)
lim = (min(x[inf].min(), y[inf].min()) / 1.5, max(x[inf].max(), y[inf].max()) * 1.5)
g = np.array(lim)
ax.plot(g, g, color=MODEL, lw=1.3, zorder=2); ax.plot(g, 2 * g, color=FLOOR, lw=0.9, ls=(0, (3, 2)), zorder=2); ax.plot(g, g / 3, color=FLOOR, lw=0.9, ls=(0, (3, 2)), zorder=2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
pre = os.path.commonprefix(loci); suf = os.path.commonprefix([l[::-1] for l in loci])[::-1]   # the part of a name that varies
short = lambda l: l[len(pre):len(l) - len(suf)] or l
for k, i in enumerate([i for i in np.argsort(-np.abs(ratio)) if (hi[i] or lo[i]) and inf[i]][:a.label]):
    ax.annotate(short(loci[i]), (x[i], y[i]), xytext=(6, -10 - 9 * k), textcoords="offset points",
                fontsize=6, color=INK, arrowprops=dict(arrowstyle="-", color=FLOOR, lw=0.5))
ax.set_xlabel("expected locus divergence for this pair (%)", fontsize=8, color=MUTED); ax.set_ylabel("observed locus divergence (%)", fontsize=8, color=MUTED)
ax.set_title(f"{a.a} vs {a.b}", fontsize=10, color=INK, loc="left", fontweight="bold", pad=16)
ax.text(0, 1.012, f"{len(loci)} loci ({int((~inf).sum())} with no difference expected or seen) · pair divergence m = {m:.2f}% · "
        f"{hi.sum()} too divergent, {lo.sum()} too similar · {ring.sum()} of {a.a}'s cells masked", transform=ax.transAxes, fontsize=7, color=MUTED, va="bottom")
bins = np.linspace(min(ratio[inf].min(), -2.5), max(ratio[inf].max(), 2.5), 41)
groups = [(lo & inf, BLUE), (hi & inf, RED), (~hi & ~lo & inf, GREY)]          # stacked, outliers at the base and grey on top; log count axis so a few loci still show
hx.hist([ratio[sel] for sel, _ in groups], bins=bins, color=[c for _, c in groups], stacked=True, edgecolor=SURF, linewidth=0.4, zorder=2)
hx.set_yscale("log"); hx.set_ylim(0.8, None)
hx.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0))); hx.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
hx.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(round(v)):,}" if v >= 1 else ""))
hx.set_title(f"{int(inf.sum())} loci with a difference expected or seen · log count axis", fontsize=7, color=MUTED, loc="left", pad=4)
hx.axvline(0, color=MODEL, lw=1.3); hx.axvline(1, color=FLOOR, lw=0.9, ls=(0, (3, 2))); hx.axvline(-np.log2(3), color=FLOOR, lw=0.9, ls=(0, (3, 2)))
hx.set_xlabel("log2 (observed / expected)", fontsize=8, color=MUTED); hx.set_ylabel("loci", fontsize=8, color=MUTED)
handles = [Line2D([], [], marker="o", ls="", color=GREY, ms=5, label="within expectation"),
           Line2D([], [], marker="o", ls="", color=RED, ms=5, label="too divergent (pair fails)"),
           Line2D([], [], marker="o", ls="", color=BLUE, ms=5, label="too similar (pair fails)"),
           Line2D([], [], marker="o", ls="", mfc="none", mec=INK, ms=6, label=f"{a.a}'s cell masked"),
           Line2D([], [], color=MODEL, lw=1.3, label="expected"), Line2D([], [], color=FLOOR, lw=0.9, ls=(0, (3, 2)), label="effect floors (2x, 1/3)")]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.01), labelcolor=INK)
for axx in (ax, hx):
    for s in ("top", "right"): axx.spines[s].set_visible(False)
    for s in ("left", "bottom"): axx.spines[s].set_color(SPINE)
    axx.tick_params(colors=MUTED, labelsize=7, length=0); axx.grid(True, color=GRID, lw=0.5); axx.set_axisbelow(True)
fig.tight_layout(rect=(0, 0.05, 1, 1))
out = a.out or f"{a.dir}/pair_{a.a}_vs_{a.b}.png"
fig.savefig(out, bbox_inches="tight"); print("wrote", out)
