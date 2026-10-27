#!/usr/bin/env python3
"""One focal taxon, every locus: its cell summarised over its comparators, outliers highlighted.

    python plot_taxon.py DIR TAXON [--out FILE] [--label N]

Reads cells.tsv and masked.csv. One point per locus: x = the median over the taxon's comparators
of the divergence the locus expects for that pair (its rate at the pair's m), y = the median over
the same comparators of the divergence observed; these medians are what the cell decision uses.
Diagonal = expectation. Colour is the cell's verdict:
greys pass (darker the more comparators it fails against), reds are masked too divergent (by how
many of the comparators it fails: a majority up to all), blues are masked too similar (hollow blue:
in a similar-side event but not masked), hollow red: ambiguous (fewer comparators than requested).
Right: the same ratios as a stacked histogram on a log count axis. --label N names the N most
extreme masked loci (the part of the name that varies across loci)."""
import argparse, csv, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("taxon"); ap.add_argument("--out"); ap.add_argument("--label", type=int, default=0)
a = ap.parse_args()
R = [r for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t") if r["taxon"] == a.taxon]
if not R: raise SystemExit(f"{a.taxon} has no tested cells in {a.dir}")
side = {}
if os.path.exists(f"{a.dir}/masked.csv"):
    for r in csv.DictReader(open(f"{a.dir}/masked.csv")):
        if r["taxon"] == a.taxon: side[r["locus"]] = r["side"]
per = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/per_locus.tsv"), delimiter="\t")}
FAILS = sorted(int(k[4:]) for k in next(iter(per.values())) if re.match(r"fail\d+$", k)) or [3, 4, 5]
NC, LOW = FAILS[-1], FAILS[0]
loci = [r["locus"] for r in R]
obs = np.array([100.0 - float(r["med_identity"]) for r in R]); exp = np.array([100.0 - float(r["med_expected"]) for r in R])
o = np.array([int(r["overlap"]) for r in R]); nf = np.array([int(r["n_fail"]) for r in R]); nsim = np.array([int(r["n_similar"]) for r in R])
flag = np.array([r["flag"] == "1" for r in R]); amb = np.array([r["ambiguous"] == "1" for r in R])
smk = np.array([side.get(l) in ("similar", "both") for l in loci]); simev = (nsim > 0) & ~smk
mclose = np.array([float(r["pair_div"]) for r in R])
eps = 0.5 / o * 100
y = np.maximum(obs, eps); x = np.maximum(exp, eps); ratio = np.log2(y / x)
inf = (obs * o / 100 >= 1) | (exp * o / 100 >= 1)

INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
MODEL, FLOOR, AMB, SIMRING = "#3b3a37", "#a5a39d", "#e34948", "#2a78d6"
def _ramp(c0, c1, n):
    u, v = np.array(matplotlib.colors.to_rgb(c0)), np.array(matplotlib.colors.to_rgb(c1))
    return [matplotlib.colors.to_hex(u + (v - u) * (i / max(n - 1, 1))) for i in range(n)]
GREYS = dict(zip(range(LOW), _ramp("#a8a6a1", "#4a4944", LOW))); REDS = dict(zip(FAILS, _ramp("#f29492", "#9c1c1c", len(FAILS))))
BLUES = {1: "#86b6ef", 2: "#3987e5", 3: "#184f95"}
plain = ~flag & ~amb & ~smk
fig, (ax, hx) = plt.subplots(1, 2, figsize=(10.5, 4.6), dpi=150, gridspec_kw=dict(width_ratios=(3, 2)))
for k, c in GREYS.items():
    s = plain & (nf == k) & inf
    if s.any(): ax.scatter(x[s], y[s], s=12, color=c, alpha=0.85, linewidths=0, zorder=2)
for k, c in REDS.items():
    s = flag & (nf == k) & inf
    if s.any(): ax.scatter(x[s], y[s], s=20, color=c, edgecolors=SURF, linewidths=0.5, zorder=4)
s = amb & inf
if s.any(): ax.scatter(x[s], y[s], s=26, facecolors="none", edgecolors=AMB, linewidths=0.8, zorder=5)
s = simev & inf
if s.any(): ax.scatter(x[s], y[s], s=30, facecolors="none", edgecolors=SIMRING, linewidths=0.8, zorder=5)
for k, c in BLUES.items():
    s = smk & (np.minimum(nsim, 3) == k) & inf
    if s.any(): ax.scatter(x[s], y[s], s=20, color=c, edgecolors=SURF, linewidths=0.5, zorder=6)
lim = (min(x[inf].min(), y[inf].min()) / 1.5, max(x[inf].max(), y[inf].max()) * 1.5); g = np.array(lim)
ax.plot(g, g, color=MODEL, lw=1.3, zorder=1)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
pre = os.path.commonprefix(loci); suf = os.path.commonprefix([l[::-1] for l in loci])[::-1]
short = lambda l: l[len(pre):len(l) - len(suf)] or l
for k, i in enumerate([i for i in np.argsort(-np.abs(ratio)) if (flag[i] or smk[i]) and inf[i]][:a.label]):
    ax.annotate(short(loci[i]), (x[i], y[i]), xytext=(6, -10 - 9 * k), textcoords="offset points", fontsize=6, color=INK, arrowprops=dict(arrowstyle="-", color=FLOOR, lw=0.5))
ax.set_xlabel("expected locus divergence, median over comparators (%)", fontsize=8, color=MUTED); ax.set_ylabel("observed locus divergence, median over comparators (%)", fontsize=8, color=MUTED)
ax.set_title(f"{a.taxon}  ·  {os.path.basename(os.path.normpath(a.dir))}", fontsize=10, color=INK, loc="left", fontweight="bold", pad=16)
ax.text(0, 1.012, f"{len(loci)} loci ({int((~inf).sum())} with no difference expected or seen) · closest comparator at m = {np.median(mclose):.2f}% (median over loci) · "
        f"{int(flag.sum())} masked too divergent, {int(smk.sum())} too similar" + (f", {int(amb.sum())} ambiguous" if amb.any() else ""), transform=ax.transAxes, fontsize=7, color=MUTED, va="bottom")
bins = np.linspace(min(ratio[inf].min(), -2.5), max(ratio[inf].max(), 2.5), 41)
# the nine fills of the scatter, stacked from the base: reds, blues, greys; within each colour from fewest pairs to most
groups = [(flag & (nf == k) & inf, REDS[k]) for k in sorted(FAILS)] + \
         [(smk & (np.minimum(nsim, 3) == k) & inf, BLUES[k]) for k in (1, 2, 3)] + \
         [(plain & (nf == k) & inf, GREYS[k]) for k in sorted(GREYS)]
hx.hist([ratio[s] for s, _ in groups], bins=bins, color=[c for _, c in groups], stacked=True, edgecolor=SURF, linewidth=0.4, zorder=2)
hx.set_yscale("log"); hx.set_ylim(0.8, None)
hx.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0))); hx.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
hx.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(round(v)):,}" if v >= 1 else ""))
hx.axvline(0, color=MODEL, lw=1.3)
hx.set_xlabel("log2 (observed / expected), medians over comparators", fontsize=8, color=MUTED); hx.set_ylabel("loci", fontsize=8, color=MUTED)
hx.set_title(f"{int(inf.sum())} loci with a difference expected or seen · log count axis", fontsize=7, color=MUTED, loc="left", pad=4)
handles = [Line2D([], [], marker="o", ls="", color=c, ms=4.5, label=("passes · " if k == 0 else "") + f"fails {k} of {NC}") for k, c in GREYS.items()] + \
          [Line2D([], [], marker="o", ls="", color=c, mec=SURF, ms=5.5, label=("masked, too divergent · " if k == LOW else "") + f"fails {k} of {NC}") for k, c in REDS.items()] + \
          [Line2D([], [], marker="o", ls="", color=BLUES[k], mec=SURF, ms=5.5, label=("masked, too similar · " if k == 1 else "") + f"{k}{'+' if k == 3 else ''} pair{'s' if k > 1 else ''}") for k in (1, 2, 3)] + \
          [Line2D([], [], marker="o", ls="", mfc="none", mec=SIMRING, ms=6, label="in a similar event, not masked"),
           Line2D([], [], marker="o", ls="", mfc="none", mec=AMB, ms=6, label="ambiguous"), Line2D([], [], color=MODEL, lw=1.3, label="expected")]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, -0.04), labelcolor=INK, columnspacing=1.2, handletextpad=0.4)
for axx in (ax, hx):
    for s_ in ("top", "right"): axx.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"): axx.spines[s_].set_color(SPINE)
    axx.tick_params(colors=MUTED, labelsize=7, length=0); axx.grid(True, color=GRID, lw=0.5); axx.set_axisbelow(True)
fig.tight_layout(rect=(0, 0.07, 1, 1))
out = a.out or f"{a.dir}/taxon_{a.taxon}.png"
fig.savefig(out, bbox_inches="tight"); print("wrote", out)
