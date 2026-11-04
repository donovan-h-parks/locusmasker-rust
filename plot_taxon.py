#!/usr/bin/env python3
"""One focal taxon, every locus: its cell summarised over its comparators, outliers highlighted.

    python plot_taxon.py DIR TAXON [--out FILE] [--label N]

Reads DIR/report.npz (written by every run of locusmasker.py 0.11.0 and later unless --no-report).
One point per locus: x = the median over the taxon's comparators of the divergence the locus expects
for that pair (its rate at the pair's m), y = the median over the same comparators of the divergence
observed; these medians are what the too-divergent test uses. Diagonal = expectation. Colour is the
cell's verdict: greys pass (darker the more comparators it fails against), reds are masked too
divergent (by how many of the comparators it fails: a majority up to all), blues are masked too
similar, the copy in a too-similar pair whose direction was identified (by how many taxa the
sequence is significantly too similar to); hollow blue: in a too-similar pair but not masked (the
donor, or the direction undetermined); hollow red: ambiguous (fewer comparators than requested).
Right: the same ratios as a stacked histogram on a log count axis. --label N names the N most
extreme masked loci (the part of the name that varies across loci). A sparse report (more than 8
million taxa x loci) carries these numbers only for sequences with a failing test: only those are
drawn, and the subtitle says so."""
import argparse, io, json, os, zipfile, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("taxon"); ap.add_argument("--out"); ap.add_argument("--label", type=int, default=0)
a = ap.parse_args()
p = os.path.join(a.dir, "report.npz")
if not os.path.exists(p):
    raise SystemExit(f"{a.dir} was written by LocusMasker 0.10 or earlier (cells.tsv, no report.npz): this script reads the report.npz of 0.11.0 and later; re-run "
                     "locusmasker.py, or use the plot scripts of the version that made the run" if os.path.exists(f"{a.dir}/cells.tsv") else
                     f"{p} not found: this script reads report.npz, which every run of locusmasker.py writes unless --no-report")
Z = zipfile.ZipFile(p); man = json.loads(Z.read("manifest.json"))
if man.get("format_version", 0) < 4: raise SystemExit(f"{p} is report format {man.get('format_version')} (LocusMasker {man.get('version', '?')}); this script reads format 4 and later (0.11.0)")
def npy(k): return np.load(io.BytesIO(Z.read(k)))
TJ = json.loads(Z.read("taxa.json")); LJ = json.loads(Z.read("loci.json"))
if a.taxon not in TJ["names"]: raise SystemExit(f"{a.taxon} is not a taxon of the run in {a.dir}")
t = TJ["names"].index(a.taxon)
NC, LOW = int(man["options"]["K_eff"]), int(man["options"]["M_eff"]); FAILS = list(range(LOW, NC + 1))   # comparators; the majority that fails
L = np.flatnonzero(npy("cells/n_comp.npy")[t] > 0)                  # the taxon's tested loci, in locus order
if not len(L): raise SystemExit(f"{a.taxon} has no tested cells in {a.dir}" + ("" if TJ["status"][t] == "tested" else f" ({TJ['status'][t]})"))
nshort = 0
if man.get("carried", {}).get("numeric_columns") == "sparse":       # numbers only for the cells with a failing test
    st, sl = npy("cells/sparse_taxon.npy"), npy("cells/sparse_locus.npy"); rows = np.flatnonzero(st == t)
    at = dict(zip(sl[rows].tolist(), rows.tolist())); keep = np.array([l in at for l in L.tolist()], bool)
    nshort = int((~keep).sum()); L = L[keep]
    if not len(L): raise SystemExit(f"the report in {a.dir} is sparse (numbers only for sequences with a failing test) and none of {a.taxon}'s has one")
    ri = np.array([at[l] for l in L.tolist()], dtype=np.intp); col = lambda k: npy(f"cells/{k}.npy")[ri]
else: col = lambda k: npy(f"cells/{k}.npy")[t, L]
loci = [LJ["names"][l] for l in L]
side = npy("cells/side.npy")[t, L]                                  # bits: 1 too divergent, 2 the copy in a too-similar pair, 4 ambiguous, 16 too divergent, not masked
obs = 100.0 - col("med_identity").astype(float); exp = 100.0 - col("med_expected").astype(float)
o = col("closest_overlap").astype(float); nf = npy("cells/n_fail.npy")[t, L].astype(int); nsim = npy("cells/n_similar.npy")[t, L].astype(int)
amb = (side & 4) > 0; flag = ((side & 17) > 0) & ~amb
smk = (side & 2) > 0; simev = (nsim > 0) & ~smk                     # nsim: the taxa the sequence is significantly too similar to
# the closest comparator's genome-wide divergence (max(m, m_res), as the too-divergent test uses it): its expected divergence
# at the locus taken back through the locus's expectation curve, P (1 - exp(-rate c(m) / P)), c(p) = -P ln(1 - p / P)
P = np.array(LJ["plateau"], dtype=float)[L]; rate = np.array([np.nan if v is None else float(v) for v in LJ["test_rate"]])[L]
with np.errstate(invalid="ignore", divide="ignore"):
    cexp = -P * np.log(1.0 - np.minimum(100.0 - col("closest_expected").astype(float), 0.98 * P) / P)
    mclose = P * (1.0 - np.exp(-cexp / (rate * P)))
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
for axis in (ax.xaxis, ax.yaxis):                                       # minor ticks labelled only when no power of 10 is in view
    axis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(0, 0.4)))
pre = os.path.commonprefix(loci); suf = os.path.commonprefix([l[::-1] for l in loci])[::-1]
short = lambda l: l[len(pre):len(l) - len(suf)] or l
for k, i in enumerate([i for i in np.argsort(-np.abs(ratio)) if (flag[i] or smk[i]) and inf[i]][:a.label]):
    ax.annotate(short(loci[i]), (x[i], y[i]), xytext=(6, -10 - 9 * k), textcoords="offset points", fontsize=6, color=INK, arrowprops=dict(arrowstyle="-", color=FLOOR, lw=0.5))
ax.set_xlabel("expected locus divergence, median over comparators (%)", fontsize=8, color=MUTED); ax.set_ylabel("observed locus divergence, median over comparators (%)", fontsize=8, color=MUTED)
TITLE = dict(fontsize=10, color=INK, loc="left", fontweight="bold")
ax.set_title(f"{a.taxon}  ·  {os.path.basename(os.path.normpath(a.dir))}", pad=16, **TITLE)
sub = (f"{len(loci)} loci ({int((~inf).sum())} with no difference expected or seen) · closest comparator at m = {np.nanmedian(mclose):.2f}% (median over loci) · ",
       f"{int(flag.sum())} masked too divergent, {int(smk.sum())} too similar" + (f", {int(amb.sum())} ambiguous" if amb.any() else "")
       + (f" · {nshort} loci without a failing test not drawn (the report is sparse)" if nshort else ""))
subt = ax.text(0, 1.012, "".join(sub), transform=ax.transAxes, fontsize=7, color=MUTED, va="bottom")
bins = np.linspace(min(ratio[inf].min(), -2.5), max(ratio[inf].max(), 2.5), 41)
# the nine fills of the scatter, stacked from the base: reds, blues, greys; within each colour from fewest pairs to most;
# a sequence masked for both reasons counts once, with the blues (its blue mark is drawn over its red one)
groups = [(flag & ~smk & (nf == k) & inf, REDS[k]) for k in sorted(FAILS)] + \
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
          [Line2D([], [], marker="o", ls="", mfc="none", mec=SIMRING, ms=6, label="in a too-similar pair, not masked"),
           Line2D([], [], marker="o", ls="", mfc="none", mec=AMB, ms=6, label="ambiguous"), Line2D([], [], color=MODEL, lw=1.3, label="expected")]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, -0.04), labelcolor=INK, columnspacing=1.2, handletextpad=0.4)
for axx in (ax, hx):
    for s_ in ("top", "right"): axx.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"): axx.spines[s_].set_color(SPINE)
    axx.tick_params(colors=MUTED, labelsize=7, length=0); axx.grid(True, color=GRID, lw=0.5); axx.set_axisbelow(True)
ax.tick_params(which="minor", labelsize=7, labelcolor=MUTED)
fig.tight_layout(rect=(0, 0.07, 1, 1))
rr = fig.canvas.get_renderer()
if subt.get_window_extent(rr).x1 > hx.get_window_extent(rr).x0 - 8:      # it would run into the histogram's title (left-aligned on its axes): two lines
    subt.set_text(sub[0].rstrip(" ·") + "\n" + sub[1]); ax.set_title(ax.get_title(loc="left"), pad=25, **TITLE); fig.tight_layout(rect=(0, 0.07, 1, 1))
out = a.out or f"{a.dir}/taxon_{a.taxon}.png"
fig.savefig(out, bbox_inches="tight"); print("wrote", out)
