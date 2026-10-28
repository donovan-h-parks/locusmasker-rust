#!/usr/bin/env python3
"""One taxon pair, every locus: observed against expected divergence, outliers highlighted.

    python plot_pair.py DIR TAXON_A TAXON_B [--out FILE] [--label N]

Reads tests.tsv.gz, which locus_outliers.py writes only with --all-tests: every test of the run, so
the pair's too-similar test at every locus where the two share >= --fit-overlap sites, and its
too-divergent test where one is among the other's comparators. Also masked.csv, and report.npz
for the pair's m. Left: one point per locus, x = the divergence the locus expects for this pair
(its rate at the pair's m; the too-similar test holds the pair to this times the row factor, at
most 1), y = what it shows; the diagonal is the expectation. Right: the same ratios as a
histogram. Colours are the pair's own verdict at each locus: red = too divergent, blue = too
similar, grey = within expectation; the subtitle says which of the two is the copy where a too-
similar pair's direction was identified. A dark ring marks loci where TAXON_A's cell was masked
(too divergent against its comparators, or the copy in a too-similar pair). The --label N most
extreme loci are named."""
import argparse, csv, gzip, json, math, os, zipfile, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--out"); ap.add_argument("--label", type=int, default=6)
a = ap.parse_args()
if a.a == a.b: raise SystemExit("give two different taxa")
tp = os.path.join(a.dir, "tests.tsv.gz")
if not os.path.exists(tp):
    raise SystemExit(f"{tp} not found: a pair's tests at every locus are in tests.tsv.gz, which locus_outliers.py writes only with --all-tests. "
                     f"Re-run it with --all-tests (the same --aln and options, a new --out) and give this script that directory")
rp = os.path.join(a.dir, "report.npz"); Z = zipfile.ZipFile(rp) if os.path.exists(rp) else None
if Z is not None:                                                    # name a taxon the run does not have before reading every test
    names = json.loads(Z.read("taxa.json"))["names"]
    for x in (a.a, a.b):
        if x not in names: raise SystemExit(f"{x} is not a taxon of the run in {a.dir}")
# every test of the pair, by locus: one row per test (divergent, similar, direction); the file is read as bytes and a line
# is split only when it holds the two names side by side, so a large file streams quickly
rows = {}; keys = {(a.a, a.b), (a.b, a.a)}; pats = tuple(f"\t{x}\t{y}\t".encode() for x, y in keys)
with gzip.open(tp, "rb") as fh:
    head = fh.readline().decode().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(head)}
    for line in fh:
        if pats[0] in line or pats[1] in line:
            f = line.decode().rstrip("\n").split("\t")
            if (f[c["taxon_a"]], f[c["taxon_b"]]) in keys: rows.setdefault(f[c["locus"]], {})[f[c["test"]]] = f
if not rows: raise SystemExit(f"{a.a} and {a.b} were never tested against each other: they share fewer than --fit-overlap sites at every locus, "
                              f"or one of them is not in the test (excluded, or not a taxon of the run)")
masked = {}
if os.path.exists(f"{a.dir}/masked.csv"):
    for r in csv.DictReader(open(f"{a.dir}/masked.csv")):
        if r["taxon"] == a.a: masked[r["locus"]] = r["reason"]
loci = sorted(rows); R = [rows[l] for l in loci]
num = lambda d, k: float((d.get("similar") or d["divergent"])[c[k]])   # the too-similar row where there is one: it exists wherever the pair is testable
o = np.array([num(d, "overlap") for d in R]); obs = 100.0 * np.array([num(d, "mismatches") for d in R]) / o
exp = 100.0 * np.array([num(d, "expected") for d in R]) / o            # expected mismatches at the pair's m (before the row factor)
hi = np.array(["divergent" in d and d["divergent"][c["fail"]] == "1" for d in R]); lo = np.array(["similar" in d and d["similar"][c["fail"]] == "1" for d in R])
copy = np.array([d["direction"][c["recipient"]] if "direction" in d else "" for d in R], dtype=object)   # the identified copy of a too-similar pair
ring = np.array([l in masked for l in loci])
def pair_m(kind):
    """the pair's genome-wide divergence: the test's expected divergence at each locus taken back through the locus's
    expectation curve, P (1 - exp(-rate c(m) / P)), c(p) = -P ln(1 - p / P); the median over loci. The too-similar test
    uses m, the too-divergent test max(m, m_res), the pair's resolution floor."""
    if Z is None: return None
    LJ = json.loads(Z.read("loci.json")); at = {n: i for i, n in enumerate(LJ["names"])}; v = []
    for l, d in zip(loci, R):
        if kind not in d or LJ["test_rate"][at[l]] is None: continue
        P_, r_ = float(LJ["plateau"][at[l]]), float(LJ["test_rate"][at[l]]); p_ = 100.0 * float(d[kind][c["expected"]]) / float(d[kind][c["overlap"]])
        v.append(P_ * (1.0 - math.exp(P_ * math.log(1.0 - min(p_, 0.98 * P_) / P_) / (r_ * P_))))
    return float(np.median(v)) if v else None
m, m_hi = pair_m("similar"), pair_m("divergent")
mtxt = ("" if m is None and m_hi is None else f"pair divergence m = {m if m is not None else m_hi:.2f}%"
        + (f" (the too-divergent test at its floor, {m_hi:.2f}%)" if m is not None and m_hi is not None and m_hi > 1.005 * m else "") + " · ")
ncopy = [int(((copy == x) & lo).sum()) for x in (a.a, a.b)]; nund = int(lo.sum()) - sum(ncopy)
dirtxt = ", ".join([f"{x} the copy at {k}" for x, k in zip((a.a, a.b), ncopy) if k] + ([f"no direction at {nund}"] if nund else []))
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
ax.plot(g, g, color=MODEL, lw=1.3, zorder=2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_aspect("equal")
for axis in (ax.xaxis, ax.yaxis):                                       # minor ticks labelled only when no power of 10 is in view
    axis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(0, 0.4)))
pre = os.path.commonprefix(loci); suf = os.path.commonprefix([l[::-1] for l in loci])[::-1]   # the part of a name that varies
short = lambda l: l[len(pre):len(l) - len(suf)] or l
for k, i in enumerate([i for i in np.argsort(-np.abs(ratio)) if (hi[i] or lo[i]) and inf[i]][:a.label]):
    ax.annotate(short(loci[i]), (x[i], y[i]), xytext=(6, -10 - 9 * k), textcoords="offset points",
                fontsize=6, color=INK, arrowprops=dict(arrowstyle="-", color=FLOOR, lw=0.5))
ax.set_xlabel("expected locus divergence for this pair (%)", fontsize=8, color=MUTED); ax.set_ylabel("observed locus divergence (%)", fontsize=8, color=MUTED)
sub = (f"{len(loci)} loci ({int((~inf).sum())} with no difference expected or seen) · {mtxt}",
       f"{hi.sum()} too divergent, {lo.sum()} too similar" + (f" ({dirtxt})" if dirtxt else "") + f" · {ring.sum()} of {a.a}'s cells masked")
TITLE = dict(fontsize=10, color=INK, loc="left", fontweight="bold")
ax.set_title(f"{a.a} vs {a.b}", pad=16, **TITLE)
subt = ax.text(0, 1.012, "".join(sub), transform=ax.transAxes, fontsize=7, color=MUTED, va="bottom")
bins = np.linspace(min(ratio[inf].min(), -2.5), max(ratio[inf].max(), 2.5), 41)
groups = [(lo & inf, BLUE), (hi & inf, RED), (~hi & ~lo & inf, GREY)]          # stacked, outliers at the base and grey on top; log count axis so a few loci still show
hx.hist([ratio[sel] for sel, _ in groups], bins=bins, color=[c for _, c in groups], stacked=True, edgecolor=SURF, linewidth=0.4, zorder=2)
hx.set_yscale("log"); hx.set_ylim(0.8, None)
hx.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0))); hx.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
hx.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(round(v)):,}" if v >= 1 else ""))
hx.set_title(f"{int(inf.sum())} loci with a difference expected or seen · log count axis", fontsize=7, color=MUTED, loc="left", pad=4)
hx.axvline(0, color=MODEL, lw=1.3)
hx.set_xlabel("log2 (observed / expected)", fontsize=8, color=MUTED); hx.set_ylabel("loci", fontsize=8, color=MUTED)
handles = [Line2D([], [], marker="o", ls="", color=GREY, ms=5, label="within expectation"),
           Line2D([], [], marker="o", ls="", color=RED, ms=5, label="too divergent (pair fails)"),
           Line2D([], [], marker="o", ls="", color=BLUE, ms=5, label="too similar (pair fails)"),
           Line2D([], [], marker="o", ls="", mfc="none", mec=INK, ms=6, label=f"{a.a}'s cell masked"),
           Line2D([], [], color=MODEL, lw=1.3, label="expected")]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.01), labelcolor=INK)
for axx in (ax, hx):
    for s in ("top", "right"): axx.spines[s].set_visible(False)
    for s in ("left", "bottom"): axx.spines[s].set_color(SPINE)
    axx.tick_params(colors=MUTED, labelsize=7, length=0); axx.grid(True, color=GRID, lw=0.5); axx.set_axisbelow(True)
ax.tick_params(which="minor", labelsize=7, labelcolor=MUTED)
fig.tight_layout(rect=(0, 0.05, 1, 1))
rr = fig.canvas.get_renderer()
if subt.get_window_extent(rr).x1 > hx.get_window_extent(rr).x0 - 8:      # it would run into the histogram's title (left-aligned on its axes): two lines
    subt.set_text(sub[0].rstrip(" ·") + "\n" + sub[1]); ax.set_title(f"{a.a} vs {a.b}", pad=25, **TITLE); fig.tight_layout(rect=(0, 0.05, 1, 1))
out = a.out or f"{a.dir}/pair_{a.a}_vs_{a.b}.png"
fig.savefig(out, bbox_inches="tight"); print("wrote", out)
