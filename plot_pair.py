#!/usr/bin/env python3
"""One taxon pair, every locus where the two were tested against each other: observed against expected divergence,
the loci where the pair failed a test highlighted.

    python plot_pair.py DIR TAXON_A TAXON_B [--out FILE] [--label N]

Reads tests.tsv.gz, which locusmasker.py writes only with --all-tests: every test of the run. A pair is tested
at a locus when one is among the other's comparators there (the too-divergent test and the too-similar comparator
test, in both directions) or one is the other's closest match at the locus outside its comparators (the too-similar
closest-match test). Also masked.csv (which of the two is the copy where a too-similar pair's direction was called,
and TAXON_A's masked cells) and report.npz for the pair's m. Left: one point per locus, x = the divergence the test
expected for the pair (the comparator tests: the locus rate at the pair's m, before the row factor; the closest-match
test: what the pair's relatives show), y = what it shows; the diagonal is the expectation. Right: the same ratios as
a histogram. Colours are the pair's verdict at each locus: red = too divergent, blue = too similar, grey = within
expectation; the subtitle says which of the two is the copy where a too-similar pair's direction was called. A dark
ring marks loci where TAXON_A's cell was masked. The --label N most extreme loci are named."""
import argparse, csv, gzip, json, math, os, zipfile, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--out"); ap.add_argument("--label", type=int, default=6)
a = ap.parse_args()
if a.a == a.b: raise SystemExit("give two different taxa")
tp = os.path.join(a.dir, "tests.tsv.gz")
if not os.path.exists(tp):
    raise SystemExit(f"{tp} not found: a pair's tests at every locus are in tests.tsv.gz, which locusmasker.py writes only with --all-tests. "
                     f"Re-run it with --all-tests (the same --aln and options, a new --out) and give this script that directory")
rp = os.path.join(a.dir, "report.npz"); Z = zipfile.ZipFile(rp) if os.path.exists(rp) else None
if Z is not None:                                                    # name a taxon the run does not have before reading every test
    names = json.loads(Z.read("taxa.json"))["names"]
    for x in (a.a, a.b):
        if x not in names: raise SystemExit(f"{x} is not a taxon of the run in {a.dir}")
# every test of the pair, by locus and test (close_relative, closest_match_relative, closest_match_outside; one row per direction
# tested); the file is read as bytes and a line is split only when it holds the two names side by side, so a large
# file streams quickly
COLS = ("test", "locus", "taxon", "comparator", "overlap", "mismatches", "expected", "row_factor", "p", "fail")
rows = {}; keys = {(a.a, a.b), (a.b, a.a)}; pats = tuple(f"\t{x}\t{y}\t".encode() for x, y in keys)
with gzip.open(tp, "rb") as fh:
    head = fh.readline().decode().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(head)}
    missing = [k for k in COLS if k not in c]
    if missing: raise SystemExit(f"{tp} has no column {', '.join(missing)}: it was written by another version of locusmasker.py (its columns: {', '.join(head)})")
    for line in fh:
        if pats[0] in line or pats[1] in line:
            f = line.decode().rstrip("\n").split("\t")
            if (f[c["taxon"]], f[c["comparator"]]) in keys: rows.setdefault(f[c["locus"]], {}).setdefault(f[c["test"]], []).append(f)
if not rows: raise SystemExit(f"{a.a} and {a.b} were never tested against each other: at no locus is one among the other's comparators or its closest match "
                              f"(they share fewer than --fit-overlap sites, or are not close relatives), or one of them is not in the test (excluded, or not a taxon of the run)")
masked = {}; copy = {}                                                 # TAXON_A's masked cells by locus; the copy of the pair where a direction was called
if os.path.exists(f"{a.dir}/masked.csv"):
    with open(f"{a.dir}/masked.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["taxon"] == a.a: masked[r["locus"]] = r["reason"]
            partner = r.get("closest_match_partner"); role = r.get("closest_match_role")
            if r["taxon"] in (a.a, a.b) and partner in (a.a, a.b) and role in ("copy", "uncertain"):
                copy[r["locus"]] = r["taxon"] if role == "copy" else "?"
loci = sorted(rows); R = [rows[l] for l in loci]
def test_row(d):
    """the row whose expectation the point shows: the too-similar comparator test (at m, before the row factor), else the
    too-divergent test, else the closest-match test"""
    return (d.get("closest_match_relative") or d.get("close_relative") or d["closest_match_outside"])[0]
num = lambda d, k: float(test_row(d)[c[k]])
o = np.array([num(d, "overlap") for d in R]); obs = 100.0 * np.array([num(d, "mismatches") for d in R]) / o
exp = 100.0 * np.array([num(d, "expected") for d in R]) / o
failed = lambda d, k: any(f[c["fail"]] == "1" for f in d.get(k, []))
hi = np.array([failed(d, "close_relative") for d in R]); lo = np.array([failed(d, "closest_match_relative") or failed(d, "closest_match_outside") for d in R])
copy_of = np.array([copy.get(l, "") for l in loci], dtype=object)
ring = np.array([l in masked for l in loci])
def pair_m(kind):
    """the pair's genome-wide divergence: the test's expected divergence at each locus taken back through the locus's
    expectation curve, P (1 - exp(-rate c(m) / P)), c(p) = -P ln(1 - p / P); the median over loci. The too-similar
    comparator test uses m, the too-divergent test max(m, m_res), the pair's resolution floor."""
    if Z is None: return None
    LJ = json.loads(Z.read("loci.json")); at = {n: i for i, n in enumerate(LJ["names"])}; v = []
    for l, d in zip(loci, R):
        if kind not in d or LJ["test_rate"][at[l]] is None: continue
        f = d[kind][0]; P_, r_ = float(LJ["plateau"][at[l]]), float(LJ["test_rate"][at[l]]); p_ = 100.0 * float(f[c["expected"]]) / float(f[c["overlap"]])
        v.append(P_ * (1.0 - math.exp(P_ * math.log(1.0 - min(p_, 0.98 * P_) / P_) / (r_ * P_))))
    return float(np.median(v)) if v else None
m, m_hi = pair_m("closest_match_relative"), pair_m("close_relative")
mtxt = ("" if m is None and m_hi is None else f"pair divergence m = {m if m is not None else m_hi:.2f}%"
        + (f" (the too-divergent test at its floor, {m_hi:.2f}%)" if m is not None and m_hi is not None and m_hi > 1.005 * m else "") + " · ")
ncopy = [int(((copy_of == x) & lo).sum()) for x in (a.a, a.b)]; nunc = int(((copy_of == "?") & lo).sum())
dirtxt = ", ".join([f"{x} the copy at {k}" for x, k in zip((a.a, a.b), ncopy) if k] + ([f"direction uncertain at {nunc}"] if nunc else []))
eps = 0.5 / o * 100                                          # half a mismatch: keeps zero divergences on a log axis
y = np.maximum(obs, eps); x = np.maximum(exp, eps); ratio = np.log2(y / x)
inf = (obs * o / 100 >= 1) | (exp * o / 100 >= 1)              # loci where at least one side expects or shows a difference

INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
GREY, RED, BLUE, MODEL, FLOOR = "#9a9893", "#e34948", "#3987e5", "#3b3a37", "#a5a39d"
col = np.where(hi, RED, np.where(lo, BLUE, GREY))
fig, (ax, hx) = plt.subplots(1, 2, figsize=(10.5, 4.6), dpi=150, gridspec_kw=dict(width_ratios=(3, 2)))
order = [i for i in np.argsort(hi | lo) if inf[i]]              # outliers drawn last; uninformative loci not drawn
if not inf.any(): raise SystemExit(f"{a.a} and {a.b} show no difference and expect none at any of the {len(loci)} loci where they were tested: nothing to plot")
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
sub = (f"{len(loci)} loci where the pair was tested ({int((~inf).sum())} with no difference expected or seen) · {mtxt}",
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
if subt.get_window_extent(rr).x1 > hx.get_window_extent(rr).x0 - 30:     # it would run into the histogram's title (left-aligned on its axes): two lines
    subt.set_text(sub[0].rstrip(" ·") + "\n" + sub[1]); ax.set_title(f"{a.a} vs {a.b}", pad=25, **TITLE); fig.tight_layout(rect=(0, 0.05, 1, 1))
out = a.out or f"{a.dir}/pair_{a.a}_vs_{a.b}.png"
fig.savefig(out, bbox_inches="tight"); print("wrote", out)
