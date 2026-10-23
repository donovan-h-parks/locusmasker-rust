#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF (and PNGs with --png).

    python plot_panels.py DIR [--plots loci,taxa,bars] [--top 100 | --all] [--cols 7] [--out DIR/plots.pdf] [--png] [--title NAME]

Pages
  loci   one panel per locus: every tested cell of the locus (a taxon) as one point
  taxa   one panel per taxon: every tested cell of the taxon (a locus) as one point
  bars   one bar per locus with a masked cell (all of them, whatever --top): too divergent (reds) then too similar (blues)

In the panels a point is a cell, placed at the median over its comparators of the divergence the
locus expects for the pair (x) and of the divergence observed (y) — the two numbers the cell
decision is made on. The diagonal is the expectation; the dashed lines are the effect floors
(2x above, 1/3 below). Shape is the verdict: circles pass, triangles pointing up are masked as too
divergent, triangles pointing down as too similar. Colour repeats the verdict and grades it:
greys by comparators failed (0, 1, 2), reds by comparators failed (3, 4, 5; an ambiguous cell,
masked with fewer comparators than requested, takes the lightest red), blues by too-similar
pairs (1, 2, 3+). Cells with no difference expected or seen sit on the
axis floor marked 0. Loci and taxa are ordered by cells masked; --top N (default 100) shows the N
most masked on each page, --all shows every locus and every taxon (a page per panel row of 7, so
a dataset with thousands of loci or taxa makes a very large PDF). The bars page always lists every
locus that has a masked cell."""
import argparse, csv, collections, math, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--plots", default="taxa,loci,bars")
ap.add_argument("--out"); ap.add_argument("--png", action="store_true"); ap.add_argument("--cols", type=int, default=7)
ap.add_argument("--top", type=int, default=100, help="panels per page: the N loci / taxa with the most masked cells")
ap.add_argument("--all", action="store_true", help="plot every locus and every taxon, however many (overrides --top)")
ap.add_argument("--title", help="dataset name for the page headers (default: the directory name)")
a = ap.parse_args(); out = a.out or f"{a.dir}/plots.pdf"
NAME = a.title or os.path.basename(os.path.normpath(a.dir))
wanted = [p.strip() for p in a.plots.split(",") if p.strip()]
bad = [p for p in wanted if p not in ("loci", "taxa", "bars")]
if bad: raise SystemExit(f"unknown plot(s): {bad}; choose from loci, taxa, bars")

# ---- data: one record per tested cell
side = {}
if os.path.exists(f"{a.dir}/masked.csv"):
    for r in csv.DictReader(open(f"{a.dir}/masked.csv")): side[(r["locus"], r["taxon"])] = r["side"]
R = list(csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t"))
if not R: raise SystemExit("nothing to plot: no cell was tested")
per = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/per_locus.tsv"), delimiter="\t")}
FAILS = sorted(int(k[4:]) for k in next(iter(per.values())) if re.match(r"fail\d+$", k)) or [3, 4, 5]
NC, LOW = FAILS[-1], FAILS[0]
LOC = np.array([r["locus"] for r in R]); TAX = np.array([r["taxon"] for r in R])
EXP = np.array([100.0 - float(r["med_expected"]) for r in R]); OBS = np.array([100.0 - float(r["med_identity"]) for r in R])
O = np.array([int(r["overlap"]) for r in R]); NF = np.array([int(r["n_fail"]) for r in R]); NS = np.array([int(r["n_similar"]) for r in R])
FLAG = np.array([r["flag"] == "1" for r in R]); AMB = np.array([r["ambiguous"] == "1" for r in R])
SMK = np.array([side.get((l, t)) in ("similar", "both") for l, t in zip(LOC, TAX)])
DMK = FLAG | AMB; PLAIN = ~DMK & ~SMK; MASKED = DMK | SMK
INF = (OBS * O / 100 >= 1) | (EXP * O / 100 >= 1)                 # a difference expected or seen

# ---- tokens. Text wears ink; marks wear the verdict colours. The three ramps are light greys,
# dark reds and the reference blues so that the families stay apart by lightness as well as hue
# (worst cross-family OKLab distance 6.2 under protanopia, 9.7 under deuteranopia); shape carries
# the verdict independently of colour.
SURF, INK, MUTED, FAINT = "#ffffff", "#0b0b0b", "#52514e", "#8c8a85"
GRID, SPINE, MODEL, FLOOR = "#ebeae6", "#d4d2cc", "#3b3a37", "#a5a39d"
GREY3 = ["#b9b7b2", "#918f8a", "#6a6862"]; RED3 = ["#d9433f", "#a8201f", "#6f0d0f"]; BLUE3 = ["#86b6ef", "#3987e5", "#184f95"]
def ramp(cols, n):
    if n <= 3: return cols[:n] if n < 3 else cols
    u = np.linspace(0, 2, n); return [matplotlib.colors.to_hex(np.array(matplotlib.colors.to_rgb(cols[int(v)])) * (1 - v % 1) + np.array(matplotlib.colors.to_rgb(cols[min(int(v) + 1, 2)])) * (v % 1)) for v in u]
GREYS = dict(zip(range(LOW), ramp(GREY3, LOW))); REDS = dict(zip(FAILS, ramp(RED3, len(FAILS)))); BLUES = {1: BLUE3[0], 2: BLUE3[1], 3: BLUE3[2]}

def fit_title(name, width_in, fs=7.5):
    cap = int(width_in * 72 / (0.56 * fs))
    if len(name) <= cap: return name, fs
    seps = [i for i, ch in enumerate(name) if ch in "_-. |" and 0.25 * len(name) < i < 0.75 * len(name)]
    i = min(seps, key=lambda i: abs(i - len(name) / 2)) + 1 if seps else len(name) // 2
    lines = [name[:i], name[i:]]; longest = max(len(l) for l in lines)
    return "\n".join(lines), fs if longest <= cap else max(5.5, fs * cap / longest)

def log_ticks(lo, hi):
    cand = [100, 70, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5, 1, 0.7, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.01, 0.005, 0.002, 0.001]
    t = [v for v in cand if lo <= v <= hi]
    while len(t) > 3: t = t[::2]
    return t

def panels(by, free=False):
    """one panel per locus (by='locus') or per taxon (by='taxon'); free = x and y ranged separately"""
    key = LOC if by == "locus" else TAX
    counts = collections.Counter(key[MASKED]); groups = collections.defaultdict(list)
    for i, k in enumerate(key): groups[k].append(i)
    order = sorted(groups, key=lambda k: (-counts[k], k)); NTOT = len(order)
    if not a.all: order = order[:a.top]
    ncol = min(a.cols, len(order)); nrow = math.ceil(len(order) / ncol); pw, ph = 2.3, 2.25
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 0.7, ph * nrow + 1.35), dpi=110, squeeze=False)
    for ax, k in zip(axs.ravel(), order):
        idx = np.array(groups[k]); inf = INF[idx]
        # the panel's own range: the central 99.8% of its informative values with a margin, the same on
        # both axes so the diagonal stays at 45 degrees; cells with a zero sit on a floor tick labelled 0
        # the axes are scaled by matplotlib from the plotted points alone (free: x and y separately;
        # otherwise one square range); a zero divergence cannot sit on a log axis, so such cells go
        # on a floor at half the smallest positive value, ticked 0
        px, py = EXP[idx][inf], OBS[idx][inf]
        zx, zy = bool((px <= 0).any()), bool((py <= 0).any())
        fx = (px[px > 0].min() / 2 if (px > 0).any() else 0.01); fy = (py[py > 0].min() / 2 if (py > 0).any() else 0.01)
        if not free: fx = fy = min(fx, fy)
        X = np.where(EXP > 0, EXP, fx); Y = np.where(OBS > 0, OBS, fy)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.margins(0.04)
        def draw(sel, marker, color, size, hollow=False, edge=None, z=3):
            s = idx[sel & inf]
            if len(s): ax.scatter(X[s], Y[s], s=size, marker=marker, facecolors="none" if hollow else color, edgecolors=edge or (SURF if not hollow else color),
                                  linewidths=0.8 if hollow else 0.5, zorder=z, rasterized=not hollow and marker == "o")
        for j, c in GREYS.items(): draw(PLAIN[idx] & (NF[idx] == j), "o", c, 12, z=2)
        draw(AMB[idx] & ~FLAG[idx], "^", REDS[LOW], 26, z=5)
        for j, c in REDS.items(): draw(FLAG[idx] & (NF[idx] == j), "^", c, 26, z=5)
        for j, c in BLUES.items(): draw(SMK[idx] & (np.minimum(NS[idx], 3) == j), "v", c, 26, z=6)
        ax.autoscale_view()                                           # from the scatters' data limits (relim would drop collections)
        (xlo, xhi), (ylo, yhi) = ax.get_xlim(), ax.get_ylim()
        if not free: xlo = ylo = min(xlo, ylo); xhi = yhi = max(xhi, yhi)
        ax.set_autoscale_on(False)
        g = np.array([min(xlo, ylo) / 10, max(xhi, yhi) * 10])
        ax.plot(g, g, color=MODEL, lw=1.1, zorder=1); ax.plot(g, 2 * g, color=FLOOR, lw=0.8, ls=(0, (3, 2)), zorder=1); ax.plot(g, g / 3, color=FLOOR, lw=0.8, ls=(0, (3, 2)), zorder=1)
        nm = int(MASKED[idx].sum()); n = len(idx); hot = nm >= 40
        title, fs = fit_title(k, pw - 0.25)
        ax.set_title(title, fontsize=fs, color=INK, fontweight="bold" if hot else "normal", pad=12, linespacing=0.95)
        ax.text(0.5, 1.012, f"{nm} of {n} masked" + (f" · {int((~inf).sum())} with no differences" if (~inf).sum() else ""), transform=ax.transAxes,
                fontsize=6.3, color=INK if nm else FAINT, ha="center", va="bottom")
        ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
        if not free: ax.set_aspect("equal")
        tx = log_ticks(xlo, xhi); ty = log_ticks(ylo, yhi)
        ax.set_xticks(tx + ([fx] if zx else [])); ax.set_xticklabels([f"{v:g}" for v in tx] + (["0"] if zx else []))
        ax.set_yticks(ty + ([fy] if zy else [])); ax.set_yticklabels([f"{v:g}" for v in ty] + (["0"] if zy else []))
        for axis in (ax.xaxis, ax.yaxis): axis.set_minor_formatter(matplotlib.ticker.NullFormatter()); axis.set_minor_locator(matplotlib.ticker.NullLocator())
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"): ax.spines[s_].set_color(SPINE); ax.spines[s_].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, labelsize=6, length=0, pad=1.5); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for ax in axs.ravel()[len(order):]: ax.axis("off")
    fig.supxlabel("expected divergence of the cell, median over its comparators (%)", fontsize=8, color=MUTED, y=0.004)
    fig.supylabel("observed divergence of the cell, median over its comparators (%)", fontsize=8, color=MUTED, x=0.004)
    H = ph * nrow + 1.35
    unit = "loci" if by == "locus" else "taxa"
    head = f"{NAME}  ·  {int(MASKED.sum()):,} of {len(R):,} cells masked ({int(DMK.sum()):,} too divergent, {int(SMK.sum()):,} too similar)" + (f"  ·  {len(order)} of {NTOT} {unit} shown" if len(order) < NTOT else f"  ·  one panel per {by}")
    sub = ("one point = a cell (" + ("a taxon at this locus" if by == "locus" else "a locus of this taxon") + "), at the median over its comparators of expected and observed divergence  ·  "
           "diagonal = each cell's own expectation  ·  dashed = the effect floors applied to each cell (2x above, 1/3 below)  ·  no line is fitted to a panel  ·  shape = verdict, colour grades it")
    fig.text(0.01, 1 - 0.14 / H, head, fontsize=11, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 1 - 0.36 / H, sub, fontsize=7.5, color=MUTED, ha="left", va="top")
    handles = [Line2D([], [], marker="o", ls="", color=c, ms=4.5, label=("passes · " if j == 0 else "") + f"fails {j} of {NC}") for j, c in GREYS.items()] + \
              [Line2D([], [], marker="^", ls="", color=c, mec=SURF, ms=6, label=("masked, too divergent · " if j == LOW else "") + f"fails {j} of {NC}") for j, c in REDS.items()] + \
              [Line2D([], [], marker="v", ls="", color=BLUES[j], mec=SURF, ms=6, label=("masked, too similar · " if j == 1 else "") + f"{j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}") for j in (1, 2, 3)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.52 / H), ncol=min(len(handles), max(4, int((pw * ncol + 0.7) / 1.5))),
               frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.2, labelcolor=INK)
    fig.tight_layout(rect=(0.012, 0.012, 1, 1 - 1.0 / H), h_pad=1.0, w_pad=0.5)
    return fig

def bars():
    counts = collections.Counter(LOC[MASKED])
    order = sorted(counts, key=lambda k: (-counts[k], k)); NTOT = len(set(LOC))      # every locus with at least one masked cell
    if not order: order = sorted(set(LOC))[:1]
    n = len(order); ys = np.arange(n)[::-1]
    lab_in = min(0.048 * max(len(l) for l in order) + 0.4, 4.5)
    fig, ax = plt.subplots(figsize=(7.0 + lab_in, 0.17 * n + 1.6), dpi=200)
    right = [(f"fails {j} of {NC}", REDS[j], np.array([int(((LOC == k) & ((FLAG & (NF == j)) | ((AMB & ~FLAG) if j == LOW else False))).sum()) for k in order])) for j in FAILS]
    left = [(f"{j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}", BLUES[j], np.array([int(((LOC == k) & SMK & (np.minimum(NS, 3) == j)).sum()) for k in order])) for j in (1, 2, 3)]
    base = np.zeros(n)
    for lab, c, v in right + left:                                   # one stacked bar: too divergent first, then too similar
        nz = v > 0; ax.barh(ys[nz], v[nz], left=base[nz], color=c, edgecolor=SURF, linewidth=0.8, height=0.72); base += v
    top = max(base.max(), 1)
    for y, t_ in zip(ys, base):
        if t_: ax.text(t_ + top * 0.012, y, f"{int(t_):,}", va="center", ha="left", fontsize=6, color=MUTED)
    ax.set_yticks(ys); ax.set_yticklabels(order, fontsize=6); ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, top * 1.12)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_xlabel("cells masked", fontsize=7.5, color=MUTED)
    ax.set_title(f"{NAME}  ·  {int(MASKED.sum()):,} cells masked, by locus  ·  {n} of {NTOT} loci have a masked cell", fontsize=9, color=INK, loc="left", fontweight="bold", pad=8)
    handles = [Patch(facecolor=c, edgecolor=SURF, label=("too divergent · " if lab.startswith(f"fails {LOW}") else "") + lab) for lab, c, _ in right] + \
              [Patch(facecolor=c, edgecolor=SURF, label=("too similar · " if lab.startswith("1 ") else "") + lab) for lab, c, _ in left]
    ax.legend(handles=handles, frameon=False, fontsize=6.5, loc="lower right", labelcolor=INK, ncol=2)
    for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color(SPINE); ax.tick_params(colors=MUTED, length=0, labelsize=6.5); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    fig.tight_layout()
    return fig

makers = {"loci": lambda: panels("locus", free=True), "taxa": lambda: panels("taxon"), "bars": bars}
with PdfPages(out) as pdf:
    for name in wanted:
        fig = makers[name]()
        pdf.savefig(fig, bbox_inches="tight")
        if a.png: fig.savefig(f"{os.path.splitext(out)[0]}_{name}.png", bbox_inches="tight")
        plt.close(fig); print("page", name)
print("wrote", out)
