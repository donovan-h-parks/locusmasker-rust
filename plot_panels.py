#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF (and PNGs with --png).

    python plot_panels.py DIR [--plots taxa,loci,bars,residuals_taxa,residuals_loci] [--top 100 | --all] [--cols 7] [--out DIR/plots.pdf] [--png] [--title NAME]

Pages
  taxa            one panel per taxon: every tested cell of the taxon (a locus) as one point
  loci            one panel per locus: every tested cell of the locus (a taxon) as one point
  bars            two histograms, one above the other: loci, and taxa, by the number of masked sequences they hold,
                  zero included; each bar is coloured by the mix of mask kinds among those units' masks; a facet's
                  y axis is broken (and says so) when its zero bar dwarfs the rest
  residuals_taxa  one strip per taxon: every cell's extra divergence beyond expectation, as corrected distance in
                  units of m (a typical locus's divergence), on an inverse-hyperbolic-sine axis (linear near 0,
                  logarithmic far out, no joint); masked cells as triangles, the rest grey; cells identical to their
                  comparators are counted in the row label rather than drawn (their residual is just minus the
                  expectation). Cells with one cause share one value, so a batch of bad gene models or one source
                  of contamination is a tight cluster and a scatter is many independent problems.
  residuals_loci  the same, one strip per locus; the shaded band at the right of a row is where a saturated
                  sequence lands, so a tight cluster inside it is foreign sequence, not a shared paralogue.

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

ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--plots", default="taxa,loci,bars,residuals_taxa,residuals_loci")
ap.add_argument("--out"); ap.add_argument("--png", action="store_true"); ap.add_argument("--cols", type=int, default=7)
ap.add_argument("--top", type=int, default=100, help="panels per page: the N loci / taxa with the most masked sequences")
ap.add_argument("--all", action="store_true", help="plot every locus and every taxon, however many (overrides --top)")
ap.add_argument("--title", help="dataset name for the page headers (default: the directory name)")
a = ap.parse_args(); out = a.out or f"{a.dir}/plots.pdf"
NAME = a.title or os.path.basename(os.path.normpath(a.dir))
wanted = [p.strip() for p in a.plots.split(",") if p.strip()]
bad = [p for p in wanted if p not in ("loci", "taxa", "bars", "residuals_taxa", "residuals_loci")]
if bad: raise SystemExit(f"unknown plot(s): {bad}; choose from taxa, loci, bars, residuals_taxa, residuals_loci")

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
RATE = np.array([float(per[l]["rate"]) if per[l]["rate"] not in ("nan", "") else np.nan for l in LOC]); PLAT = np.array([float(per[l]["plateau"]) for l in LOC])
def _corr(p, P): return -P * np.log(1.0 - np.minimum(np.maximum(p, 0.0), 0.98 * P) / P)
with np.errstate(invalid="ignore", divide="ignore"):
    EXM = (_corr(OBS, PLAT) - _corr(EXP, PLAT)) / RATE                # extra divergence beyond expectation, in units of m
    CEIL = (_corr(0.98 * PLAT, PLAT) - _corr(EXP, PLAT)) / RATE        # where a saturated sequence would land

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

def panels(by):
    """one panel per locus (by='locus') or per taxon (by='taxon'); each panel's axes autoscale from its points"""
    key = LOC if by == "locus" else TAX
    counts = collections.Counter(key[MASKED]); groups = collections.defaultdict(list)
    for i, k in enumerate(key): groups[k].append(i)
    order = sorted(groups, key=lambda k: (-counts[k], k)); NTOT = len(order)
    if not a.all: order = order[:a.top]
    ncol = min(a.cols, len(order)); nrow = math.ceil(len(order) / ncol); pw, ph = 2.3, 2.25
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 0.7, ph * nrow + 1.5), dpi=110, squeeze=False)
    for ax, k in zip(axs.ravel(), order):
        idx = np.array(groups[k]); inf = INF[idx]
        # the panel's own range: the central 99.8% of its informative values with a margin, the same on
        # both axes so the diagonal stays at 45 degrees; cells with a zero sit on a floor tick labelled 0
        # the axes are scaled by matplotlib from every plotted point, x and y separately; a zero divergence
        # cannot sit on a log axis, so such cells go on a floor at half the smallest positive value, ticked 0
        px, py = EXP[idx][inf], OBS[idx][inf]
        zx, zy = bool((px <= 0).any()), bool((py <= 0).any())
        fx = (px[px > 0].min() / 2 if (px > 0).any() else 0.01); fy = (py[py > 0].min() / 2 if (py > 0).any() else 0.01)
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
        ax.set_autoscale_on(False)
        g = np.array([min(xlo, ylo) / 10, max(xhi, yhi) * 10])
        ax.plot(g, g, color=MODEL, lw=1.1, zorder=1); ax.plot(g, 2 * g, color=FLOOR, lw=0.8, ls=(0, (3, 2)), zorder=1); ax.plot(g, g / 3, color=FLOOR, lw=0.8, ls=(0, (3, 2)), zorder=1)
        nm = int(MASKED[idx].sum()); n = len(idx)
        title, fs = fit_title(k, pw - 0.25)
        ax.set_title(title, fontsize=fs, color=INK, pad=12, linespacing=0.95)
        ax.text(0.5, 1.012, f"{nm} of {n} sequences masked ({100 * nm / n:.1f}%)", transform=ax.transAxes,
                fontsize=6.3, color=INK if nm else FAINT, ha="center", va="bottom")
        ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
        tx = log_ticks(fx * 10 ** 0.45 if zx else xlo, xhi); ty = log_ticks(fy * 10 ** 0.45 if zy else ylo, yhi)   # keep ticks clear of the 0 floor tick
        ax.set_xticks(tx + ([fx] if zx else [])); ax.set_xticklabels([f"{v:g}" for v in tx] + (["0"] if zx else []))
        ax.set_yticks(ty + ([fy] if zy else [])); ax.set_yticklabels([f"{v:g}" for v in ty] + (["0"] if zy else []))
        for axis in (ax.xaxis, ax.yaxis): axis.set_minor_formatter(matplotlib.ticker.NullFormatter()); axis.set_minor_locator(matplotlib.ticker.NullLocator())
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"): ax.spines[s_].set_color(SPINE); ax.spines[s_].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, labelsize=6, length=0, pad=1.5); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for ax in axs.ravel()[len(order):]: ax.axis("off")
    fig.supxlabel("expected divergence of the cell, median over its comparators (%)", fontsize=8, color=MUTED, y=0.004)
    fig.supylabel("observed divergence of the cell, median over its comparators (%)", fontsize=8, color=MUTED, x=0.004)
    H = ph * nrow + 1.5
    unit = "loci" if by == "locus" else "taxa"
    head = f"{NAME}  ·  {'Locus' if by == 'locus' else 'Taxon'} plot"
    sub = f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)"
    note = ("each point represents a single sequence (a single locus in a single taxon), diagonal line is where the observed is the same as the expected "
            "divergence, dashed lines are the thresholds used for defining outliers (2x above, 1/3 below)" + (f"  ·  the {len(order)} most masked of {NTOT} {unit} are shown" if len(order) < NTOT else ""))
    fig.text(0.01, 1 - 0.14 / H, head, fontsize=11, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 1 - 0.36 / H, sub, fontsize=9, color=INK, ha="left", va="top")
    fig.text(0.01, 1 - 0.55 / H, note, fontsize=7.5, color=MUTED, ha="left", va="top")
    handles = [Line2D([], [], marker="o", ls="", color=c, ms=4.5, label=("passes · " if j == 0 else "") + f"fails {j} of {NC}") for j, c in GREYS.items()] + \
              [Line2D([], [], marker="^", ls="", color=c, mec=SURF, ms=6, label=("masked, too divergent · " if j == LOW else "") + f"fails {j} of {NC}") for j, c in REDS.items()] + \
              [Line2D([], [], marker="v", ls="", color=BLUES[j], mec=SURF, ms=6, label=("masked, too similar · " if j == 1 else "") + f"{j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}") for j in (1, 2, 3)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.70 / H), ncol=min(len(handles), max(4, int((pw * ncol + 0.7) / 1.5))),
               frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.2, labelcolor=INK)
    fig.tight_layout(rect=(0.012, 0.012, 1, 1 - 1.15 / H), h_pad=1.0, w_pad=0.5)
    return fig

def hist_facet(fig, slot, by):
    """one histogram: units (loci or taxa) by the number of masked sequences they hold, zero included, bars
    coloured by the mix of mask kinds among those units' masks; a broken y axis when the zero bar dwarfs the
    rest. Draws into gridspec slot `slot`; returns (top axes or None, bottom axes, legend handles)."""
    key = LOC if by == "locus" else TAX
    nmask = collections.Counter(key[MASKED]); units = sorted(set(key)); NTOT = len(units)
    k_of = {u: nmask.get(u, 0) for u in units}; kmax = max(k_of.values()) if k_of else 0
    if kmax <= 25: bins = [(k, k, str(k)) for k in range(kmax + 1)]
    else:
        edges = [0, 1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10 ** 9]; bins = []
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            hi_i = hi_ - 1 if hi_ > lo_ + 1 else lo_
            if lo_ > kmax: break
            bins.append((lo_, hi_i, str(lo_) if hi_i == lo_ else (f"{lo_}–{min(hi_i, kmax)}" if hi_i < 10 ** 8 else f"{lo_}+")))
    cats = [(f"fails {j} of {NC}", REDS[j], lambda idx, j=j: DMK[idx] & (np.minimum(np.maximum(NF[idx], LOW), NC) == j)) for j in FAILS] + \
           [(f"{j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}", BLUES[j], lambda idx, j=j: SMK[idx] & (np.minimum(NS[idx], 3) == j)) for j in (1, 2, 3)]
    heights = []; comp = []
    for lo_, hi_, lab in bins:
        sel_units = [u for u in units if lo_ <= k_of[u] <= hi_]; n = len(sel_units); heights.append(n)
        if lo_ == 0 or n == 0: comp.append([n] + [0] * len(cats)); continue
        idx = np.flatnonzero(np.isin(key, sel_units) & MASKED)
        shares = np.array([float(f(idx).sum()) for _, _, f in cats]); shares = shares / max(shares.sum(), 1)
        comp.append([0] + list(n * shares))
    comp = np.array(comp); x = np.arange(len(bins))
    big = len(bins) > 1 and heights[0] > 5 * max(heights[1:] + [1])
    if big:
        sub = slot.subgridspec(2, 1, height_ratios=(1, 2.6), hspace=0.08); axt = fig.add_subplot(sub[0]); axb = fig.add_subplot(sub[1], sharex=axt); axes = (axt, axb)
    else:
        axb = fig.add_subplot(slot); axt = None; axes = (axb,)
    cols = [GREY3[0]] + [c for _, c, _ in cats]
    for ax in axes:
        base = np.zeros(len(bins))
        for ci_, c in enumerate(cols):
            v = comp[:, ci_]; nz = v > 0
            ax.bar(x[nz], v[nz], bottom=base[nz], color=c, edgecolor=SURF, linewidth=0.8, width=0.72); base += v
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"): ax.spines[s_].set_color(SPINE)
        ax.tick_params(colors=MUTED, length=0, labelsize=7); ax.yaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    unit = "loci" if by == "locus" else "taxa"
    if big:
        top_lo = max(heights[1:]) * 1.25; axb.set_ylim(0, top_lo); axt.set_ylim(heights[0] * 0.75, heights[0] * 1.12)
        axt.spines["bottom"].set_visible(False); axb.spines["top"].set_visible(False); axt.tick_params(labelbottom=False)
        d = 0.012
        for ax, y0 in ((axt, 0), (axb, 1)):
            ax.plot((-d, +d), (y0 - 2 * d, y0 + 2 * d), transform=ax.transAxes, color=INK, clip_on=False, lw=0.9)
        axt.text(0.0, 0.02, "y axis broken: the top panel is the zero bar only, on its own scale", transform=axt.transAxes, fontsize=6.5, color="#b0171f", va="bottom", ha="left")
        axt.text(x[0], heights[0] * 1.02, f"{heights[0]:,}", ha="center", va="bottom", fontsize=7, color=INK)
        for xi, h in zip(x[1:], heights[1:]):
            if h: axb.text(xi, h + top_lo * 0.01, f"{h:,}", ha="center", va="bottom", fontsize=7, color=INK)
        axt.set_ylabel(unit, fontsize=8, color=MUTED)
    else:
        top = max(heights + [1]); axb.set_ylim(0, top * 1.12)
        for xi, h in zip(x, heights):
            if h: axb.text(xi, h + top * 0.01, f"{h:,}", ha="center", va="bottom", fontsize=7, color=INK)
    axb.set_xticks(x); axb.set_xticklabels([b[2] for b in bins], fontsize=7)
    axb.set_xlabel(f"sequences masked in the {by}", fontsize=8, color=MUTED); axb.set_ylabel(unit, fontsize=8, color=MUTED)
    (axt or axb).set_title(f"by {by}  ·  {len(nmask)} of {NTOT} {unit} have a masked sequence", fontsize=8.5, color=INK, loc="left", pad=6)
    handles = [Patch(facecolor=GREY3[0], edgecolor=SURF, label=f"no masked sequence")] + \
              [Patch(facecolor=c, edgecolor=SURF, label=("too divergent · " if lab.startswith(f"fails {LOW}") else "too similar · " if lab.startswith("1 ") else "") + lab) for lab, c, _ in cats]
    return axt, axb, handles

def bars():
    """the two marginals side by side: loci and taxa by the number of masked sequences they hold"""
    fig = plt.figure(figsize=(8.5, 9.6), dpi=150)
    gs = fig.add_gridspec(2, 1, hspace=0.32, left=0.09, right=0.98, top=0.855, bottom=0.06)
    _, axl, handles = hist_facet(fig, gs[0], "locus"); _, axr, _ = hist_facet(fig, gs[1], "taxon")
    fig.text(0.09, 0.975, f"{NAME}  ·  Number of sequences masked per locus and per taxon", fontsize=10, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.09, 0.95, f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)  ·  bar colour = mix of mask kinds among those loci's or taxa's masks", fontsize=8.5, color=INK, ha="left", va="top")
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.09, 0.93), ncol=4, frameon=False, fontsize=6.5, labelcolor=INK, columnspacing=1.2, handletextpad=0.4)
    return [fig]

ROWS_RESID = 150

def residuals(by):
    """one strip per taxon (by='taxon') or per locus (by='locus'), ROWS_RESID rows per page, most masked first"""
    key = LOC if by == "locus" else TAX
    ok = np.isfinite(EXM) & (RATE > 0)
    counts = collections.Counter(key[MASKED])
    order = sorted(set(key[ok]), key=lambda k: (-counts[k], k)); NTOT = len(order)
    if not a.all: order = order[:a.top]
    pages = [order[i:i + ROWS_RESID] for i in range(0, len(order), ROWS_RESID)]
    return [residual_page(by, chunk, i + 1, len(pages), NTOT, ok) for i, chunk in enumerate(pages)]

def residual_page(by, order, page, npages, NTOT, ok):
    key = LOC if by == "locus" else TAX
    n = len(order); ys = np.arange(n)[::-1]
    pre = os.path.commonprefix(order); suf = os.path.commonprefix([k[::-1] for k in order])[::-1] if n > 1 else ""
    short = (lambda k: k[len(pre):len(k) - len(suf)] or k) if by == "locus" and n > 1 else (lambda k: k)
    labels = []
    fig, ax = plt.subplots(figsize=(10.5, 0.27 * n + 1.8), dpi=150)
    rng = np.random.default_rng(0)
    for y, k in zip(ys, order):
        sel = (key == k) & ok
        ident = sel & (OBS <= 0); drawn = sel & (OBS > 0)
        labels.append(f"{short(k)}   {int(MASKED[sel].sum())} masked · {int(ident.sum()):,} identical")
        jit = rng.uniform(-0.28, 0.28, int(drawn.sum())); idx = np.flatnonzero(drawn)
        plain = PLAIN[idx]; npl = int(plain.sum())
        size, alpha = (10, 0.6) if npl <= 200 else (5, 0.4) if npl <= 1000 else (3, 0.25)
        ax.scatter(EXM[idx][plain], y + jit[plain], s=size, color=GREY3[0], alpha=alpha, linewidths=0, rasterized=True, zorder=1)
        for j, c in REDS.items():
            s_ = DMK[idx] & (np.minimum(np.maximum(NF[idx], LOW), NC) == j)
            if s_.any(): ax.scatter(EXM[idx][s_], y + jit[s_], s=22, marker="^", color=c, edgecolors=SURF, linewidths=0.4, zorder=3)
        for j, c in BLUES.items():
            s_ = SMK[idx] & (np.minimum(NS[idx], 3) == j)
            if s_.any(): ax.scatter(EXM[idx][s_], y + jit[s_], s=22, marker="v", color=c, edgecolors=SURF, linewidths=0.4, zorder=3)
        if by == "locus" and drawn.any():
            c_lo = float(np.percentile(CEIL[idx], 10))
            ax.fill_betweenx([y - 0.4, y + 0.4], c_lo, 1e5, color="#eeedea", zorder=0); ax.plot([c_lo, c_lo], [y - 0.4, y + 0.4], color="#c9c7c2", lw=0.8, zorder=0)
    vals = EXM[np.isin(key, order) & ok & (OBS > 0)]
    xlo = min(float(vals.min()) * 1.3, -30) if len(vals) else -30; xhi = max(float(vals.max()) * 1.3, 300) if len(vals) else 300
    ax.set_xscale("asinh", linear_width=1.0); ax.set_xlim(xlo, xhi)
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([-1000, -100, -10, -1, 0, 1, 10, 100, 1000])); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.axvline(0, color=MODEL, lw=1.0, zorder=2)
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=6.5); ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("extra divergence beyond expectation, in units of a typical locus's divergence (m)   ◀ too similar   |   too divergent ▶", fontsize=8, color=MUTED)
    unit = "taxa" if by == "taxon" else "loci"
    head = f"{NAME}  ·  Residuals by {by}" + (f"  ·  page {page} of {npages}" if npages > 1 else "")
    sub = f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)" + (f"  ·  the {min(a.top, NTOT) if not a.all else NTOT} most masked of {NTOT} {unit}" if not a.all and NTOT > a.top else f"  ·  all {NTOT} {unit}")
    note = ("each point is one sequence; triangles are masked (reds too divergent, blues too similar), grey are the " + by + "'s other sequences\n"
            "sequences identical to their comparators are counted in the row label, not drawn; axis linear near 0, logarithmic far out" + ("; shaded = where a saturated sequence would land" if by == "locus" else ""))
    ax.set_title(f"{head}\n{sub}\n{note}", fontsize=8.5, color=INK, loc="left", linespacing=1.35)
    for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color(SPINE); ax.tick_params(colors=MUTED, length=0, labelsize=6.5); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    fig.tight_layout()
    return fig

makers = {"loci": lambda: [panels("locus")], "taxa": lambda: [panels("taxon")], "bars": bars,
          "residuals_taxa": lambda: residuals("taxon"), "residuals_loci": lambda: residuals("locus")}
with PdfPages(out) as pdf:
    for name in wanted:
        figs = makers[name]()
        for i, fig in enumerate(figs):
            pdf.savefig(fig, bbox_inches="tight")
            if a.png: fig.savefig(f"{os.path.splitext(out)[0]}_{name}{'' if len(figs) == 1 else f'_{i + 1}'}.png", bbox_inches="tight")
            plt.close(fig)
        print("page", name, f"({len(figs)} pages)" if len(figs) > 1 else "")
print("wrote", out)
