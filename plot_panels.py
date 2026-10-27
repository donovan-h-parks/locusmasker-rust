#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF (and PNGs with --png).

    python plot_panels.py DIR [--plots taxa,loci,bars,matches,residuals_taxa,residuals_loci] [--top 100 | --all] [--cols 7] [--out DIR/plots.pdf] [--png] [--title NAME]
                          [--alignment ALN [--partitions FILE] [--datatype dna|aa] [--match-taxa 10]]

Pages
  taxa            one panel per taxon: every tested cell of the taxon (a locus) as one point
  loci            one panel per locus: every tested cell of the locus (a taxon) as one point
  bars            two histograms, one above the other: loci, and taxa, by the number of masked sequences they hold,
                  zero as its own bar and the positive counts in equal-width bins (Freedman-Diaconis width, so
                  width 1 on small datasets); each bar is coloured by the mix of mask kinds among those units'
                  masks; a facet's x axis is broken when empty bins separate a few far units, and its y axis
                  when the zero bar dwarfs the rest (both marked)
  matches         for the 10 taxa with the most masked sequences (--match-taxa): one panel each, a sorted bar chart
                  of the masked sequences by the taxon each matches best at its locus; needs the alignment
                  (--alignment). Identity is counted over the masked sequence's own sites and a site the other
                  taxon lacks counts as a mismatch (so a fragmentary taxon cannot win by lacking the bad part);
                  ties go to the closest relative genome-wide. Bars are stacked by mask kind; the black tick is
                  the count expected if the masked sequences matched like the taxon's unmasked ones; the names
                  of the 5 closest relatives genome-wide are bold. Masked sequences that match the same taxa as
                  the unmasked ones (bars on their ticks) are damaged copies of the right sequence, such as a
                  wrong exon in a gene model; bars shifted onto a taxon that is not a close relative point to
                  contamination from it, or to a shared paralogue.
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
decision is made on. The diagonal is the expectation. Shape is the verdict: circles pass, triangles pointing up are masked as too
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

ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--plots", default="taxa,loci,bars,matches,residuals_taxa,residuals_loci")
ap.add_argument("--out"); ap.add_argument("--png", action="store_true"); ap.add_argument("--cols", type=int, default=7)
ap.add_argument("--top", type=int, default=100, help="panels per page: the N loci / taxa with the most masked sequences")
ap.add_argument("--all", action="store_true", help="plot every locus and every taxon, however many (overrides --top)")
ap.add_argument("--title", help="dataset name for the page headers (default: the directory name)")
ap.add_argument("--alignment", help="the alignment the run was made from (NEXUS, FASTA or PHYLIP); needed for the matches page")
ap.add_argument("--partitions", help="RAxML-style partition file, as given to the tool, when the alignment has no charsets")
ap.add_argument("--datatype", choices=["dna", "aa"], help="as given to the tool, when the alignment does not declare it")
ap.add_argument("--match-taxa", type=int, default=10, help="matches page: the N taxa with the most masked sequences")
a = ap.parse_args(); out = a.out or f"{a.dir}/plots.pdf"
NAME = a.title or os.path.basename(os.path.normpath(a.dir))
wanted = [p.strip() for p in a.plots.split(",") if p.strip()]
bad = [p for p in wanted if p not in ("loci", "taxa", "bars", "matches", "residuals_taxa", "residuals_loci")]
if bad: raise SystemExit(f"unknown plot(s): {bad}; choose from taxa, loci, bars, matches, residuals_taxa, residuals_loci")
if "matches" in wanted and not a.alignment:
    print("page matches skipped: it needs the alignment (--alignment, and --partitions if the alignment has no charsets)"); wanted.remove("matches")

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
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 0.7, ph * nrow + 1.7), dpi=110, squeeze=False)
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
        ax.plot(g, g, color=MODEL, lw=1.1, zorder=1)
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
    H = ph * nrow + 1.7
    unit = "loci" if by == "locus" else "taxa"
    head = f"{NAME}  ·  {'Locus' if by == 'locus' else 'Taxon'} plot"
    sub = f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)"
    note = ("each point represents a single sequence (a single locus in a single taxon), diagonal line is where the observed is the same as the expected "
            "divergence" + (f"  ·  the {len(order)} most masked of {NTOT} {unit} are shown" if len(order) < NTOT else ""))
    fig.text(0.01, 1 - 0.14 / H, head, fontsize=11, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 1 - 0.36 / H, sub, fontsize=9, color=INK, ha="left", va="top")
    fig.text(0.01, 1 - 0.55 / H, note, fontsize=7.5, color=MUTED, ha="left", va="top")
    handles = [Line2D([], [], marker="o", ls="", color=c, ms=4.5, label=f"not masked · fails {j} of {NC}") for j, c in GREYS.items()] + \
              [Line2D([], [], marker="^", ls="", color=c, mec=SURF, ms=6, label=f"masked, too divergent · fails {j} of {NC}") for j, c in REDS.items()] + \
              [Line2D([], [], marker="v", ls="", color=BLUES[j], mec=SURF, ms=6, label=f"masked, too similar · {j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}") for j in (1, 2, 3)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.68 / H), ncol=3, frameon=False, fontsize=7, handletextpad=0.3, columnspacing=2.0, labelcolor=INK)
    fig.tight_layout(rect=(0.012, 0.012, 1, 1 - 1.35 / H), h_pad=1.0, w_pad=0.5)
    return fig

def hist_facet(fig, slot, by):
    """one histogram: units (loci or taxa) by the number of masked sequences they hold. Zero is its own bar;
    the positive counts fall into equal-width integer bins whose width follows the Freedman-Diaconis rule
    (so a small dataset gets width 1 and a large one a handful of bins). Bars are coloured by the mix of
    mask kinds among those units' masks. The x axis is broken (two panels, marked) when a run of empty bins
    separates a few far units from the rest; the y axis is broken when the zero bar dwarfs the rest.
    Draws into gridspec slot `slot`; returns the legend handles."""
    key = LOC if by == "locus" else TAX
    nmask = collections.Counter(key[MASKED]); units = sorted(set(key)); NTOT = len(units)
    k_of = np.array([nmask.get(u, 0) for u in units]); kmax = int(k_of.max()) if len(k_of) else 0
    cats = [(f"not masked · fails {j} of {NC}", GREYS[j], lambda idx, j=j: PLAIN[idx] & (NF[idx] == j)) for j in sorted(GREYS)] + \
           [(f"masked, too divergent · fails {j} of {NC}", REDS[j], lambda idx, j=j: DMK[idx] & (np.minimum(np.maximum(NF[idx], LOW), NC) == j)) for j in FAILS] + \
           [(f"masked, too similar · {j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}", BLUES[j], lambda idx, j=j: SMK[idx] & (np.minimum(NS[idx], 3) == j)) for j in (1, 2, 3)]
    pos = k_of[k_of > 0]
    if len(pos) >= 2:
        iqr = float(np.percentile(pos, 75) - np.percentile(pos, 25)); w = 2 * iqr / len(pos) ** (1 / 3)
        if w < 1: w = kmax / max(1.0, 1 + math.log2(len(pos)))          # Sturges when the quartiles coincide
        w = max(1, int(math.ceil(w)))
    else: w = 1
    bins = [(0, 0)] + [(lo_, min(lo_ + w - 1, kmax)) for lo_ in range(1, kmax + 1, w)]
    heights = np.array([int(((k_of >= lo_) & (k_of <= hi_)).sum()) for lo_, hi_ in bins])
    comp = np.zeros((len(bins), len(cats)))                            # bars with masks: the kinds of mask among their masked sequences; the zero bar: how many comparators its sequences fail
    for b, (lo_, hi_) in enumerate(bins):
        if heights[b] == 0: continue
        sel_units = [u for u, kk in zip(units, k_of) if lo_ <= kk <= hi_]
        idx = np.flatnonzero(np.isin(key, sel_units))
        shares = np.array([float(f(idx).sum()) for _, _, f in cats])
        if lo_ > 0: shares[:len(GREYS)] = 0
        shares = shares / max(shares.sum(), 1)
        comp[b] = heights[b] * shares
    centres = np.array([(lo_ + hi_) / 2 for lo_, hi_ in bins]); widths = np.array([(hi_ - lo_ + 1) * 0.9 for lo_, hi_ in bins])
    # x break: a run of at least three empty bins after which only a few occupied bins remain
    occ = [b for b in range(len(bins)) if heights[b] > 0]; xbreak = None
    for b1, b2 in zip(occ[:-1], occ[1:]):
        if b2 - b1 > 3 and sum(1 for b in occ if b >= b2) <= 3: xbreak = (b1, b2); break
    ybreak = len(bins) > 1 and heights[0] > 5 * max(int(heights[1:].max()), 1)
    ncol = 2 if xbreak else 1; nrow = 2 if ybreak else 1
    wr = (1,) if not xbreak else (xbreak[0] + 1, max(len(bins) - xbreak[1], 1)); wr = (wr[0], max(wr[1], 0.25 * wr[0])) if xbreak else wr
    sub = slot.subgridspec(nrow, ncol, height_ratios=(1, 2.6) if ybreak else (1,), width_ratios=wr, hspace=0.08, wspace=0.05)
    ax = np.empty((nrow, ncol), dtype=object)
    for r in range(nrow):
        for c in range(ncol):
            ax[r, c] = fig.add_subplot(sub[r, c], sharex=ax[0, c] if r else None, sharey=ax[r, 0] if c else None)
    cols = [c for _, c, _ in cats]
    for r in range(nrow):
        for c in range(ncol):
            a_x = ax[r, c]; base = np.zeros(len(bins))
            for ci_, col in enumerate(cols):
                v = comp[:, ci_]; nz = v > 0
                a_x.bar(centres[nz], v[nz], bottom=base[nz], color=col, edgecolor=SURF, linewidth=0.8, width=widths[nz]); base += v
            for s_ in ("top", "right"): a_x.spines[s_].set_visible(False)
            for s_ in ("left", "bottom"): a_x.spines[s_].set_color(SPINE)
            a_x.tick_params(colors=MUTED, length=0, labelsize=7); a_x.yaxis.grid(True, color=GRID, lw=0.5); a_x.set_axisbelow(True)
            a_x.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5, integer=True))
    edge = lambda b: (bins[b][0] - 0.5, bins[b][1] + 0.5)
    blabel = lambda b: str(bins[b][0]) if bins[b][0] == bins[b][1] else f"{bins[b][0]}–{bins[b][1]}"
    def label_bins(a_x, bs):                                          # every bar carries its range; thinned and turned when crowded
        bs = [b for b in bs if heights[b] > 0 or len(bs) <= 12]
        a_x.set_xticks([centres[b] for b in bs]); a_x.set_xticklabels([blabel(b) for b in bs], fontsize=6.5 if len(bs) <= 12 else 6, rotation=0 if len(bs) <= 12 else 90)
    if xbreak:
        ax[0, 0].set_xlim(-0.6, edge(xbreak[0])[1] + 0.6 * w); ax[0, 1].set_xlim(edge(xbreak[1])[0] - 0.6 * w, edge(len(bins) - 1)[1] + 0.6 * w)
        for r in range(nrow):
            ax[r, 1].spines["left"].set_visible(False); ax[r, 1].tick_params(labelleft=False); ax[r, 0].spines["right"].set_visible(False)
            d = 0.012
            for a_x, x0 in ((ax[r, 0], 1), (ax[r, 1], 0)):
                a_x.plot((x0 - d, x0 + d), (-2 * d, 2 * d), transform=a_x.transAxes, color=INK, clip_on=False, lw=0.9)
        label_bins(ax[nrow - 1, 0], range(0, xbreak[0] + 1)); label_bins(ax[nrow - 1, 1], range(xbreak[1], len(bins)))
    else:
        ax[0, 0].set_xlim(-0.6, edge(len(bins) - 1)[1] + 0.6 * w); label_bins(ax[nrow - 1, 0], range(len(bins)))
    if ybreak:
        top_lo = max(int(heights[1:].max()), 1) * 1.25; ax[1, 0].set_ylim(0, top_lo); ax[0, 0].set_ylim(heights[0] * 0.75, heights[0] * 1.12)
        for c in range(ncol):
            ax[0, c].spines["bottom"].set_visible(False); ax[1, c].spines["top"].set_visible(False); ax[0, c].tick_params(labelbottom=False)
            d = 0.012
            for a_x, y0 in ((ax[0, c], 0), (ax[1, c], 1)):
                a_x.plot((-d, +d), (y0 - 2 * d, y0 + 2 * d), transform=a_x.transAxes, color=INK, clip_on=False, lw=0.9)
        ax[0, 0].text(0, heights[0] * 1.02, f"{heights[0]:,}", ha="center", va="bottom", fontsize=7, color=INK)
        lab_from, lab_row = 1, 1
    else:
        ax[0, 0].set_ylim(0, max(int(heights.max()), 1) * 1.12); lab_from, lab_row = 0, 0
    for b in range(lab_from, len(bins)):
        if heights[b]:
            col_ = 1 if (xbreak and b >= xbreak[1]) else 0
            ax[lab_row, col_].text(centres[b], heights[b] + ax[lab_row, 0].get_ylim()[1] * 0.012, f"{heights[b]:,}", ha="center", va="bottom", fontsize=7, color=INK)
    unit = "loci" if by == "locus" else "taxa"
    ax[nrow - 1, 0].set_xlabel(f"sequences masked in the {by}", fontsize=8, color=MUTED); ax[nrow - 1, 0].set_ylabel(unit, fontsize=8, color=MUTED)
    if ybreak: ax[0, 0].set_ylabel(unit, fontsize=8, color=MUTED)
    ax[0, 0].set_title(f"by {by}  ·  {len(nmask)} of {NTOT} {unit} have a masked sequence", fontsize=8.5, color=INK, loc="left", pad=6)
    return [Patch(facecolor=c, edgecolor=SURF, label=lab) for lab, c, _ in cats]

def bars():
    """the two marginals side by side: loci and taxa by the number of masked sequences they hold"""
    fig = plt.figure(figsize=(8.5, 9.6), dpi=150)
    gs = fig.add_gridspec(2, 1, hspace=0.32, left=0.09, right=0.98, top=0.845, bottom=0.06)
    handles = hist_facet(fig, gs[0], "locus"); hist_facet(fig, gs[1], "taxon")
    fig.text(0.09, 0.975, f"{NAME}  ·  Number of sequences masked per locus and per taxon", fontsize=10, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.09, 0.95, f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)", fontsize=8.5, color=INK, ha="left", va="top")
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.09, 0.93), ncol=3, frameon=False, fontsize=6.5, labelcolor=INK, columnspacing=2.0, handletextpad=0.5)
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

def load_alignment():
    """the alignment as the tool encodes it: state index per site, -1 for anything outside the alphabet"""
    import locus_outliers as lo
    names, seqs, loci, fmt, declared, nexus = lo.read_input(a.alignment, a.partitions)
    seqs = [bytes(q) if isinstance(q, (bytes, bytearray)) else q.encode() for q in seqs]     # the readers return bytes
    dtype = a.datatype or declared
    if dtype is None:
        probe = b"".join(seqs)[:200000].upper().decode("ascii", "replace"); letters = [c for c in probe if c.isalpha()]
        dtype = "dna" if letters and sum(c in "ACGTUN" for c in letters) / len(letters) >= 0.95 else "aa"
    alphabet = lo.DNA if dtype == "dna" else lo.PROTEIN
    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = lut[np.frombuffer(b"".join(seqs), dtype=np.uint8).reshape(len(names), len(seqs[0]))]
    return names, X, dict(loci)

def best_matches(names, X, sites_of, taxon, min_sites=30):
    """for every tested cell of the taxon: the taxon it matches best at that locus, and the identity.
    Per-site divergence compared on the same sites (as locus_outliers): the closest taxon among those sharing the
    most sites with the focal is the champion; any taxon closer over its own shared sites challenges it on the sites
    all three have (at least min_sites) and wins only if it is closer there; ties go to the closest relative
    genome-wide. Identity is over the sites the pair shares."""
    ix = {t: i for i, t in enumerate(names)}; ti = ix[taxon]
    m = np.load(f"{a.dir}/m.npy"); mt = open(f"{a.dir}/taxa.txt").read().split(); mi = {t: i for i, t in enumerate(mt)}
    mrow = np.array([m[mi[t], mi[taxon]] if t in mi and t != taxon else np.nan for t in names])
    order = np.argsort(np.where(np.isnan(mrow), np.inf, mrow)); rank = np.empty(len(names), int); rank[order] = np.arange(len(names))
    close = [names[i] for i in order[:min(5, len(names) - 1)] if not np.isnan(mrow[i])]
    idx = np.flatnonzero(TAX == taxon); best = np.full(len(idx), "", dtype=object); ident = np.full(len(idx), np.nan)
    for k, i in enumerate(idx):
        sites = sites_of.get(LOC[i])
        if sites is None: continue
        Y = X[:, sites]; f = Y[ti]; fp = f >= 0
        both = (Y >= 0) & fp[None, :]; shared = both.sum(1); same = ((Y == f[None, :]) & both).sum(1)
        ok = shared >= min_sites; ok[ti] = False
        if not ok.any(): continue
        with np.errstate(invalid="ignore", divide="ignore"): dvg = np.where(ok, (shared - same) / np.maximum(shared, 1), np.inf)
        tier = ok & (shared >= shared[ok].max()); lo_ = dvg[tier].min()
        tied = np.flatnonzero(tier & (dvg == lo_)); c = int(tied[np.argmin(rank[tied])])
        for t in sorted(np.flatnonzero(dvg < lo_).tolist(), key=lambda t: (-shared[t], dvg[t], rank[t])):
            common = fp & (Y[c] >= 0) & (Y[t] >= 0)
            if int(common.sum()) < min_sites: continue
            dt_ = int(((Y[t] != f) & common).sum()); dc_ = int(((Y[c] != f) & common).sum())
            if dt_ < dc_ or (dt_ == dc_ and rank[t] < rank[c]): c = t
        best[k] = names[c]; ident[k] = 100.0 * (1 - dvg[c])
    return idx, best, ident, close

def matches():
    """one panel per taxon among the most masked: its masked sequences by the taxon they match best"""
    counts = collections.Counter(TAX[MASKED]); order = [t for t, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))][:a.match_taxa]
    if not order: return []
    names, X, sites_of = load_alignment()
    order = [t for t in order if t in names]
    MAXBARS = 20; ncol = 2 if len(order) > 1 else 1; nrow = math.ceil(len(order) / ncol); pw, ph = 6.4, 2.75
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 0.8, ph * nrow + 1.9), dpi=110, squeeze=False); H = ph * nrow + 1.9
    cats = [(REDS[j], lambda s, j=j: DMK[s] & (np.minimum(np.maximum(NF[s], LOW), NC) == j)) for j in FAILS] + \
           [(BLUES[j], lambda s, j=j: SMK[s] & (np.minimum(NS[s], 3) == j)) for j in (1, 2, 3)]
    for ax, t in zip(axs.ravel(), order):
        idx, best, ident, close = best_matches(names, X, sites_of, t)
        ok = best != ""; mk = MASKED[idx] & ok; um = ~MASKED[idx] & ok
        cm = collections.Counter(best[mk]); cu = collections.Counter(best[um]); nm, nu = int(mk.sum()), int(um.sum())
        taxa_sorted = sorted({u for u in cm} | {u for u in cu if cu[u] / max(nu, 1) * nm >= 1}, key=lambda u: (-cm[u], -cu[u]))   # a taxon the unmasked sequences predict stays on the axis even with no bar
        if len(taxa_sorted) > MAXBARS:
            rest = taxa_sorted[MAXBARS - 1:]; taxa_sorted = taxa_sorted[:MAXBARS - 1] + [f"other ({len(rest)} taxa)"]
            cm[taxa_sorted[-1]] = sum(cm[u] for u in rest); cu[taxa_sorted[-1]] = sum(cu[u] for u in rest)
            for u in rest: best[mk & (best == u)] = taxa_sorted[-1]
        xs = np.arange(len(taxa_sorted)); bottom = np.zeros(len(xs))
        for c, f in cats:
            h = np.array([float((mk & f(idx) & (best == u)).sum()) for u in taxa_sorted])
            if h.any(): ax.bar(xs, h, bottom=bottom, width=0.82, color=c, edgecolor=SURF, lw=0.4)
            bottom += h
        expd = np.array([cu[u] / max(nu, 1) * nm for u in taxa_sorted])
        ax.hlines(expd, xs - 0.41, xs + 0.41, color=INK, lw=1.1, zorder=4)
        top = max(bottom.max(), expd.max(), 1)
        for x_, o_ in zip(xs, bottom): ax.text(x_, o_ + top * 0.015, f"{int(o_)}", ha="center", va="bottom", fontsize=5.8, color=INK)
        ax.set_xticks(xs); ax.set_xticklabels(taxa_sorted, rotation=90, fontsize=5.8, color=INK)
        for lab, u in zip(ax.get_xticklabels(), taxa_sorted):
            if u in close: lab.set_fontweight("bold")
        ax.set_xlim(-0.7, max(len(xs) - 0.3, 0.7)); ax.set_ylim(0, top * 1.18)
        n_all = int(ok.sum()); pm = 100 * np.mean([u in close for u in best[mk]]) if nm else 0; pu = 100 * np.mean([u in close for u in best[um]]) if nu else 0
        title, fs = fit_title(t, pw - 0.5, 8)
        ax.set_title(title, fontsize=fs, color=INK, pad=22, linespacing=0.95)
        ax.text(0.5, 1.065, f"{nm} of {n_all} sequences masked ({100 * nm / max(n_all, 1):.1f}%)", transform=ax.transAxes, fontsize=6.3, color=INK, ha="center", va="bottom")
        ax.text(0.5, 1.012, f"best match is one of the 5 closest relatives for {pm:.0f}% of masked and {pu:.0f}% of unmasked sequences  ·  "
                            f"identity to it: masked {np.nanmedian(ident[mk]) if nm else float('nan'):.1f}%, unmasked {np.nanmedian(ident[um]) if nu else float('nan'):.1f}%",
                transform=ax.transAxes, fontsize=5.8, color=MUTED, ha="center", va="bottom")
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True, nbins=5))
        for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
        for s_ in ("left", "bottom"): ax.spines[s_].set_color(SPINE); ax.spines[s_].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, labelsize=6, length=0, pad=1.5); ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.set_ylabel("masked sequences", fontsize=6.5, color=MUTED)
    for ax in axs.ravel()[len(order):]: ax.axis("off")
    head = f"{NAME}  ·  Best match of each masked sequence"
    sub = f"{int(MASKED.sum()):,} of {len(R):,} sequences masked ({100 * MASKED.sum() / len(R):.1f}%)"
    note = (f"the {len(order)} most masked taxa; each masked sequence is matched against every other taxon at its locus (identity over the masked sequence's sites, a site the other taxon lacks counting as a mismatch; "
            "ties to the closest relative genome-wide); bars count masked sequences by the taxon they match best; the tick is the count expected if they matched like the taxon's unmasked sequences; bold: the 5 closest relatives genome-wide")
    fig.text(0.01, 1 - 0.14 / H, head, fontsize=11, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 1 - 0.36 / H, sub, fontsize=9, color=INK, ha="left", va="top")
    fig.text(0.01, 1 - 0.55 / H, note, fontsize=7.2, color=MUTED, ha="left", va="top", wrap=True)
    handles = [Patch(facecolor=c, edgecolor=SURF, label=f"masked, too divergent · fails {j} of {NC}") for j, c in REDS.items()] + \
              [Patch(facecolor=BLUES[j], edgecolor=SURF, label=f"masked, too similar · {j}{'+' if j == 3 else ''} pair{'s' if j > 1 else ''}") for j in (1, 2, 3)] + \
              [Line2D([], [], color=INK, lw=1.1, label="count expected if the masked sequences matched like the unmasked ones")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.78 / H), ncol=3, frameon=False, fontsize=7, handletextpad=0.5, columnspacing=2.0, labelcolor=INK)
    fig.tight_layout(rect=(0.012, 0.012, 1, 1 - 1.5 / H), h_pad=2.2, w_pad=1.0)
    return [fig]

makers = {"loci": lambda: [panels("locus")], "taxa": lambda: [panels("taxon")], "bars": bars, "matches": matches,
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
