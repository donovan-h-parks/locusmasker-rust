#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF.

    python plot_panels.py DIR [--plots log,linear,per_locus] [--out DIR/plots.pdf] [--png] [--top N] [--title NAME]

Pages (default all three, in this order; choose any subset with --plots):
  log       one panel per locus, both axes log scales on (100 - identity), labelled as identity
  linear    the same panels on plain axes
  per_locus a bar per locus of the cells flagged, split by how many comparators they failed

One point is one cell: the taxon's sequence at the locus against its closest comparator (x = that
pair's divergence shown as identity, 100 - m).
Greys: cells that pass (darker the more comparators they fail against). Reds: cells masked on the
divergent side (a majority up to all comparators failed, light to dark). Blues: cells masked on
the similar side (by number of too-similar pairs). Hollow red: ambiguous (fewer comparators than
requested, fails the ones it has). Hollow blue: in a similar-side event but not masked.
Solid line: the locus expectation, drawn over the range of the pairs it was fitted on. Dashed:
the effect floors (below: too divergent; above: too similar); the per-pair k*sqrt(E) arm is not
drawable. Axis ranges follow the dataset's cells and are shared by every panel. --png also
writes each page as a PNG next to the PDF."""
import argparse, csv, collections, math, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--plots", default="log,linear,per_locus")
ap.add_argument("--out"); ap.add_argument("--png", action="store_true"); ap.add_argument("--cols", type=int, default=7)
ap.add_argument("--highlight", type=int, default=40, help="loci with at least this many flags get a bold title")
ap.add_argument("--top", type=int, default=0, help="plot only the N loci with the most masked cells (0 = all)")
ap.add_argument("--title", help="dataset name for the page headers (default: the directory name)")
a = ap.parse_args(); out = a.out or f"{a.dir}/plots.pdf"
NAME = a.title or os.path.basename(os.path.normpath(a.dir))
wanted = [p.strip() for p in a.plots.split(",") if p.strip()]
bad = [p for p in wanted if p not in ("log", "linear", "per_locus")]
if bad: raise SystemExit(f"unknown plot(s): {bad}; choose from log, linear, per_locus")

# ---- data
per = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/per_locus.tsv"), delimiter="\t")}
cells = collections.defaultdict(list); trends = collections.defaultdict(list)
SIMMASK = set()
if os.path.exists(f"{a.dir}/masked.csv"):
    for r in csv.DictReader(open(f"{a.dir}/masked.csv")):
        if r.get("side") in ("similar", "both"): SIMMASK.add((r["locus"], r["taxon"]))
for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t"):
    x = 100.0 - float(r["pair_div"])
    cells[r["locus"]].append((x, float(r["identity"]), int(r["n_fail"]), int(r["flag"]), int(r["ambiguous"]),
                              int(r.get("n_similar", 0)), int((r["locus"], r["taxon"]) in SIMMASK)))
for r in csv.DictReader(open(f"{a.dir}/trends.tsv"), delimiter="\t"):
    trends[r["locus"]].append((100.0 - float(r["pair_div"]), float(r["trend"]), float(r["cutoff"]), float(r["cutoff_hi"])))
order = sorted((m for m in per if m in cells), key=lambda m: (-(int(per[m]["flagged"]) + int(per[m]["ambiguous"])), m))
if not order: raise SystemExit("nothing to plot: no locus was tested")
NTOT = len(order)
if a.top: order = order[:a.top]
FAILS = sorted(int(k[4:]) for k in next(iter(per.values())) if re.match(r"fail\d+$", k)) or [3, 4, 5]
NC, LOW = FAILS[-1], FAILS[0]            # comparators requested, and the fails that make a majority

# ---- tokens (light surface); text wears ink, marks wear the series colours
SURF, INK, MUTED, FAINT = "#ffffff", "#0b0b0b", "#52514e", "#8c8a85"
GRID, SPINE = "#ebeae6", "#d4d2cc"
MODEL, FLOOR = "#3b3a37", "#a5a39d"                     # expectation line; effect floors
def _ramp(c0, c1, n):
    u, v = np.array(matplotlib.colors.to_rgb(c0)), np.array(matplotlib.colors.to_rgb(c1))
    return [matplotlib.colors.to_hex(u + (v - u) * (i / max(n - 1, 1))) for i in range(n)]
GREYS = dict(zip(range(LOW), _ramp("#a8a6a1", "#4a4944", LOW)))          # passing cells, by comparators failed
REDS = dict(zip(FAILS, _ramp("#f29492", "#9c1c1c", len(FAILS))))         # divergent-side masks
BLUES = {1: "#86b6ef", 2: "#3987e5", 3: "#184f95"}                       # similar-side masks, by pairs
AMB, SIMRING = "#e34948", "#2a78d6"

def fit_title(name, width_in, fs=8.0):
    """a locus name that fits its panel: wrap once at a separator near the middle, then shrink"""
    cap = int(width_in * 72 / (0.56 * fs))
    if len(name) <= cap: return name, fs
    seps = [i for i, ch in enumerate(name) if ch in "_-. |" and 0.25 * len(name) < i < 0.75 * len(name)]
    i = min(seps, key=lambda i: abs(i - len(name) / 2)) + 1 if seps else len(name) // 2
    lines = [name[:i], name[i:]]
    longest = max(len(l) for l in lines)
    return "\n".join(lines), fs if longest <= cap else max(5.5, fs * cap / longest)

def log_ticks(lo, hi):
    """tick positions on a log divergence axis within [lo, hi], labelled as identity; the floor at
    lo carries the label 100, so no tick is placed within a quarter decade of it"""
    cand = [50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    t = [v for v in cand if lo * 10 ** 0.45 <= v <= hi]
    while len(t) > 5: t = t[::2]
    return t

def panels(LOG):
    ncol = min(a.cols, len(order)); nrow = math.ceil(len(order) / ncol)
    pw, ph = 2.3, 2.2
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 0.7, ph * nrow + 1.25), dpi=110, sharex=True, sharey=True, squeeze=False)
    allc = np.concatenate([np.array(cells[m]) for m in order])
    if LOG:
        dx, dy = 100.0 - allc[:, 0], 100.0 - allc[:, 1]
        floored_x, floored_y = bool((dx <= 0).any()), bool((dy <= 0).any())   # identity-100 cells sit on a floor tick
        fx = (dx[dx > 0].min() if (dx > 0).any() else 1.0) / (2.0 if floored_x else 1.0)
        Tx = lambda v: np.maximum(100.0 - np.asarray(v, float), fx)
        tmin = dy[dy > 0].min() if (dy > 0).any() else 1.0            # the y range must also hold every drawn expectation
        for m in order:
            t = np.array(trends[m]); g = 100.0 - t[:, 0]; xs_ = Tx(np.array(cells[m])[:, 0])
            sel = (g >= xs_.min() / 1.3) & (g <= xs_.max() * 1.3) & (100.0 - t[:, 1] > 0)
            if sel.any(): tmin = min(tmin, (100.0 - t[sel, 1]).min())
        fy = tmin / (2.0 if floored_y else 1.0)
        Ty = lambda v: np.maximum(100.0 - np.asarray(v, float), fy)
        xlim = (dx.max() * 1.4, fx / 1.4); ylim = (max(dy.max(), dx.max()) * 1.4, fy / 1.4)
    else:
        Tx = Ty = lambda v: np.asarray(v, float)
        xlim = (max(allc[:, 0].min() - 2, -1), 100.6); ylim = (max(min(allc[:, 1].min(), allc[:, 0].min()) - 3, -2), 101)
    for ax, m in zip(axs.ravel(), order):
        d = np.array(cells[m]); p = per[m]
        x, y, nf, fl, amb, sim, smk = d[:, 0], d[:, 1], d[:, 2].astype(int), d[:, 3] == 1, d[:, 4] == 1, d[:, 5] > 0, d[:, 6] == 1
        nsim = d[:, 5].astype(int)
        plain = ~fl & ~amb
        for k, c in GREYS.items():
            s = plain & (nf == k)
            if s.any(): ax.scatter(Tx(x[s]), Ty(y[s]), s=7, color=c, alpha=0.85, linewidths=0, rasterized=True, zorder=2)
        for k, c in REDS.items():
            s = fl & (nf == k)
            if s.any(): ax.scatter(Tx(x[s]), Ty(y[s]), s=16, color=c, edgecolors=SURF, linewidths=0.5, zorder=4)
        if amb.any(): ax.scatter(Tx(x[amb]), Ty(y[amb]), s=18, facecolors="none", edgecolors=AMB, linewidths=0.8, zorder=5)
        s = sim & ~smk
        if s.any(): ax.scatter(Tx(x[s]), Ty(y[s]), s=22, facecolors="none", edgecolors=SIMRING, linewidths=0.8, zorder=5)
        for k, c in BLUES.items():
            s = smk & (np.minimum(nsim, 3) == k)
            if s.any(): ax.scatter(Tx(x[s]), Ty(y[s]), s=16, color=c, edgecolors=SURF, linewidths=0.5, zorder=6)
        # the expectation and its floors, over the span of the pairs this locus was fitted on
        t = np.array(trends[m]); t = t[np.argsort(t[:, 0])]
        if LOG:
            g = np.log10(np.maximum(100.0 - t[:, 0], 1e-6)); lo, hi = np.log10(Tx(x).min()), np.log10(Tx(x).max())
            if hi - lo < math.log10(1.5): lo -= math.log10(1.25); hi += math.log10(1.25)
            xs = np.sort(np.concatenate([[lo], g[(g > lo) & (g < hi)], [hi]])); keep = np.argsort(g)
            col = lambda j: 10 ** np.interp(xs, g[keep], np.log10(Ty(t[keep, j])))
            px = 10 ** xs
        else:
            g = t[:, 0]; lo, hi = x.min(), x.max()
            if hi - lo < 1.0: lo -= 0.5; hi += 0.5
            xs = np.sort(np.concatenate([[lo], g[(g > lo) & (g < hi)], [hi]]))
            col = lambda j: np.interp(xs, g, t[:, j]); px = xs
        dash = col(2)
        clip = (lambda v: np.where(v > fy, v, np.nan)) if LOG else (lambda v: v)
        ax.plot(px, col(1), color=MODEL, lw=1.3, solid_capstyle="round", zorder=3)
        ax.plot(px, clip(dash), color=FLOOR, lw=0.9, ls=(0, (3, 2)), zorder=3)
        ax.plot(px, clip(col(3)), color=FLOOR, lw=0.9, ls=(0, (3, 2)), zorder=3)
        nfl = int(p["flagged"]); namb = int(p["ambiguous"]); hot = nfl + namb >= a.highlight
        title, fs = fit_title(m, pw - 0.25, 7.5)
        note = f"{nfl} of {p['tested']} flagged" + (f" · {namb} ambiguous" if namb else "")
        ax.set_title(title, fontsize=fs, color=INK, fontweight="bold" if hot else "normal", pad=12, linespacing=0.95)
        ax.text(0.5, 1.012, note, transform=ax.transAxes, fontsize=6.3, color=INK if nfl + namb else FAINT, ha="center", va="bottom")
        if LOG:
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*xlim); ax.set_ylim(*ylim)
            tx = log_ticks(fx, xlim[0]); ty = log_ticks(fy, ylim[0])
            ax.set_xticks(tx + ([fx] if floored_x else [])); ax.set_xticklabels([f"{100 - v:g}" for v in tx] + (["100"] if floored_x else []))
            ax.set_yticks(ty + ([fy] if floored_y else [])); ax.set_yticklabels([f"{100 - v:g}" for v in ty] + (["100"] if floored_y else []))
            for axis in (ax.xaxis, ax.yaxis): axis.set_minor_formatter(matplotlib.ticker.NullFormatter()); axis.set_minor_locator(matplotlib.ticker.NullLocator())
        else:
            ax.set_xlim(*xlim); ax.set_ylim(*ylim)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4)); ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        for s in ("left", "bottom"): ax.spines[s].set_color(SPINE); ax.spines[s].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, labelsize=6.5, length=0, pad=2); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for ax in axs.ravel()[len(order):]: ax.axis("off")
    xlab = "identity to closest comparator, 100 − m (%)"
    fig.supxlabel(xlab, fontsize=8, color=MUTED, y=0.004); fig.supylabel("locus identity (%)", fontsize=8, color=MUTED, x=0.004)
    tot = sum(int(p["flagged"]) for p in per.values()); tamb = sum(int(p["ambiguous"]) for p in per.values()); tested = sum(int(p["tested"]) for p in per.values())
    head = f"{NAME}  ·  {tot:,} of {tested:,} cells flagged" + (f", {tamb:,} ambiguous" if tamb else "") + (f"  ·  {len(order)} of {NTOT} loci shown" if len(order) < NTOT else "")
    sub = ("one point = a taxon's sequence at the locus vs its closest comparator  ·  solid = expected identity over the pairs it was fitted on  ·  "
           "dashed = effect floors (below: too divergent, above: too similar)" + ("  ·  both axes log on 100 − identity, labelled as identity" if LOG else ""))
    H = ph * nrow + 1.25
    fig.text(0.01, 1 - 0.14 / H, head, fontsize=11, color=INK, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 1 - 0.36 / H, sub, fontsize=7.5, color=MUTED, ha="left", va="top")
    handles = [Line2D([], [], marker="o", ls="", color=c, ms=4, label=("passes · " if k == 0 else "") + f"fails {k} of {NC}") for k, c in GREYS.items()] + \
              [Line2D([], [], marker="o", ls="", color=c, mec=SURF, ms=5.5, label=("masked, too divergent · " if k == LOW else "") + f"fails {k} of {NC}") for k, c in REDS.items()] + \
              [Line2D([], [], marker="o", ls="", mfc="none", mec=AMB, ms=6, label=f"ambiguous: < {NC} comparators, fails those it has")] + \
              ([Line2D([], [], marker="o", ls="", color=BLUES[k], mec=SURF, ms=5.5, label=("masked, too similar · " if k == 1 else "") + f"{k}{'+' if k == 3 else ''} pair{'s' if k > 1 else ''}") for k in (1, 2, 3)]
               + [Line2D([], [], marker="o", ls="", mfc="none", mec=SIMRING, ms=6, label="in a similar-side event, not masked")])
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.52 / H), ncol=min(len(handles), max(4, int((pw * ncol + 0.7) / 1.45))),
               frameon=False, fontsize=7, handletextpad=0.3, columnspacing=1.2, labelcolor=INK)
    fig.tight_layout(rect=(0.012, 0.012, 1, 1 - 0.98 / H), h_pad=1.0, w_pad=0.5)
    return fig

def per_locus():
    rows = [per[m] for m in order]; n = len(rows)
    lab_in = min(0.048 * max(len(r["locus"]) for r in rows) + 0.4, 4.5)
    fig, ax = plt.subplots(figsize=(5.6 + lab_in, 0.17 * n + 1.5), dpi=200)
    ys = np.arange(n)[::-1]; left = np.zeros(n)
    series = [(f"fail{k}", REDS[k], f"fails {k} of {NC}") for k in FAILS] + [("ambiguous", SURF, "ambiguous")]
    for key, c, lab in series:
        v = np.array([int(r[key]) for r in rows]); nz = v > 0                 # zero-width bars would still show their edge
        ax.barh(ys[nz], v[nz], left=left[nz], color=c, edgecolor=AMB if key == "ambiguous" else SURF, linewidth=0.8, height=0.72); left += v
    top = max(left.max(), 1)
    for y, tot in zip(ys, left):
        if tot: ax.text(tot + top * 0.012, y, f"{int(tot):,}", va="center", fontsize=6, color=MUTED)
    ax.set_yticks(ys); ax.set_yticklabels([r["locus"] for r in rows], fontsize=6)
    ax.set_xlabel("cells masked on the divergent side", fontsize=7.5, color=MUTED); ax.set_xlim(0, top * 1.12); ax.set_ylim(-0.7, n - 0.3)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_title(f"{NAME}  ·  {int(left.sum()):,} cells masked on the divergent side, by locus" + (f" (top {n} of {NTOT})" if n < NTOT else ""),
                 fontsize=9, color=INK, loc="left", fontweight="bold", pad=8)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=c, edgecolor=AMB if key == "ambiguous" else SURF, linewidth=0.8, label=lab) for key, c, lab in series],
              frameon=False, fontsize=6.5, loc="lower right", labelcolor=INK)
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(SPINE); ax.tick_params(colors=MUTED, length=0, labelsize=6.5); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    fig.tight_layout()
    return fig

makers = {"log": lambda: panels(True), "linear": lambda: panels(False), "per_locus": per_locus}
with PdfPages(out) as pdf:
    for name in wanted:
        fig = makers[name]()
        pdf.savefig(fig, bbox_inches="tight")
        if a.png: fig.savefig(f"{os.path.splitext(out)[0]}_{name}.png", bbox_inches="tight")
        plt.close(fig); print("page", name)
print("wrote", out)
