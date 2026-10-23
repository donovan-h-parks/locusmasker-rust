#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF.

    python plot_panels.py DIR [--plots log,linear,per_locus] [--out DIR/plots.pdf] [--png]

Pages (default all three, in this order; choose any subset with --plots):
  log       one panel per locus, both axes log scales on (100 - identity) labelled as identity
  linear    the same panels on plain axes
  per_locus a bar per locus of the cells flagged, split by how many comparators they failed

In the panels one point is one cell: the taxon against its closest comparator (x = that pair's
identity: alignment-wide in band mode, 100 - m in z mode). Greys by how many of the cell's
comparators it fails against on the divergent side (below a majority: light to dark); reds for
flagged cells (a majority up to all: light to dark), drawn on top; a hollow red circle is an
ambiguous cell (fewer comparators than requested, fails against the ones it has); a hollow blue
circle is a cell involved in a similar-side event. Solid line = locus expectation; dashed = the
cutoff (band mode: trend + locus median residual - band; z mode: the min-ratio effect floor
only — the per-pair k*sqrt(E) arm is not drawable). --png also writes each page as a PNG."""
import argparse, csv, collections, math, os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--plots", default="log,linear,per_locus")
ap.add_argument("--out"); ap.add_argument("--png", action="store_true"); ap.add_argument("--cols", type=int, default=7)
ap.add_argument("--highlight", type=int, default=40, help="loci with at least this many flags get a red title")
a = ap.parse_args(); out = a.out or f"{a.dir}/plots.pdf"
wanted = [p.strip() for p in a.plots.split(",") if p.strip()]
bad = [p for p in wanted if p not in ("log", "linear", "per_locus")]
if bad: raise SystemExit(f"unknown plot(s): {bad}; choose from log, linear, per_locus")
per = {r["locus"]: r for r in csv.DictReader(open(f"{a.dir}/per_locus.tsv"), delimiter="\t")}
cells = collections.defaultdict(list); trends = collections.defaultdict(list)
ZMODE = "pair_div" in open(f"{a.dir}/cells.tsv").readline()
for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t"):
    x = float(r["ani"]) if "ani" in r else 100.0 - float(r["pair_div"])
    cells[r["locus"]].append((x, float(r["identity"]), int(r["n_fail"]), int(r["flag"]), int(r["ambiguous"]),
                              int(r.get("n_similar", 0))))
for r in csv.DictReader(open(f"{a.dir}/trends.tsv"), delimiter="\t"):
    x = float(r["ani"]) if "ani" in r else 100.0 - float(r["pair_div"])
    trends[r["locus"]].append((x, float(r["trend"]), float(r["cutoff"]) if "cutoff" in r else float("nan")))
order = sorted((m for m in per if m in cells), key=lambda m: (-(int(per[m]["flagged"]) + int(per[m]["ambiguous"])), m))
if not order: raise SystemExit("nothing to plot: no locus was tested")
INK, MUTED, GRID, LINE, AMB = "#0b0b0b", "#52514e", "#e1e0d9", "#2a78d6", "#b0171f"
FAILS = sorted(int(k[4:]) for k in (next(iter(per.values())) if per else {}) if re.match(r"fail\d+$", k)) or [3, 4, 5]
NC, LOW = FAILS[-1], FAILS[0]            # comparators requested, and the fails that make a majority
def _ramp(c0, c1, n):
    u, v = np.array(matplotlib.colors.to_rgb(c0)), np.array(matplotlib.colors.to_rgb(c1))
    return [matplotlib.colors.to_hex(u + (v - u) * (i / max(n - 1, 1))) for i in range(n)]
GREYS = dict(zip(range(LOW), _ramp("#c9c7c1", "#55534e", LOW)))
REDS = dict(zip(FAILS, _ramp("#f2a19f", "#b0171f", len(FAILS))))

def panels(LOG):
    T = (lambda v: np.maximum(100.0 - np.asarray(v, float), 0.1)) if LOG else (lambda v: np.asarray(v, float))
    ncol = min(a.cols, len(order)); nrow = math.ceil(len(order) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.2 * ncol + 0.5, 2.05 * nrow + 0.8), dpi=110, sharex=True, sharey=True, squeeze=False)
    for ax, m in zip(axs.ravel(), order):
        d = np.array(cells[m]); p = per[m]
        x, y, nf, fl, amb, sim = d[:, 0], d[:, 1], d[:, 2].astype(int), d[:, 3] == 1, d[:, 4] == 1, d[:, 5] > 0
        plain = ~fl & ~amb
        for k, c in GREYS.items():
            s = plain & (nf == k); ax.scatter(T(x[s]), T(y[s]), s=1.5, color=c, alpha=0.5 if k == 0 else 0.9, linewidths=0, rasterized=True, zorder=2)
        for k, c in REDS.items():
            s = fl & (nf == k); ax.scatter(T(x[s]), T(y[s]), s=5, color=c, alpha=0.95, linewidths=0, zorder=4)
        ax.scatter(T(x[amb]), T(y[amb]), s=9, facecolors="none", edgecolors=AMB, linewidths=0.7, zorder=5)
        ax.scatter(T(x[sim]), T(y[sim]), s=12, facecolors="none", edgecolors=LINE, linewidths=0.7, zorder=5)
        t = np.array(trends[m]); t = t[t[:, 0] <= 99.6]
        dash = t[:, 2] if np.isfinite(t[0, 2]) else t[:, 1] + float(p["med"]) - float(p["band"])
        ax.plot(T(t[:, 0]), T(t[:, 1]), color=LINE, lw=1.0, zorder=3); ax.plot(T(t[:, 0]), T(dash), color=LINE, lw=0.9, ls=(0, (3, 2)), zorder=3)
        nfl = int(p["flagged"]); namb = int(p["ambiguous"]); hot = nfl + namb >= a.highlight
        ax.set_title(m, fontsize=8.5, color=REDS[NC] if hot else INK, fontweight="bold" if hot else "normal", pad=3)
        scale = f"band {float(p['band']):.0f} pts" if "band" in p else f"cut z {float(p['cut']):.1f}"
        txt = f"{nfl} flagged  ({LOW}..{NC} of {NC}: " + "/".join(p[f"fail{j}"] for j in FAILS) + ")" + (f", {namb} ambiguous" if namb else "") + f"\n{scale} · n {p['tested']}"
        ax.text(0.03, 0.90 if LOG else 0.04, txt, transform=ax.transAxes, fontsize=6.0, color=INK if nfl + namb else MUTED, va="top" if LOG else "bottom")
        if LOG:
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(45, 0.08); ax.set_ylim(110, 0.08)
            tx = [40, 20, 10, 5, 2, 1, 0.1]; ax.set_xticks(tx); ax.set_xticklabels([f"{100 - v:g}" for v in tx])
            ty = [100, 50, 20, 10, 5, 2, 1, 0.1]; ax.set_yticks(ty); ax.set_yticklabels([f"{100 - v:g}" for v in ty])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.tick_params(which="minor", length=0)
        else:
            ax.set_xlim(55, 100.5); ax.set_ylim(-2, 101); ax.set_xticks([60, 80, 100]); ax.set_yticks([0, 50, 100])
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        for s in ("left", "bottom"): ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=7, length=0); ax.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
    for ax in axs.ravel()[len(order):]: ax.axis("off")
    xlab = "identity to closest comparator, 100 − m (%)" if ZMODE else "alignment-wide identity to closest comparator (%)"
    for ax in axs[-1]: ax.set_xlabel(xlab, fontsize=7.5, color=MUTED)
    for ax in axs[:, 0]: ax.set_ylabel("locus identity (%)", fontsize=7.5, color=MUTED)
    tot = sum(int(p["flagged"]) for p in per.values()); tamb = sum(int(p["ambiguous"]) for p in per.values()); tested = sum(int(p["tested"]) for p in per.values())
    fig.suptitle(f"{tot:,} of {tested:,} cells flagged" + (f", {tamb} ambiguous" if tamb else "") + " · one point = taxon vs its closest comparator · solid = locus trend, dashed = cutoff"
                 + (" · log axes on 100 − identity, labelled as identity" if LOG else ""), fontsize=11, color=INK, x=0.01, ha="left", y=0.996)
    handles = [Line2D([], [], marker="o", ls="", color=c, ms=4, label=f"fails {k} of {NC}") for k, c in GREYS.items()] + \
              [Line2D([], [], marker="o", ls="", color=c, ms=5, label=f"flagged, fails {k} of {NC}") for k, c in REDS.items()] + \
              [Line2D([], [], marker="o", ls="", mfc="none", mec=AMB, ms=6, label=f"ambiguous: < {NC} comparators, fails those it has")] + \
              ([Line2D([], [], marker="o", ls="", mfc="none", mec=LINE, ms=6, label="in a similar-side event")] if ZMODE else [])
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 0.988), ncol=7, frameon=False, fontsize=7.5, handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout(rect=(0, 0, 1, 0.968), h_pad=0.9, w_pad=0.6)
    return fig

def per_locus():
    rows = [per[m] for m in order]; n = len(rows)
    fig, ax = plt.subplots(figsize=(6.4, 0.16 * n + 1.4), dpi=200)
    ys = np.arange(n)[::-1]; left = np.zeros(n)
    for key, c, lab in [(f"fail{k}", REDS[k], f"fails {k} of {NC}") for k in FAILS] + [("ambiguous", "#ffffff", "ambiguous")]:
        v = np.array([int(r[key]) for r in rows]); ax.barh(ys, v, left=left, color=c, edgecolor=AMB if key == "ambiguous" else "none", linewidth=0.6, height=0.7, label=lab); left += v
    top = max(left.max(), 1)
    for y, tot in zip(ys, left):
        if tot: ax.text(tot + top * 0.01, y, f"{int(tot)}", va="center", fontsize=6, color=MUTED)
    ax.set_yticks(ys); ax.set_yticklabels([r["locus"] for r in rows], fontsize=5.8)
    ax.set_xlabel("cells flagged", fontsize=8, color=MUTED); ax.set_xlim(0, top * 1.1)
    ax.set_title(f"{int(left.sum()):,} flagged cells by locus", fontsize=9.5, color=INK, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(colors=MUTED, length=0); ax.xaxis.grid(True, color=GRID, lw=0.5); ax.set_axisbelow(True)
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
