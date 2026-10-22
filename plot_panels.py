#!/usr/bin/env python3
"""Figures from locus_outliers.py output, written to one PDF.

    python plot_panels.py DIR [--plots log,linear,per_locus] [--out DIR/plots.pdf] [--png]

Pages (default all three, in this order; choose any subset with --plots):
  log       one panel per locus, both axes log scales on (100 - identity) labelled as identity
  linear    the same panels on plain axes
  per_locus a bar per locus of the cells flagged, split by how many comparators they failed

In the panels one point is one cell: the taxon against its closest comparator. Greys by how
many of the cell's comparators it fails against (0, 1, 2: light to dark); reds for flagged cells
(3, 4, 5 of 5: light to dark), drawn on top; a hollow red circle is an ambiguous cell (fewer than
five comparators, fails against the ones it has). Solid line = locus trend, dashed = trend minus
the band. --png also writes each page as a PNG next to the PDF."""
import argparse, csv, collections, math, os, numpy as np, matplotlib
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
for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t"):
    cells[r["locus"]].append((float(r["ani"]), float(r["identity"]), int(r["n_fail"]), int(r["flag"]), int(r["ambiguous"])))
for r in csv.DictReader(open(f"{a.dir}/trends.tsv"), delimiter="\t"): trends[r["locus"]].append((float(r["ani"]), float(r["trend"])))
order = sorted(per, key=lambda m: (-(int(per[m]["flagged"]) + int(per[m]["ambiguous"])), m))
INK, MUTED, GRID, LINE = "#0b0b0b", "#52514e", "#e1e0d9", "#2a78d6"
GREYS = {0: "#c9c7c1", 1: "#8f8d87", 2: "#55534e"}; REDS = {3: "#f2a19f", 4: "#e35d5b", 5: "#b0171f"}; AMB = "#b0171f"
NC = max(int(r["n_comp"]) for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t")) if per else 5

def panels(LOG):
    T = (lambda v: np.maximum(100.0 - np.asarray(v, float), 0.1)) if LOG else (lambda v: np.asarray(v, float))
    ncol = min(a.cols, len(order)); nrow = math.ceil(len(order) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.2 * ncol + 0.5, 2.05 * nrow + 0.8), dpi=110, sharex=True, sharey=True, squeeze=False)
    for ax, m in zip(axs.ravel(), order):
        d = np.array(cells[m]); p = per[m]; band = float(p["band"])
        x, y, nf, fl, amb = d[:, 0], d[:, 1], d[:, 2].astype(int), d[:, 3] == 1, d[:, 4] == 1
        plain = ~fl & ~amb
        for k, c in GREYS.items():
            s = plain & (nf == k); ax.scatter(T(x[s]), T(y[s]), s=1.5, color=c, alpha=0.5 if k == 0 else 0.9, linewidths=0, rasterized=True, zorder=2)
        for k, c in REDS.items():
            s = fl & (nf == k); ax.scatter(T(x[s]), T(y[s]), s=5, color=c, alpha=0.95, linewidths=0, zorder=4)
        ax.scatter(T(x[amb]), T(y[amb]), s=9, facecolors="none", edgecolors=AMB, linewidths=0.7, zorder=5)
        t = np.array(trends[m]); t = t[t[:, 0] <= 99.6]
        ax.plot(T(t[:, 0]), T(t[:, 1]), color=LINE, lw=1.0, zorder=3); ax.plot(T(t[:, 0]), T(t[:, 1] - band), color=LINE, lw=0.9, ls=(0, (3, 2)), zorder=3)
        nfl = int(p["flagged"]); namb = int(p["ambiguous"]); hot = nfl + namb >= a.highlight
        ax.set_title(m, fontsize=8.5, color=REDS[5] if hot else INK, fontweight="bold" if hot else "normal", pad=3)
        txt = f"{nfl} flagged  (3/4/5 of {NC}: {p['fail3']}/{p['fail4']}/{p['fail5']})" + (f", {namb} ambiguous" if namb else "") + f"\nband {band:.0f} pts · n {p['tested']}"
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
    for ax in axs[-1]: ax.set_xlabel("alignment-wide identity to closest comparator (%)", fontsize=7.5, color=MUTED)
    for ax in axs[:, 0]: ax.set_ylabel("locus identity (%)", fontsize=7.5, color=MUTED)
    tot = sum(int(p["flagged"]) for p in per.values()); tamb = sum(int(p["ambiguous"]) for p in per.values()); tested = sum(int(p["tested"]) for p in per.values())
    fig.suptitle(f"{tot:,} of {tested:,} cells flagged" + (f", {tamb} ambiguous" if tamb else "") + " · one point = taxon vs its closest comparator · solid = locus trend, dashed = cutoff"
                 + (" · log axes on 100 − identity, labelled as identity" if LOG else ""), fontsize=11, color=INK, x=0.01, ha="left", y=0.996)
    handles = [Line2D([], [], marker="o", ls="", color=GREYS[k], ms=4, label=f"fails {k} of {NC}") for k in (0, 1, 2)] + \
              [Line2D([], [], marker="o", ls="", color=REDS[k], ms=5, label=f"flagged, fails {k} of {NC}") for k in (3, 4, 5)] + \
              [Line2D([], [], marker="o", ls="", mfc="none", mec=AMB, ms=6, label=f"ambiguous: < {NC} comparators, fails those it has")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 0.988), ncol=7, frameon=False, fontsize=7.5, handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout(rect=(0, 0, 1, 0.968), h_pad=0.9, w_pad=0.6)
    return fig

def per_locus():
    rows = [per[m] for m in order]; n = len(rows)
    fig, ax = plt.subplots(figsize=(6.4, 0.16 * n + 1.4), dpi=200)
    ys = np.arange(n)[::-1]; left = np.zeros(n)
    for key, c, lab in (("fail3", REDS[3], f"fails 3 of {NC}"), ("fail4", REDS[4], f"fails 4 of {NC}"), ("fail5", REDS[5], f"fails 5 of {NC}"), ("ambiguous", "#ffffff", "ambiguous")):
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
