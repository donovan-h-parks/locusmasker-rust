#!/usr/bin/env python3
"""Where along each masked sequence do its differences from its closest relatives sit?
    python bad_exon.py DIR ALIGNMENT TAXON [--out FILE]
A difference: a site where the focal taxon differs from every one of its 4 closest relatives that has the site.
Windows: runs of differences no more than GAP sites apart with at least MINW differences."""
import argparse, csv, re, sys, collections, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("aln"); ap.add_argument("taxon"); ap.add_argument("--out"); ap.add_argument("--gap", type=int, default=12); ap.add_argument("--minw", type=int, default=4)
a = ap.parse_args(); GAP, MINW = a.gap, a.minw
def read_nexus_simple(path):
    names, seqs, charsets, inmat = [], [], {}, False
    for line in open(path):
        s = line.strip()
        if not inmat:
            if s.lower() == "matrix": inmat = True; continue
            m = re.match(r"charset\s+(\S+)\s*=\s*(.+);", s, re.I)
            if m:
                sites = []
                for part in m.group(2).split():
                    r = re.match(r"^(\d+)-(\d+)$", part)
                    if r: sites.extend(range(int(r.group(1)) - 1, int(r.group(2))))
                charsets[m.group(1)] = np.array(sites, dtype=int)
            continue
        if s == ";": inmat = False; continue
        if not s: continue
        p = s.split(None, 1); names.append(p[0]); seqs.append(p[1].replace(" ", ""))
    return names, seqs, charsets
names, seqs, cs = read_nexus_simple(a.aln)
X = np.frombuffer("".join(s.upper() for s in seqs).encode(), dtype=np.uint8).reshape(len(names), -1)
miss = np.isin(X, np.frombuffer(b"-?XN*", dtype=np.uint8)); ix = {t: i for i, t in enumerate(names)}
m = np.load(f"{a.dir}/m.npy"); mt = open(f"{a.dir}/taxa.txt").read().split(); mi = {t: i for i, t in enumerate(mt)}
row = m[mi[a.taxon]]; rel = [mt[j] for j in np.argsort(np.where(np.isnan(row), np.inf, row))[:4]]
masked = [r for r in csv.DictReader(open(f"{a.dir}/masked.csv")) if r["taxon"] == a.taxon]
tested = [r["locus"] for r in csv.DictReader(open(f"{a.dir}/cells.tsv"), delimiter="\t") if r["taxon"] == a.taxon]
mset = {r["locus"] for r in masked}
def analyse(loc):
    s = cs[loc]; f = X[ix[a.taxon], s]; fp = ~miss[ix[a.taxon], s]
    R = np.stack([X[ix[t], s] for t in rel]); P = np.stack([~miss[ix[t], s] for t in rel])
    comparable = fp & (P.sum(0) >= 2)
    diff = comparable & ~((R == f[None, :]) & P).any(0)
    pos = np.flatnonzero(diff); ncomp = int(comparable.sum())
    # windows: runs of differences
    wins = []
    if len(pos):
        start = prev = pos[0]; cnt = 1
        for p_ in pos[1:]:
            if p_ - prev <= GAP: prev = p_; cnt += 1
            else:
                if cnt >= MINW: wins.append((start, prev, cnt))
                start = prev = p_; cnt = 1
        if cnt >= MINW: wins.append((start, prev, cnt))
    inwin = np.zeros(len(s), bool)
    for s0, s1, _ in wins: inwin[s0:s1 + 1] = True
    n_in = int((diff & inwin).sum()); n_out = int((diff & ~inwin).sum())
    comp_out = comparable & ~inwin; comp_in = comparable & inwin
    id_out = 100 * (1 - (diff & comp_out).sum() / comp_out.sum()) if comp_out.sum() else np.nan
    id_in = 100 * (1 - (diff & comp_in).sum() / comp_in.sum()) if comp_in.sum() else np.nan
    return dict(locus=loc, L=len(s), ncomp=ncomp, ndiff=len(pos), nwin=len(wins), winlen=int(inwin.sum()), n_in=n_in, n_out=n_out,
                id_out=id_out, id_in=id_in, pos=pos, wins=wins, fp=fp)
res = [analyse(l) for l in tested]; byloc = {r["locus"]: r for r in res}
M = [r for r in res if r["locus"] in mset]; U = [r for r in res if r["locus"] not in mset]
def classify(r):
    if r["ndiff"] < MINW: return "too few differences"
    frac_in = r["n_in"] / r["ndiff"]
    if r["nwin"] >= 1 and frac_in >= 0.8 and r["winlen"] <= 0.5 * r["ncomp"] and r["id_out"] >= 98: return "windowed, clean elsewhere"
    if r["nwin"] >= 1 and frac_in >= 0.8 and r["winlen"] <= 0.5 * r["ncomp"]: return "windowed, but not clean elsewhere"
    return "spread along the sequence"
cls = collections.Counter(classify(r) for r in M); clu = collections.Counter(classify(r) for r in U)
print(f"{a.taxon}: {len(M)} masked, {len(U)} unmasked; closest relatives used: {rel}")
print("masked:  ", dict(cls)); print("unmasked:", dict(clu))
W = [r for r in M if classify(r) == "windowed, clean elsewhere"]
print(f"windowed & clean: windows per sequence {collections.Counter(r['nwin'] for r in W)}; window length (sites) median {np.median([r['winlen'] for r in W]):.0f}, 10th {np.percentile([r['winlen'] for r in W], 10):.0f}, 90th {np.percentile([r['winlen'] for r in W], 90):.0f}; "
      f"identity inside windows median {np.nanmedian([r['id_in'] for r in W]):.0f}%, outside median {np.nanmedian([r['id_out'] for r in W]):.1f}%")
# unmasked cells that also carry a window: how many, and how big
UW = [r for r in U if classify(r) == "windowed, clean elsewhere"]
print(f"unmasked sequences with a window: {len(UW)}; their window length median {np.median([r['winlen'] for r in UW]) if UW else 0:.0f}, sequence length median {np.median([r['L'] for r in UW]) if UW else 0:.0f}")
with open(f"{a.dir}/badexon_{a.taxon}.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["locus", "masked", "class", "L", "ncomp", "ndiff", "nwin", "winlen", "id_in", "id_out", "windows"])
    for r in res: w.writerow([r["locus"], r["locus"] in mset, classify(r), r["L"], r["ncomp"], r["ndiff"], r["nwin"], r["winlen"], f"{r['id_in']:.1f}", f"{r['id_out']:.1f}", ";".join(f"{s0}-{s1}" for s0, s1, _ in r["wins"])])
# ------------------------------------------------------------------ barcode figure: one row per masked sequence
INK, MUTED, GRID, SPINE, SURF = "#0b0b0b", "#52514e", "#ebeae6", "#d4d2cc", "#ffffff"
RED, BLUE, GREY = "#a8201f", "#3987e5", "#c9c7c2"
order = {"windowed, clean elsewhere": 0, "windowed, but not clean elsewhere": 1, "spread along the sequence": 2, "too few differences": 3}
Ms = sorted(M, key=lambda r: (order[classify(r)], -r["winlen"] / max(r["ncomp"], 1) if order[classify(r)] < 2 else -r["ndiff"] / max(r["ncomp"], 1)))
ncol = 2; nrow = int(np.ceil(len(Ms) / ncol)); Lmax = max(r["L"] for r in Ms)
fig, axs = plt.subplots(1, ncol, figsize=(13, 0.052 * nrow + 1.6), dpi=150, squeeze=False); fig.patch.set_facecolor(SURF)
fig.subplots_adjust(left=0.03, right=0.99, top=1 - 1.05 / (0.052 * nrow + 1.6), bottom=0.5 / (0.052 * nrow + 1.6), wspace=0.06)
for c in range(ncol):
    ax = axs[0, c]; ax.set_facecolor(SURF)
    chunk = Ms[c * nrow:(c + 1) * nrow]
    for i, r in enumerate(chunk):
        y = nrow - i
        fp = np.flatnonzero(r["fp"])
        if len(fp): ax.hlines(y, fp.min(), fp.max() + 1, color=GREY, lw=1.6)   # the sequence, where it has data
        for s0, s1, _ in r["wins"]: ax.hlines(y, s0, s1 + 1, color=BLUE, lw=3.2, alpha=0.45)
        if len(r["pos"]): ax.vlines(r["pos"] + 0.5, y - 0.4, y + 0.4, color=RED, lw=0.5)
    ax.set_xlim(-5, Lmax + 5); ax.set_ylim(0.3, nrow + 0.7); ax.set_yticks([])
    for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
    ax.spines["bottom"].set_color(SPINE); ax.tick_params(colors=MUTED, labelsize=7, length=0); ax.set_xlabel("site within the locus (amino acids)", color=MUTED, fontsize=8)
    # class boundaries
    prev = None
    for i, r in enumerate(chunk):
        k = classify(r)
        if k != prev:
            ax.text(Lmax, nrow - i + 0.55, k, ha="right", va="bottom", fontsize=6.5, color=MUTED, style="italic"); prev = k
import os
ds = os.path.basename(os.path.normpath(a.dir))
fig.text(0.03, 0.985, f"{ds}  ·  {a.taxon}  ·  where each masked sequence differs from its closest relatives", fontsize=11, fontweight="bold", color=INK, ha="left", va="top")
H = 0.052 * nrow + 1.6
fig.text(0.03, 0.985 - 0.32 / H, f"{len(M)} masked sequences, one row each, drawn to length  ·  grey: the sequence where it has data  ·  red tick: a site where {a.taxon} differs from all of {', '.join(rel)} that have the site", fontsize=7.2, color=INK, ha="left", va="top")
fig.text(0.03, 0.985 - 0.52 / H, f"blue: a window of at least {MINW} such sites no more than {GAP} apart  ·  " + "  ·  ".join(f"{k}: {cls[k]}" for k in order if cls[k]) +
         f"  ·  windows: length median {np.median([r['winlen'] for r in W]):.0f} sites, identity inside {np.nanmedian([r['id_in'] for r in W]):.0f}%, outside {np.nanmedian([r['id_out'] for r in W]):.1f}%", fontsize=7.2, color=INK, ha="left", va="top")
out = a.out or f"{a.dir}/badexon_{a.taxon}.png"; fig.savefig(out, facecolor=SURF); print("wrote", out)
