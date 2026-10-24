#!/usr/bin/env python3
"""Is a masked sequence the reverse complement of its gene? Needs the nucleotide alignment that the protein
alignment was translated from, codon-matched (loci named by a shared ENSG id, as in Wu_2018).

    python strand.py DIR --aa AA.nex --dna DNA.nex [--example ENSG00000161544 --taxon MACAC_FAS]
                     [--windows DIR/badexon_TAXON.tsv --ref HOMO_SAP] [--out DIR/strand.tsv]

For every row of DIR/masked.csv: the ungapped DNA of the masked sequence against the ungapped DNA of up to
five of its closest relatives genome-wide (by m.npy) that have the locus, 12-mer sharing forward and
reverse-complement, the median over references (so one relative carrying the same error cannot decide):
reverse if the reverse-complement share beats the forward share by 0.05, forward if the opposite,
undecided otherwise. A control of 1,500 unmasked sequences is scored the same way. --example prints one
locus in full (the taxon's forward translation and the three reverse-complement frames next to the
reference protein). --windows runs the same test on the windows of a bad_exon.py table: do the windows
share 12-mers with the reference's coding sequence in either orientation?"""
import argparse, csv, re, sys, collections, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--aa", required=True); ap.add_argument("--dna", required=True)
ap.add_argument("--example"); ap.add_argument("--taxon", default="MACAC_FAS"); ap.add_argument("--windows"); ap.add_argument("--ref", default="HOMO_SAP"); ap.add_argument("--out")
a = ap.parse_args(); S = a.dir
def read_nexus_rows(path, want=None):
    names, seqs, charsets, inmat = [], {}, {}, False
    for line in open(path):
        s = line.strip()
        if not inmat:
            if s.lower() == "matrix": inmat = True; continue
            m = re.match(r"charset\s+(\S+)\s*=\s*(.+);", s, re.I)
            if m:
                r = re.match(r"^(\d+)-(\d+)$", m.group(2).strip()); charsets[m.group(1)] = (int(r.group(1)) - 1, int(r.group(2)))
            continue
        if s == ";": inmat = False; continue
        if not s: continue
        p = s.split(None, 1); names.append(p[0])
        if want is None or p[0] in want: seqs[p[0]] = p[1].replace(" ", "").upper()
    return names, seqs, charsets
names, dna, cs = read_nexus_rows(a.dna)
L = len(next(iter(dna.values()))); print(f"DNA alignment: {len(names)} taxa, {L:,} sites, {len(cs)} charsets", file=sys.stderr)
CODON = {}
bases = "TCAG"; aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a_ in enumerate(bases):
    for j, b_ in enumerate(bases):
        for k, c_ in enumerate(bases): CODON[a_ + b_ + c_] = aas[16 * i + 4 * j + k]
COMP = str.maketrans("ACGTN", "TGCAN")
def revcomp(s): return s.translate(COMP)[::-1]
def translate(s): return "".join(CODON.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def ungap(s): return re.sub(r"[^ACGT]", "", s)
# 1. one locus in full, when asked
def kshare(a_, b_, k=15):
    if len(a_) < k or len(b_) < k: return float("nan")
    kb = {b_[i:i + k] for i in range(len(b_) - k + 1)}
    return sum(1 for i in range(len(a_) - k + 1) if a_[i:i + k] in kb) / (len(a_) - k + 1)
if a.example:
    loc = [k for k in cs if a.example in k][0]; a0, a1 = cs[loc]
    h = ungap(dna[a.ref][a0:a1]); f = ungap(dna[a.taxon][a0:a1])
    print(f"{loc}: {a.ref} CDS {len(h)} nt, {a.taxon} {len(f)} nt")
    print(f"  {a.ref} protein, forward:      ", translate(h)[:60])
    print(f"  {a.taxon} forward frame 1:   ", translate(f)[:60])
    for fr in range(3): print(f"  {a.taxon} revcomp frame {fr + 1}:   ", translate(revcomp(f)[fr:])[:60])
    print(f"  15-mers of {a.taxon} found in {a.ref}: forward {100 * kshare(f, h):.0f}%, reverse complement {100 * kshare(revcomp(f), h):.0f}%")
# 2. every masked sequence in Wu: forward vs reverse-complement 15-mer sharing with its closest relative that has data
ens = {re.search(r"(ENSG\d+)", k).group(1): k for k in cs}
def dna_of(aa_locus): return cs[ens[re.search(r"(ENSG\d+)", aa_locus).group(1)]]
m = np.load(f"{S}/m.npy"); mt = open(f"{S}/taxa.txt").read().split(); mi = {t: i for i, t in enumerate(mt)}
orders = {t: np.argsort(np.where(np.isnan(m[mi[t]]), np.inf, m[mi[t]])) for t in mt}
def classify(l, t, seq=None, k=12, margin=0.05, nref=5):
    """forward or reverse complement, judged against the taxon's closest relatives that have data: the median k-mer
    sharing over up to nref references, so that one relative carrying the same error does not decide"""
    a0, a1 = dna_of(l); f = seq if seq is not None else ungap(dna[t][a0:a1])
    if len(f) < 4 * k: return None
    fw, rv, refs = [], [], []
    for j in orders[t]:
        if mt[j] == t: continue
        g = ungap(dna[mt[j]][a0:a1])
        if len(g) >= 60: fw.append(kshare(f, g, k)); rv.append(kshare(revcomp(f), g, k)); refs.append(mt[j])
        if len(refs) == nref: break
    if not refs: return None
    fwm, rvm = float(np.median(fw)), float(np.median(rv))
    return (",".join(refs), fwm, rvm, "reverse" if rvm > fwm + margin else "forward" if fwm > rvm + margin else "undecided")
masked = list(csv.DictReader(open(f"{S}/masked.csv")))
res = collections.defaultdict(lambda: collections.Counter()); rows = []
for r in masked:
    c = classify(r["locus"], r["taxon"])
    kind = c[3] if c else "too short"; res[r["taxon"]][kind] += 1
    if c: rows.append((r["locus"], r["taxon"], c[0], f"{c[1]:.2f}", f"{c[2]:.2f}", kind, r["side"]))
tot = collections.Counter()
print("\nmasked sequences by taxon: reverse complement / forward / undecided, judged by median 12-mer sharing with up to 5 closest relatives")
for t, v in sorted(res.items(), key=lambda kv: -kv[1]["reverse"]):
    tot.update(v)
    if v["reverse"] or sum(v.values()) >= 40: print(f"  {t:10s} reverse {v['reverse']:4d}  forward {v['forward']:4d}  undecided {v['undecided']:4d}  too short {v['too short']:3d}")
print(f"  {'all':10s} reverse {tot['reverse']:4d}  forward {tot['forward']:4d}  undecided {tot['undecided']:4d}  too short {tot['too short']:3d}")
with open(a.out or f"{S}/strand.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["locus", "taxon", "reference", "fwd_15mer_share", "rev_15mer_share", "kind", "side"]); w.writerows(rows)
# control: unmasked cells
cells = [r for r in csv.DictReader(open(f"{S}/cells.tsv"), delimiter="\t")]
mk = {(r["locus"], r["taxon"]) for r in masked}
rng = np.random.default_rng(0); um = [r for r in cells if (r["locus"], r["taxon"]) not in mk]; um = [um[i] for i in rng.choice(len(um), 1500, replace=False)]
ctrl = collections.Counter((classify(r["locus"], r["taxon"]) or [None] * 4)[3] or "too short" for r in um)
print("control, 1,500 unmasked sequences:", dict(ctrl))
# 3. the windows of a bad_exon.py table: on the reverse strand, forward, or foreign to the reference's CDS?
if a.windows:
    aa_cs = {}
    for line in open(a.aa):
        mm = re.match(r"\s*charset\s+(\S+)\s*=\s*(\d+)-(\d+);", line, re.I)
        if mm: aa_cs[mm.group(1)] = (int(mm.group(2)) - 1, int(mm.group(3)))
    bad = [r for r in csv.DictReader(open(a.windows), delimiter="\t") if r["masked"] == "True" and r["class"] == "windowed, clean elsewhere"]
    wtax = re.search(r"badexon_(.+)\.tsv$", a.windows).group(1) if re.search(r"badexon_(.+)\.tsv$", a.windows) else a.taxon
    wres = collections.Counter(); checked = 0; lens_ok = 0
    for r in bad:
        l = r["locus"]; b0, b1 = aa_cs[l]; d0, d1 = dna_of(l)
        if (d1 - d0) != 3 * (b1 - b0): continue
        lens_ok += 1
        for wnd in r["windows"].split(";"):
            s0, s1 = map(int, wnd.split("-")); segu = ungap(dna[wtax][d0 + 3 * s0: d0 + 3 * (s1 + 1)])
            if len(segu) < 48: continue
            checked += 1
            g = ungap(dna[a.ref][d0:d1]); fw, rv = kshare(segu, g, 12), kshare(revcomp(segu), g, 12)
            wres["reverse" if rv > fw + 0.1 else "forward" if fw > rv + 0.1 else "neither"] += 1
    print(f"\n{wtax} windows (loci with codon-matched DNA: {lens_ok} of {len(bad)}; windows of >= 48 nt checked: {checked}): 12-mers of the window found in the {a.ref} CDS, forward vs reverse complement:", dict(wres))
