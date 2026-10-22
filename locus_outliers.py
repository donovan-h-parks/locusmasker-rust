#!/usr/bin/env python3
"""Per-locus outlier detection and masking in a multi-locus alignment.

    python locus_outliers.py --aln alignment.nex --out DIR
    python locus_outliers.py --aln alignment.fasta --partitions partitions.txt --out DIR
    python locus_outliers.py --aln alignment.phy   --partitions partitions.txt --out DIR

Input: a NEXUS alignment whose SETS block defines the loci, one charset per locus (read with
Biopython's Bio.Nexus), or a FASTA / PHYLIP / NEXUS alignment without charsets plus a RAxML-style
partition file ("DNA, locus1 = 1-100, 250-300" and "WAG, locus2 = 101-249\\3" forms; the model
prefix is optional). A NEXUS file with charsets and a partition file together is an error: give
one or the other. Sites are 1-based and inclusive. Data type is taken from the NEXUS file when
it states one, otherwise inferred from the residues; DNA is analysed on A, C, G, T only (U read
as T), protein on the 20 standard amino acids; everything else is missing.

Method (README.md): every taxon is compared, locus by locus, with its --n-comp closest
comparators (odd, at least 3, default 5; closest by alignment-wide identity). Each pair gives a
residual, the pair's locus identity minus the identity the locus trend predicts at the pair's
alignment-wide identity. A cell (taxon, locus) is an outlier when the median residual over its
comparators is below the band, i.e. it fails against a majority of them. Cells with fewer than
--n-comp comparators that fail against the ones they have are reported as ambiguous and masked
too. Masking blanks the cell's whole stretch of the locus, ambiguity codes included.

Outputs (in --out; <stem> is the alignment file's name without extension)
  <stem>.masked.<ext>  the alignment in its input format with every outlier cell replaced by gaps
  masked.csv           every masked cell with the evidence: locus, taxon, sites_masked, n_comp,
                       n_fail, median_residual, ambiguous, and for the closest comparator:
                       comparator, ani, overlap, identity, expected, residual
  cells.tsv            every tested cell (same columns plus flag, for unmasked cells too)
  pairs.tsv            every comparison: locus, taxon, rank, comparator, ani, overlap, identity,
                       expected, residual, fail
  per_locus.tsv        per locus: sites, taxa, tested, pairs, rate, med, spread, band, flagged,
                       fail<j> for j = majority .. n_comp, ambiguous
  trends.tsv           the locus trend on a grid, for plotting
  ani.npy, taxa.txt    taxon x taxon alignment-wide identity (%), and the row order
Options
  --n-comp 5 (odd, >= 3)  --fit-overlap 20  --min-band 25  --k 5  --bin 0.1  --bin-min 40  --seed 1
"""
import argparse, csv, os, re, sys, time
import numpy as np

PROTEIN = "ACDEFGHIKLMNPQRSTVWY"
DNA = "ACGT"


# ------------------------------------------------------------------ input
def read_fasta(txt):
    names, seqs, cur = [], [], []
    for line in txt.splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur))
            names.append(line[1:].split()[0]); cur = []
        else:
            cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    return names, seqs


def read_phylip(txt):
    lines = [l for l in txt.splitlines() if l.strip()]
    n, L = map(int, lines[0].split()); names, seqs = [], []
    for line in lines[1:]:
        p = line.split(None, 1)
        if len(p) == 2 and len(names) < n: names.append(p[0]); seqs.append(p[1].replace(" ", ""))
        elif names: seqs[-1] += line.replace(" ", "")     # interleaved continuation
    if any(len(s) != L for s in seqs):                      # interleaved with names only in the first block
        blocks = [l for l in lines[1:]]; names, seqs = [], []
        for i, line in enumerate(blocks):
            p = line.split(None, 1)
            if i < n: names.append(p[0]); seqs.append(p[1].replace(" ", "") if len(p) == 2 else "")
            else: seqs[i % n] += line.replace(" ", "")
    assert len(names) == n and all(len(s) == L for s in seqs), "cannot parse PHYLIP"
    return names, seqs


def parse_raxml_partitions(path, L):
    """name -> 0-based site index array; ranges 'a-b', 'a-b\\s' (step) and comma lists"""
    loci = []
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line or "=" not in line: continue
        lhs, rhs = line.split("=", 1)
        name = lhs.split(",")[-1].strip()
        sites = []
        for part in rhs.replace(";", "").split(","):
            part = part.strip()
            if not part: continue
            m = re.match(r"^(\d+)\s*-\s*(\d+)(?:\s*\\\s*(\d+))?$", part)
            if m:
                a, b, s = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
                sites.extend(range(a - 1, b, s))
            elif re.match(r"^\d+$", part):
                sites.append(int(part) - 1)
            else:
                sys.exit(f"cannot parse partition range {part!r} in {line!r}")
        arr = np.array(sorted(set(sites)), dtype=int)
        if arr.max() >= L: sys.exit(f"partition {name} runs past the alignment ({arr.max() + 1} > {L})")
        loci.append((name, arr))
    return loci


def infer_dtype(seqs):
    """DNA when at least 90% of the non-missing characters are A, C, G, T or U; else protein"""
    res = "".join(s.upper() for s in seqs[: min(len(seqs), 100)]); res = re.sub(r"[-?XN.~*]", "", res)
    if not res: return "protein"
    return "dna" if sum(res.count(c) for c in "ACGTU") / len(res) >= 0.9 else "protein"


def read_input(path, partitions):
    txt = open(path).read(); head = txt.lstrip()[:7].upper()
    nexus = None; dtype = None
    if head.startswith("#NEXUS"):
        from Bio.Nexus import Nexus
        nexus = Nexus.Nexus(path)
        names = list(nexus.taxlabels); seqs = [str(nexus.matrix[t]) for t in names]; fmt = "nexus"
        if nexus.charsets and partitions:
            sys.exit("the NEXUS file has charsets and a partition file was given: use one or the other, not both")
        if nexus.charsets:
            loci = [(k, np.array(sorted(v), dtype=int)) for k, v in nexus.charsets.items()]
        elif partitions:
            loci = parse_raxml_partitions(partitions, len(seqs[0]))
        else:
            sys.exit("the NEXUS file has no charsets: add a SETS block or give --partitions")
        dt = (nexus.datatype or "").lower()
        if dt in ("dna", "rna", "nucleotide"): dtype = "dna"
        elif dt == "protein": dtype = "protein"
    else:
        if not partitions: sys.exit("a FASTA or PHYLIP alignment needs --partitions")
        names, seqs = read_fasta(txt) if head.startswith(">") else read_phylip(txt)
        fmt = "fasta" if head.startswith(">") else "phylip"
        loci = parse_raxml_partitions(partitions, len(seqs[0]))
    L = len(seqs[0]); assert all(len(s) == L for s in seqs), "sequences differ in length"
    if not loci: sys.exit("no loci defined")
    if dtype is None: dtype = infer_dtype(seqs)
    return names, seqs, loci, fmt, dtype, nexus


# ------------------------------------------------------------------ output
def write_masked(path, out_dir, fmt, names, seqs, nexus):
    stem, ext = os.path.splitext(os.path.basename(path)); out = os.path.join(out_dir, f"{stem}.masked{ext or '.' + fmt}")
    if fmt == "nexus":
        from Bio.Seq import Seq
        for t, s in zip(names, seqs): nexus.matrix[t] = Seq(s)
        nexus.write_nexus_data(filename=out, interleave=False, append_sets=True)
    elif fmt == "fasta":
        with open(out, "w") as fh:
            for t, s in zip(names, seqs): fh.write(f">{t}\n{s}\n")
    else:
        with open(out, "w") as fh:
            fh.write(f"{len(names)} {len(seqs[0])}\n")
            for t, s in zip(names, seqs): fh.write(f"{t} {s}\n")
    return out


# ------------------------------------------------------------------ method
def trend(xf, yf, x, width, nmin):
    """running median of yf along xf evaluated at x: centred windows of `width`, widened x1.5
    until they hold nmin points; linear interpolation between centres, flat outside"""
    order = np.argsort(xf); xs, ys = xf[order], yf[order]
    lo, hi = xs[0], xs[-1]
    cx, cy = [], []
    for c in np.arange(lo + width / 2, hi, width):
        h = width / 2
        while True:
            sel = (xs >= c - h) & (xs <= c + h)
            if sel.sum() >= nmin or h > (hi - lo): break
            h *= 1.5
        if sel.sum() >= max(5, nmin // 4): cx.append(c); cy.append(np.median(ys[sel]))
    if not cx: return np.full(len(np.atleast_1d(x)), np.median(yf))
    return np.interp(x, np.array(cx), np.array(cy))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True); ap.add_argument("--partitions"); ap.add_argument("--out", required=True)
    ap.add_argument("--n-comp", type=int, default=5); ap.add_argument("--fit-overlap", type=int, default=20)
    ap.add_argument("--k", type=float, default=5.0); ap.add_argument("--min-band", type=float, default=25.0)
    ap.add_argument("--bin", type=float, default=0.1); ap.add_argument("--bin-min", type=int, default=40)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); K = a.n_comp
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    M = (K + 1) // 2                                   # fails that force the median below: a majority
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    names, seqs, loci, fmt, dtype, nexus = read_input(a.aln, a.partitions)
    if len(loci) < 3:
        sys.exit(f"only {len(loci)} locus defined: each locus is judged against the rest of the "
                 f"alignment, which needs at least 3 loci to mean anything")
    n, L = len(names), len(seqs[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    print(f"{n} taxa, {L} sites, {len(loci)} loci, {dtype}, {fmt}", file=sys.stderr)
    open(f"{a.out}/taxa.txt", "w").write("\n".join(names) + "\n")

    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = lut[np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L)]
    present = (X >= 0).astype(np.float32)
    shared = present @ present.T; same = np.zeros((n, n), dtype=np.float32)
    for i in range(len(alphabet)):
        A = (X == i).astype(np.float32); same += A @ A.T
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(shared > 0, 100.0 * same / shared, np.nan)
    np.fill_diagonal(ani, np.nan); np.save(f"{a.out}/ani.npy", ani)
    print(f"alignment-wide identities done, {time.time() - t0:.0f}s", file=sys.stderr)

    rng = np.random.default_rng(a.seed)
    ani_j = np.where(np.isnan(ani), -np.inf, ani + rng.random((n, n)) * 1e-6)
    lg = lambda v: np.log10(np.maximum(100.0 - v, 0.05))
    CELL = "locus\ttaxon\tn_comp\tcomparator\tani\tfocal_sites\toverlap\tidentity\texpected\tresidual\tmedian_residual\tn_fail\tflag\tambiguous\n"
    pairs = open(f"{a.out}/pairs.tsv", "w"); pairs.write("locus\ttaxon\trank\tcomparator\tani\toverlap\tidentity\texpected\tresidual\tfail\n")
    cells = open(f"{a.out}/cells.tsv", "w"); cells.write(CELL)
    per = open(f"{a.out}/per_locus.tsv", "w")
    per.write("locus\tsites\ttaxa\ttested\tpairs\trate\tmed\tspread\tband\tflagged\t"
              + "\t".join(f"fail{j}" for j in range(M, K + 1)) + "\tambiguous\n")
    trends = open(f"{a.out}/trends.tsv", "w"); trends.write("locus\tani\ttrend\n")
    mcsv = open(f"{a.out}/masked.csv", "w", newline=""); mw = csv.writer(mcsv)
    mw.writerow(["locus", "taxon", "sites_masked", "n_comp", "n_fail", "median_residual", "ambiguous", "comparator", "ani", "overlap", "identity", "expected", "residual"])
    grid = 100.0 - np.logspace(np.log10(60), np.log10(0.05), 200)
    masked = [bytearray(s.encode()) for s in seqs]; nmasked = 0
    for name, sites in loci:
        Xm = X[:, sites]; Pm = (Xm >= 0).astype(np.float32); nres = Pm.sum(1).astype(int)
        ov = Pm @ Pm.T
        need = np.minimum(a.fit_overlap, nres)[:, None]
        score = np.where(ov >= need, ani_j, -np.inf); np.fill_diagonal(score, -np.inf)
        kk_ = min(K, n - 1)
        top = np.argpartition(-score, kk_ - 1, axis=1)[:, :kk_]
        top = np.take_along_axis(top, np.argsort(-np.take_along_axis(score, top, 1), axis=1), 1)
        valid = np.isfinite(np.take_along_axis(score, top, 1)) & (nres > 0)[:, None]
        ncomp = valid.sum(1); idx = np.flatnonzero(ncomp > 0)
        if len(idx) < 2 * a.bin_min:
            print(f"{name}: only {len(idx)} testable cells, skipped", file=sys.stderr); continue
        gi = np.repeat(idx, kk_); ci = top[idx].ravel(); vi = valid[idx].ravel(); rk = np.tile(np.arange(kk_), len(idx))
        gi, ci, rk = gi[vi], ci[vi], rk[vi]
        both = (Xm[gi] >= 0) & (Xm[ci] >= 0); o = both.sum(1); kk = ((Xm[gi] == Xm[ci]) & both).sum(1)
        ident = 100.0 * kk / o; an = ani[gi, ci]; fit = o >= a.fit_overlap
        if fit.sum() < a.bin_min:
            print(f"{name}: only {int(fit.sum())} pairs with >= {a.fit_overlap} shared sites, skipped", file=sys.stderr); continue
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = float(np.median((100.0 - ident[fit]) / np.maximum(100.0 - an[fit], 0.1)))
        expected = trend(lg(an[fit]), ident[fit], lg(an), a.bin, a.bin_min); resid = ident - expected
        med = np.median(resid[fit]); spread = 1.4826 * np.median(np.abs(resid[fit] - med))
        band = max(a.min_band, a.k * spread); fail = resid < med - band
        for j in range(len(gi)):
            pairs.write(f"{name}\t{names[gi[j]]}\t{rk[j] + 1}\t{names[ci[j]]}\t{an[j]:.2f}\t{o[j]}\t{ident[j]:.1f}\t{expected[j]:.1f}\t{resid[j]:.1f}\t{int(fail[j])}\n")
        starts = np.r_[0, np.cumsum(ncomp[idx])[:-1]]
        fcount = {j: 0 for j in range(M, K + 1)}; namb = 0
        for g, s0, kc in zip(idx, starts, ncomp[idx]):
            r = resid[s0:s0 + kc]; mres = float(np.median(r)); nf = int(fail[s0:s0 + kc].sum())
            below = mres < med - band
            flag = int(below and kc >= K); amb = int(below and kc < K)
            if flag: fcount[nf] += 1
            namb += amb
            cells.write(f"{name}\t{names[g]}\t{kc}\t{names[ci[s0]]}\t{an[s0]:.2f}\t{nres[g]}\t{o[s0]}\t{ident[s0]:.1f}\t{expected[s0]:.1f}\t{resid[s0]:.1f}"
                        f"\t{mres:.1f}\t{nf}\t{flag}\t{amb}\n")
            if flag or amb:
                row = masked[g]; cnt = 0
                for site in sites:
                    if row[site] != ord("-"): cnt += 1
                    row[site] = ord("-")
                nmasked += 1
                mw.writerow([name, names[g], cnt, kc, nf, f"{mres:.1f}", amb, names[ci[s0]], f"{an[s0]:.2f}", o[s0], f"{ident[s0]:.1f}", f"{expected[s0]:.1f}", f"{resid[s0]:.1f}"])
        for gx, gy in zip(grid, trend(lg(an[fit]), ident[fit], lg(grid), a.bin, a.bin_min)): trends.write(f"{name}\t{gx:.2f}\t{gy:.2f}\n")
        nflag = sum(fcount.values()); fstr = "/".join(str(fcount[j]) for j in range(M, K + 1))
        per.write(f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t{len(idx)}\t{len(gi)}\t{rate:.3f}\t{med:.2f}\t{spread:.2f}\t{band:.1f}\t{nflag}\t"
                  + "\t".join(str(fcount[j]) for j in range(M, K + 1)) + f"\t{namb}\n")
        print(f"{name}: {len(idx)} cells, {len(gi)} pairs, band {band:.0f}, flagged {nflag} (fails {M}..{K} of {K}: {fstr}), ambiguous {namb}", file=sys.stderr)
    for fh in (pairs, cells, per, trends, mcsv): fh.close()
    out = write_masked(a.aln, a.out, fmt, names, [m.decode() for m in masked], nexus)
    print(f"masked {nmasked} cells -> {out}\ndone {time.time() - t0:.0f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
