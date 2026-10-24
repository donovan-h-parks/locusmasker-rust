#!/usr/bin/env python3
"""Per-locus outlier detection and masking in a multi-locus alignment.

    python locus_outliers.py --aln alignment.nex --datatype dna --out DIR
    python locus_outliers.py --aln alignment.fasta --partitions partitions.txt --datatype aa --out DIR
    python locus_outliers.py --aln alignment.phy   --partitions partitions.txt --datatype dna --out DIR

Input: a NEXUS alignment whose SETS block defines the loci, one charset per locus (read with
Biopython's Bio.Nexus, so only what Bio.Nexus accepts: whitespace-separated site lists, quoted
names with punctuation), or a FASTA / PHYLIP / NEXUS alignment without charsets plus a
RAxML-style partition file ("DNA, locus1 = 1-100, 250-300" and "WAG, locus2 = 101-249\\3"
forms; the model prefix is optional). Sites are 1-based and inclusive; a site may belong to one
locus only. --datatype is required: dna is analysed on A, C, G, T (U read as T), aa on the 20
standard amino acids; everything else is missing. This is a phylogenomic test, and one number — 10 — sets its
scope: the alignment must have at least 10 loci, and a taxon is only considered (tested, or
used as a comparator) when it has at least 10 loci with >= --fit-overlap residues; taxa below
that are left untouched in the output and listed in excluded_taxa.txt. At least 5 taxa.

Method:
  m         every taxon pair gets a divergence m: the site-weighted median across loci of the
            pair's rate-normalised divergence, each locus weighted by the sites the pair shares
            in it. Pairs overlapping in fewer than 5 loci have no m and are never compared
            (there is no pooled-identity fallback: it self-contaminates).
  floor     a pair's divergence is known only to the resolution of its loci: a locus with
            zero differences over n shared sites cannot exclude a divergence of about 1/n.
            Each pair therefore also gets a resolution floor m_res, the same site-weighted
            median over its loci of half a mismatch per locus (rate-normalised). A pair whose
            m is at or below its floor is "unresolved": it never votes on a locus rate, and
            on the too-divergent side its expectation is evaluated at m_res rather than m (the
            too-similar side uses m: a pair cannot be more similar than its loci allow). For
            a pair with abundant sequence the floor is negligible and nothing changes.
  rates     a locus rate is the median, over a set of taxon pairs, one equal vote per pair,
            of c(d + half a mismatch) / c(m), where d is the pair's mismatch fraction at the
            locus, c(p) = -P ln(1 - p/P) the multiple-hit correction with P the locus's
            random-identity plateau (from its residue frequencies), and m the pair's
            divergence (or its floor); rates are rescaled so the median locus has rate 1.
            Rates are estimated twice. Before m exists, each taxon's --n-comp nearest
            neighbours by the raw site-weighted median of locus divergences vote; these
            provisional rates normalise the loci when m is computed. In the test, each
            taxon's --n-comp closest comparators by m vote; that rate gives the expectation.
            Nearest-neighbour pairs are used because the test only ever compares close
            relatives, and a median over all pairs is a median over the deep branches of
            the tree, where a locus's rate can differ several-fold from its rate near the
            tips.
  compare   each taxon is compared, locus by locus, with its --n-comp closest taxa by m (odd,
            auto-shrunk to the largest odd <= n-1) among those overlapping it by >=
            min(--fit-overlap, its residues). Ties are broken by pooled identity, then taxon
            order, so runs are reproducible.
  expect    the expected corrected locus distance of a pair is r x c(max(m, m_res)), r the
            locus's test-stage rate; back-transformed through the saturation curve to an
            expected mismatch fraction for the count test. A straight line through the
            origin in corrected-distance space, one parameter per locus, 50% breakdown over
            the comparator pairs, nothing to overfit.
  dispersion the scatter of locus divergence between close relatives is wider than counting
            noise, and differently so for different loci. Each locus gets a dispersion phi:
            the robust spread (1.4826 MAD) of log(observed/expected) over its resolved
            comparator pairs with at least 10 expected mismatches, squared, minus the
            counting part 1/E, floored at 0; it is shrunk toward the dataset-wide value
            (the pair-weighted median over loci) with the weight of 30 pairs, so small
            datasets use one dispersion and large ones get per-locus values.
  test      exact tail probabilities, no normal approximation. With E the expected mismatches
            among the o shared sites, the observed count is binomial(o, E/o), or
            beta-binomial with the same mean and variance E + phi E^2 when the locus is
            dispersed. A comparison fails the DIVERGENT side when observed > --min-ratio x E
            and P(X >= observed) < alpha / 2N; the SIMILAR side when observed < E/3 and
            P(X <= observed) < alpha / 2N. N is the number of distinct comparisons in the run
            (both sides): Bonferroni, so the bar scales with dataset size. The effect floors
            keep statistically significant but biologically trivial excesses in very long
            loci from being masked; the tail test is the decision.
  decide    a cell (taxon x locus) is DIVERGENT-masked when a majority of its comparators fail
            the divergent side (fewer comparators than requested failing a majority = masked
            as ambiguous). Any single similar-side pair flags a SIMILAR event; its direction
            is attributed by a likelihood ratio over the sites the pair shares with each of
            up to 10 third parties (which party's expectations fit the shared sequence):
            |2 ln LR| > 10 (Kass & Raftery "very strong") masks the recipient, otherwise both
            cells. --mask-similar-only / --mask-different-only restrict masking to one side;
            similar.tsv is written regardless. Masking blanks the cell's whole stretch of the
            locus.

Outputs (in --out; <stem> is the alignment file's name without extension)
  <stem>.masked.<ext>  the alignment in its input format with every masked cell blanked
  masked.csv           ONE row per masked cell: side (divergent / similar / both), evidence
                       (divergent rows: the median over the cell's comparators of identity,
                       expected identity and excess mismatches, plus the closest comparator;
                       similar rows: the partner and the pair's own values), lambda
  cells.tsv            every tested cell: closest-comparator columns plus the medians over
                       comparators (med_identity, med_expected, med_excess), n_fail, n_similar, and
                       excess_m: the cell's extra divergence beyond expectation as corrected distance
                       in units of m (a typical locus's divergence) — a foreign or paralogous sequence
                       adds a constant here whatever its comparators, so cells with one cause share it
  pairs.tsv            every comparison
  similar.tsv          every similar-side event: pair, deficit, lambda, donor/recipient call
  per_locus.tsv        one row per locus, skipped loci included (status column)
  trends.tsv           expectation and both effect floors (too divergent / too similar) on a grid
  dispersion.tsv       per locus: raw dispersion estimate, pairs it used, dispersion applied
  ani.npy, m.npy       alignment-wide identity (%) and pair divergence m (points); m_hi.npy is
                       max(m, m_res), the divergence the too-divergent test is evaluated at
                       (pair_div in pairs.tsv, cells.tsv and masked.csv reports this value);
                       taxa.txt
Options
  --datatype dna|aa (required)  --n-comp 5 (odd >= 3)  --fit-overlap 20
  --alpha 0.05  --min-ratio 2  --mask-similar-only  --mask-different-only
  (constants in the code: HALF = 0.5 mismatch for the resolution floor and the rate
  pseudocount; dispersion from pairs with >= 10 expected mismatches, shrinkage weight 30)
"""
import argparse, collections, csv, math, os, re, sys, time
import numpy as np

PROTEIN = "ACDEFGHIKLMNPQRSTVWY"
DNA = "ACGT"
HALF = 0.5              # mismatches imputed to a locus that shows none (half the detection limit);
                        # also the pseudocount on the locus side of the rate ratio
PHI_MIN_E = 10.0        # a pair contributes to the locus dispersion only if it expects this many mismatches
PHI_N0 = 30.0           # shrinkage weight: a locus's dispersion is pulled toward the dataset's with this many pairs' worth
MIN_LOCI_FOR_M = 5      # pairs overlapping in fewer loci have no m and are never compared
MIN_LOCI_PER_TAXON = 10 # a taxon needs this many loci with >= fit-overlap residues to take part
MIN_FIT_PAIRS = 8       # a locus needs this many fit pairs to be tested (or to get a rate)
LR_BAR = 10.0           # |2 ln LR| above which a similar event's direction is attributed


# ------------------------------------------------------------------ input
def read_fasta(txt):
    names, seqs, cur = [], [], []
    for line in txt.splitlines():
        if line.startswith(">"):
            if cur: seqs.append("".join(cur))
            name = line[1:].split()
            if not name: sys.exit("FASTA: a header line ('>') has no name")
            names.append(name[0]); cur = []
        else:
            cur.append("".join(line.split()))
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
    if len(names) != n or any(len(s) != L for s in seqs):
        sys.exit("cannot parse PHYLIP: relaxed one-line or interleaved layouts are supported; "
                 "try FASTA plus a RAxML-style partition file")
    return names, seqs


def parse_raxml_partitions(path, L):
    """name -> 0-based site index array; ranges 'a-b', 'a-b\\s' (step) and comma lists;
    every site may belong to one partition only"""
    loci, seen = [], set()
    covered = np.zeros(L, dtype=bool)
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line or "=" not in line: continue
        lhs, rhs = line.split("=", 1)
        name = lhs.split(",")[-1].strip()
        if name in seen: sys.exit(f"partition name {name!r} appears twice")
        seen.add(name)
        sites = []
        for part in rhs.replace(";", "").split(","):
            part = part.strip()
            if not part: continue
            m = re.match(r"^(\d+)\s*-\s*(\d+)(?:\s*\\\s*(\d+))?$", part)
            if m:
                a, b, s = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
                if a < 1 or b < a or s < 1:
                    sys.exit(f"bad range {part!r} in partition {name}: sites are 1-based and a-b needs a <= b")
                sites.extend(range(a - 1, b, s))
            elif re.match(r"^\d+$", part):
                if int(part) < 1: sys.exit(f"bad site {part!r} in partition {name}: sites are 1-based")
                sites.append(int(part) - 1)
            else:
                sys.exit(f"cannot parse partition range {part!r} in {line!r}")
        if not sites: sys.exit(f"partition {name} has no sites")
        arr = np.array(sorted(set(sites)), dtype=int)
        if arr.max() >= L: sys.exit(f"partition {name} runs past the alignment ({arr.max() + 1} > {L})")
        if covered[arr].any(): sys.exit(f"partition {name} overlaps an earlier partition: a site may belong to one locus only")
        covered[arr] = True
        loci.append((name, arr))
    return loci


def read_input(path, partitions):
    """returns names, seqs, loci, format, datatype declared by a NEXUS file (or None), nexus"""
    head = open(path).read(64).lstrip()[:7].upper()
    nexus = None; declared = None
    if head.startswith("#NEXUS"):
        from Bio.Nexus import Nexus
        try:
            nexus = Nexus.Nexus(path)
        except Exception as e:
            sys.exit(f"Bio.Nexus could not parse {path}:\n  {e}\n"
                     "This tool reads NEXUS with Biopython and accepts only what Bio.Nexus accepts. Common causes: "
                     "commas inside CHARSET site lists (NEXUS wants whitespace: 'CHARSET x = 1-10 21-30;'), "
                     "taxon names with '-' or other punctuation that are not quoted, or unusual block formatting. "
                     "Check for commas in the SETS block, and/or supply the alignment as FASTA with a RAxML-style "
                     "partition file (--partitions).")
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
        if dt in ("dna", "rna", "nucleotide"): declared = "dna"
        elif dt == "protein": declared = "aa"
    else:
        if not partitions: sys.exit("a FASTA or PHYLIP alignment needs --partitions")
        txt = open(path).read()
        names, seqs = read_fasta(txt) if head.startswith(">") else read_phylip(txt)
        fmt = "fasta" if head.startswith(">") else "phylip"
        loci = parse_raxml_partitions(partitions, len(seqs[0]))
    L = len(seqs[0])
    bad = [t for t, s in zip(names, seqs) if len(s) != L]
    if bad: sys.exit(f"sequences differ in length (e.g. {bad[0]}): the input must be an alignment")
    dup = sorted(t for t, c in collections.Counter(names).items() if c > 1)
    if dup: sys.exit(f"duplicate taxon names: {dup[:5]}{' ...' if len(dup) > 5 else ''}")
    if not loci: sys.exit("no loci defined")
    return names, seqs, loci, fmt, declared, nexus


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


# ------------------------------------------------------------------ shared machinery
def wmedian_stack(vals, wts):
    """weighted median of vals over axis 0, nan-aware; (L, P) -> (P,)"""
    ok = np.isfinite(vals) & (wts > 0)
    v = np.where(ok, vals, np.inf); w = np.where(ok, wts, 0.0)
    order = np.argsort(v, axis=0)
    v_s = np.take_along_axis(v, order, 0); w_s = np.take_along_axis(w, order, 0)
    cw = np.cumsum(w_s, axis=0, dtype=np.float64); tot = cw[-1]
    idx = np.minimum((cw < tot[None, :] / 2).sum(0), v.shape[0] - 1)
    out = np.take_along_axis(v_s, idx[None, :], 0)[0]
    return np.where(tot > 0, out, np.nan)


_lgamma = np.vectorize(math.lgamma, otypes=[float])


def bb_tails(o, n, E, phi, upper):
    """log tail probability of the mismatch count among n shared sites: binomial with mean E, or, when the
    locus is overdispersed, beta-binomial with mean E and variance E + phi E^2 (the beta-binomial's
    intraclass correlation follows from phi in closed form). upper -> log P(X >= o), else log P(X <= o).
    Vectorised over pairs; the sum runs over the lower tail, which is short on the similar side."""
    o = np.asarray(o, float); n = np.asarray(n, float); E = np.asarray(E, float); out = np.empty(len(o))
    for i in range(len(o)):
        oi, ni, Ei = o[i], n[i], min(max(E[i], 1e-9), n[i] * (1 - 1e-9))
        d = Ei / ni
        kmax = int(oi) if not upper else int(oi) - 1                 # P(X <= kmax) summed exactly
        k = np.arange(0, max(kmax, 0) + 1, dtype=float)
        lchoose = _lgamma(ni + 1) - _lgamma(k + 1) - _lgamma(ni - k + 1)
        if phi < 1e-8 or ni < 2:
            lp = lchoose + k * math.log(d) + (ni - k) * math.log(1 - d)
        else:
            rho = min((d + phi * ni * d) / ((1 - d) * (ni - 1)), 1 - 1e-6)
            a_, b_ = d * (1 / rho - 1), (1 - d) * (1 / rho - 1)
            lp = lchoose + _lgamma(k + a_) + _lgamma(ni - k + b_) - _lgamma(ni + a_ + b_) - (math.lgamma(a_) + math.lgamma(b_) - math.lgamma(a_ + b_))
        mx = lp.max(); cdf = math.exp(mx) * float(np.exp(lp - mx).sum()) if kmax >= 0 else 0.0
        cdf = min(cdf, 1.0)
        out[i] = math.log(max(1.0 - cdf, 1e-300)) if upper else math.log(max(cdf, 1e-300))
    return out


def wmedian1(v, w):
    """weighted median of a 1-D array"""
    o = np.argsort(v); v, w = v[o], w[o]
    cw = np.cumsum(w); return float(v[min(int((cw < cw[-1] / 2).sum()), len(v) - 1)]) if cw[-1] > 0 else float(np.median(v))


def pair_median(vals, W, min_loci, fallback):
    """site-weighted median of vals across loci (axis 0), each locus weighted by the sites the
    pair shares in it. Pairs with fewer than min_loci contributing loci get `fallback`; also
    returns the count."""
    ok = np.isfinite(vals) & (W > 0)
    nl = ok.sum(0)
    v = wmedian_stack(vals, np.where(ok, W, 0.0))
    return np.where(nl >= min_loci, v, fallback), nl


def top_comparators(score, kk, nres):
    """rank comparators by score (higher = closer); returns top indices, validity, counts"""
    top = np.argpartition(-score, kk - 1, axis=1)[:, :kk]
    top = np.take_along_axis(top, np.argsort(-np.take_along_axis(score, top, 1), axis=1), 1)
    valid = np.isfinite(np.take_along_axis(score, top, 1)) & (nres > 0)[:, None]
    return top, valid, valid.sum(1)


def flatten_pairs(idx, top, valid, kk):
    gi = np.repeat(idx, kk); ci = top[idx].ravel(); vi = valid[idx].ravel(); rk = np.tile(np.arange(kk), len(idx))
    return gi[vi], ci[vi], rk[vi]


def norm_quantile(p):
    """inverse standard normal CDF, via bisection on erfc (no scipy)"""
    lo, hi = 0.0, 50.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.5 * math.erfc(mid / math.sqrt(2)) > 1 - p: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def lbinom(m_, o_, q):
    """binomial log-pmf, continuous in m_"""
    q = min(max(q, 1e-9), 1 - 1e-9)
    return (math.lgamma(o_ + 1) - math.lgamma(m_ + 1) - math.lgamma(o_ - m_ + 1)
            + m_ * math.log(q) + (o_ - m_) * math.log(1 - q))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True); ap.add_argument("--partitions"); ap.add_argument("--out", required=True)
    ap.add_argument("--datatype", required=True, choices=("dna", "aa"))
    ap.add_argument("--n-comp", type=int, default=5); ap.add_argument("--fit-overlap", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.05); ap.add_argument("--min-ratio", type=float, default=2.0)
    ap.add_argument("--mask-similar-only", action="store_true"); ap.add_argument("--mask-different-only", action="store_true")
    a = ap.parse_args(); K = a.n_comp
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    if a.mask_similar_only and a.mask_different_only:
        sys.exit("--mask-similar-only and --mask-different-only are mutually exclusive")
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    names, seqs, loci, fmt, declared, nexus = read_input(a.aln, a.partitions)
    if declared and declared != a.datatype:
        sys.exit(f"--datatype {a.datatype} was given but the NEXUS file declares {declared}")
    dtype = a.datatype
    if len(loci) < 10:
        sys.exit(f"only {len(loci)} loci defined: this is a phylogenomic test and needs at least 10")
    n, L = len(names), len(seqs[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    if n < 5: sys.exit(f"only {n} taxa: the test needs at least 5")
    K_eff = min(K, n - 1); K_eff -= 1 - K_eff % 2                    # largest odd <= min(K, n-1)
    M_eff = (K_eff + 1) // 2                                          # fails that make a majority
    if K_eff < K: print(f"n-comp shrunk to {K_eff} ({n} taxa)", file=sys.stderr)
    print(f"{n} taxa, {L} sites, {len(loci)} loci, {dtype}, {fmt}", file=sys.stderr)
    open(f"{a.out}/taxa.txt", "w").write("\n".join(names) + "\n")

    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = lut[np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(n, L)]

    # ---- pass 1: per-locus shared/same counts; alignment-wide identity falls out of the sums
    nloc = len(loci)
    S = np.empty((nloc, n, n), np.float32); W = np.empty((nloc, n, n), np.float32)
    plateaus = np.empty(nloc)
    for li, (name, sites) in enumerate(loci):
        Xm = X[:, sites]; Pm = (Xm >= 0).astype(np.float32)
        sh = Pm @ Pm.T; sm = np.zeros((n, n), np.float32)
        for i in range(len(alphabet)):
            Ai = (Xm == i).astype(np.float32); sm += Ai @ Ai.T
        S[li], W[li] = sm, sh
        cnt = np.bincount(Xm[Xm >= 0], minlength=len(alphabet)).astype(float)
        f = cnt / max(cnt.sum(), 1.0)
        plateaus[li] = max(100.0 * (1.0 - float((f ** 2).sum())), 1.0)
    present = np.stack([W[li].diagonal() >= a.fit_overlap for li in range(nloc)]).sum(0)
    thin = np.flatnonzero(present < MIN_LOCI_PER_TAXON)
    open(f"{a.out}/excluded_taxa.txt", "w").write("".join(f"{names[i]}\t{int(present[i])}\n" for i in thin))
    if len(thin):
        S[:, thin, :] = 0; S[:, :, thin] = 0; W[:, thin, :] = 0; W[:, :, thin] = 0   # out of the test entirely
        print(f"{len(thin)} taxa with fewer than {MIN_LOCI_PER_TAXON} loci present are excluded from the test "
              f"(listed in excluded_taxa.txt, left unmasked)", file=sys.stderr)
    tot_same = S.sum(0, dtype=np.float64); tot_shared = W.sum(0, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(tot_shared > 0, 100.0 * tot_same / tot_shared, np.nan)
    np.fill_diagonal(ani, np.nan); np.save(f"{a.out}/ani.npy", ani)
    print(f"per-locus identities done, {time.time() - t0:.0f}s", file=sys.stderr)

    # ---- pair divergence m: capped site-weighted median across loci of rate-normalised divergence
    with np.errstate(invalid="ignore", divide="ignore"):
        Dv = np.where(W > 0, 100.0 * (1.0 - S / np.maximum(W, 1)), np.nan)
    Wm = W.copy()
    Wm[:, np.arange(n), np.arange(n)] = 0.0                            # self pairs carry no vote
    nan_fb = np.full(n * n, np.nan)
    # ONE rate per locus, estimated once and used everywhere (to normalise loci inside m and as the
    # expectation in the test). It is estimated where the test looks: each taxon's nearest neighbours
    # by the provisional pair divergence m0 (the raw site-weighted median over loci, available before rates exist), one equal
    # vote per neighbour pair resolved beyond its floor, of the multiple-hit-corrected ratio
    # c(d + half a mismatch) / c(m0), m0 the provisional pair divergence (raw site-weighted median).
    Wfl = 100.0 * HALF / np.maximum(W, 1.0)                          # the imputed mismatches per locus, as a percentage: the resolution floor
    m0, _ = pair_median(Dv.reshape(nloc, -1), Wm.reshape(nloc, -1), MIN_LOCI_FOR_M, nan_fb)
    m0_fl, _ = pair_median(Wfl.reshape(nloc, -1), Wm.reshape(nloc, -1), MIN_LOCI_FOR_M, nan_fb)
    m0, m0_fl = m0.reshape(n, n), m0_fl.reshape(n, n)
    m0_hi = np.maximum(m0, m0_fl); res0 = m0 > m0_fl
    dist = np.where(np.isfinite(m0), m0, np.inf); np.fill_diagonal(dist, np.inf)   # neighbours by the provisional pair divergence (a median over loci, not swayed by one long fast locus)
    Kn = min(K_eff, n - 1)
    nb = np.argsort(dist, axis=1)[:, :Kn]
    gi0 = np.repeat(np.arange(n), Kn); ci0 = nb.ravel()
    okn = np.isfinite(dist[gi0, ci0]) & res0[gi0, ci0]
    gi0, ci0 = gi0[okn], ci0[okn]
    m0n = m0_hi[gi0, ci0]
    rate = np.full(nloc, np.nan)
    for li in range(nloc):
        o_, d_ = W[li][gi0, ci0], Dv[li][gi0, ci0]
        val = np.isfinite(d_) & (o_ >= a.fit_overlap)
        if val.sum() < MIN_FIT_PAIRS: continue
        P_ = plateaus[li]
        corr = lambda p: -P_ * np.log(1.0 - np.minimum(p, 0.98 * P_) / P_)
        rate[li] = float(np.median(corr(d_[val] + 100.0 * HALF / o_[val]) / corr(m0n[val])))
    lsites = np.array([len(s) for _, s in loci], dtype=float)
    ok_r = np.isfinite(rate) & (rate > 0)
    if ok_r.sum() >= 3:
        scale = wmedian_stack(rate[ok_r, None], lsites[ok_r, None])[0]
        rate = rate / max(scale, 1e-6)
    rr = np.where(ok_r, rate, np.nan)
    m, nl = pair_median((Dv / rr[:, None, None]).reshape(nloc, -1),
                        np.where(ok_r[:, None, None], Wm, 0.0).reshape(nloc, -1),
                        MIN_LOCI_FOR_M, nan_fb)
    m = m.reshape(n, n)
    np.fill_diagonal(m, np.nan); np.save(f"{a.out}/m.npy", m)
    # resolution floor of the pair: the median over its loci of half a mismatch per locus (rate-normalised);
    # the divergence used on the too-divergent side is the larger of the measured median and this floor
    m_fl, _ = pair_median((Wfl / rr[:, None, None]).reshape(nloc, -1),
                          np.where(ok_r[:, None, None], Wm, 0.0).reshape(nloc, -1), MIN_LOCI_FOR_M, nan_fb)
    m_fl = m_fl.reshape(n, n)
    m_hi = np.maximum(m, m_fl); np.fill_diagonal(m_hi, np.nan); np.save(f"{a.out}/m_hi.npy", m_hi)
    m_xs = np.maximum(m - m_fl, 0.0)                                  # divergence resolved beyond the floor: the rate-vote weight
    keep = np.ones(n, bool); keep[thin] = False
    nlk = nl.reshape(n, n)[np.ix_(keep, keep)][np.triu_indices(int(keep.sum()), 1)]
    fb = float((nlk < MIN_LOCI_FOR_M).mean()) if len(nlk) else 0.0
    print(f"pair divergences done, {time.time() - t0:.0f}s"
          + (f"; note: {fb:.0%} of taxon pairs overlap in <{MIN_LOCI_FOR_M} loci and are never compared"
             if fb > 0.05 else ""), file=sys.stderr)
    del Dv, Wfl, Wm

    # deterministic tie-break for comparator ranking: pooled divergence, then taxon order
    tie = np.nan_to_num(100.0 - ani, nan=100.0) * 1e-9 + np.arange(n)[None, :] * 1e-12
    masked = [bytearray(s.encode()) for s in seqs]

    nmasked = run_z(a, names, loci, X, S, W, m, m_hi, m_xs, plateaus, rate, tie, masked, K_eff, M_eff)

    out = write_masked(a.aln, a.out, fmt, names, [mb.decode() for mb in masked], nexus)
    print(f"masked {nmasked} cells -> {out}\ndone {time.time() - t0:.0f}s", file=sys.stderr)


# ------------------------------------------------------------------ the counting test
def run_z(a, names, loci, X, S, W, m, m_hi, m_xs, plateaus, rate, tie, masked, K_eff, M_eff):
    n = len(names); nloc = len(loci)
    mj = np.where(np.isfinite(m), m, np.inf) + tie
    mask_hi = not a.mask_similar_only
    mask_lo = not a.mask_different_only
    views = [np.frombuffer(mb, dtype=np.uint8) for mb in masked]        # writable views for blanking

    # pre-pass: count distinct comparisons so the Bonferroni bar is known before any test runs
    Npairs = 0; sel = []
    for li in range(nloc):
        ov = W[li]
        nres = ov.diagonal().astype(int)
        need = np.minimum(a.fit_overlap, nres)[:, None]
        score = np.where((ov >= need) & np.isfinite(mj), -mj, -np.inf)
        np.fill_diagonal(score, -np.inf)
        top, valid, ncomp = top_comparators(score, K_eff, nres)
        sel.append((top, valid, ncomp, nres))
        idx = np.flatnonzero(ncomp > 0)
        if len(idx):
            gi, ci, _ = flatten_pairs(idx, top, valid, K_eff)
            Npairs += len(np.unique(np.minimum(gi, ci) * n + np.maximum(gi, ci)))
    kq = norm_quantile(1 - a.alpha / max(2 * Npairs, 1))              # two one-sided tests per pair
    print(f"{Npairs} distinct comparisons, two-sided; Bonferroni k = {kq:.2f} at alpha {a.alpha}", file=sys.stderr)

    pairs = open(f"{a.out}/pairs.tsv", "w")
    pairs.write("locus\ttaxon\trank\tcomparator\tpair_div\toverlap\tidentity\texpected\texcess\tfail\tfail_similar\n")
    cells = open(f"{a.out}/cells.tsv", "w")
    cells.write("locus\ttaxon\tn_comp\tcomparator\tpair_div\tfocal_sites\toverlap\tidentity\texpected\texcess"
                "\tmed_identity\tmed_expected\tmed_excess\tn_fail\tn_similar\tflag\tambiguous\texcess_m\n")
    per = open(f"{a.out}/per_locus.tsv", "w")
    per.write("locus\tsites\ttaxa\ttested\tpairs\trate\tplateau\tcut\tflagged\t"
              + "\t".join(f"fail{j}" for j in range(M_eff, K_eff + 1)) + "\tambiguous\tsimilar\tstatus\n")
    disp = open(f"{a.out}/dispersion.tsv", "w"); disp.write("locus\tphi_hat\tpairs\tphi\n")
    trends = open(f"{a.out}/trends.tsv", "w"); trends.write("locus\tpair_div\ttrend\tcutoff\tcutoff_hi\n")
    simf = open(f"{a.out}/similar.tsv", "w")
    simf.write("locus\ttaxon_a\ttaxon_b\tpair_div\toverlap\tidentity\texpected\tdeficit\tlambda\tdonor\trecipient\tmasked\n")
    grid_m = np.logspace(np.log10(0.05), np.log10(80.0), 200)
    mrows = {}                                                        # (locus, taxon) -> masked.csv row

    def mask_cell(name, g, sites, side, evidence):
        key = (name, g)
        if key in mrows:
            if mrows[key][3] != side: mrows[key][3] = "both"
            return
        v = views[g][sites]; cnt = int((v != 45).sum()); views[g][sites] = 45   # 45 = '-'
        mrows[key] = [name, names[g], cnt, side] + evidence

    def locus_pairs(li):
        """the comparator pairs of a locus and the locus's rate; None when the locus is skipped"""
        ov, sm = W[li], S[li]
        top, valid, ncomp, nres = sel[li]
        idx = np.flatnonzero(ncomp > 0)
        if len(idx) == 0: return "skipped-no-cells"
        gi, ci, rk = flatten_pairs(idx, top, valid, K_eff)
        o = ov[gi, ci]; same = sm[gi, ci]; mis = o - same
        ident = 100.0 * same / o; dob = 100.0 - ident
        mp = m[gi, ci]; mph = m_hi[gi, ci]; mpx = m_xs[gi, ci]
        fit = (o >= a.fit_overlap) & (mpx > 0)
        if fit.sum() < MIN_FIT_PAIRS: return f"skipped-few-resolved-pairs({int(fit.sum())})"
        P_ = plateaus[li]
        corr = lambda p: -P_ * np.log(1.0 - np.minimum(np.asarray(p, float), 0.98 * P_) / P_)
        pc = 100.0 * HALF / o
        rate_c = max(float(np.median((corr(dob + pc) / corr(mph))[fit])), 1e-3)
        def expect(mv):
            mvc = np.maximum(np.asarray(mv, float), 1e-4)
            return np.maximum(P_ * (1.0 - np.exp(-rate_c * corr(mvc) / P_)), 1e-4)
        return dict(idx=idx, gi=gi, ci=ci, rk=rk, o=o, same=same, mis=mis, ident=ident, dob=dob, mp=mp, mph=mph, fit=fit,
                    ov=ov, sm=sm, ncomp=ncomp, nres=nres, expect=expect, corr=corr, rate_c=rate_c)

    # dispersion pre-pass: the robust spread of observed/expected among each locus's resolved comparator
    # pairs with enough expected mismatches, with the Poisson part removed: Var(O/E) = 1/E + phi
    phi_hat = np.full(nloc, np.nan); phi_n = np.zeros(nloc)
    for li in range(nloc):
        L = locus_pairs(li)
        if isinstance(L, str): continue
        E = L["o"] * L["expect"](L["mph"]) / 100.0
        use = L["fit"] & (E >= PHI_MIN_E)
        if use.sum() < 10: continue
        r = np.log((L["mis"][use] + HALF) / E[use])
        s_ = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        phi_hat[li] = max(s_ ** 2 - float(np.median(1.0 / E[use])), 0.0); phi_n[li] = use.sum()
    okp = phi_n > 0
    phi_bar = wmedian1(phi_hat[okp], phi_n[okp]) if okp.any() else 0.0
    print(f"dispersion: dataset-wide phi {phi_bar:.4f} (robust SD of O/E beyond Poisson {math.sqrt(phi_bar):.2f}); "
          f"per-locus phi 10/50/90%: {np.nanpercentile(phi_hat[okp], [10, 50, 90]).round(4).tolist() if okp.any() else 'n/a'}", file=sys.stderr)

    for li, (name, sites) in enumerate(loci):
        L = locus_pairs(li)
        skip = L if isinstance(L, str) else None
        if not skip:
            idx, gi, ci, rk, o, same, mis, ident, dob, mp, mph, fit, ov, sm, ncomp, nres, expect, corr, rate_c = \
                (L[k] for k in ("idx", "gi", "ci", "rk", "o", "same", "mis", "ident", "dob", "mp", "mph", "fit", "ov", "sm", "ncomp", "nres", "expect", "corr", "rate_c"))
            phi = (phi_n[li] * (phi_hat[li] if okp[li] else 0.0) + PHI_N0 * phi_bar) / (phi_n[li] + PHI_N0)
        if skip:
            per.write(f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t0\t0\t"
                      f"{rate[li]:.3f}\t{plateaus[li]:.1f}\t{kq:.2f}\t0\t"
                      + "\t".join("0" for _ in range(M_eff, K_eff + 1)) + f"\t0\t0\t{skip}\n")
            print(f"{name}: {skip}", file=sys.stderr); continue

        # the count test with a quasi-Poisson variance: Var(O) = E + phi E^2, phi the locus's dispersion
        # (shrunk toward the dataset's); phi = 0 is the plain Poisson test
        dhat = expect(mph)                                          # divergent side: the divergence we cannot exclude
        dlo = expect(mp)                                            # similar side: the divergence we can
        E = o * dhat / 100.0
        E_lo = o * dlo / 100.0
        # exact binomial tails over the shared sites (beta-binomial when the locus is overdispersed), evaluated
        # only for the pairs that clear the effect floors; a pair fails when its tail probability is below alpha / 2N
        excess = mis - E
        logbar = math.log(a.alpha / max(2 * Npairs, 1))
        fail_hi = excess > (a.min_ratio - 1) * E
        if fail_hi.any():
            fail_hi[fail_hi] = bb_tails(mis[fail_hi], o[fail_hi], E[fail_hi], phi, True) < logbar
        defi = E_lo - mis
        fail_lo = defi > (2.0 / 3.0) * E_lo
        if fail_lo.any():
            fail_lo[fail_lo] = bb_tails(mis[fail_lo], o[fail_lo], E_lo[fail_lo], phi, False) < logbar

        pairs.writelines(f"{name}\t{names[gi[j]]}\t{rk[j] + 1}\t{names[ci[j]]}\t{mph[j]:.3f}\t{int(o[j])}"
                         f"\t{ident[j]:.2f}\t{100 - dhat[j]:.2f}\t{excess[j]:.1f}\t{int(fail_hi[j])}\t{int(fail_lo[j])}\n"
                         for j in range(len(gi)))

        # ---- divergent side: majority over comparators
        starts = np.r_[0, np.cumsum(ncomp[idx])[:-1]]
        fcount = {j: 0 for j in range(M_eff, K_eff + 1)}; namb = 0
        crow = []
        for g, s0, kc in zip(idx, starts, ncomp[idx]):
            slc = slice(s0, s0 + kc)
            nf = int(fail_hi[slc].sum()); nsim = int(fail_lo[slc].sum())
            flag = int(kc >= K_eff and nf >= M_eff)
            amb = int(kc < K_eff and nf >= kc // 2 + 1)
            if flag: fcount[nf] += 1
            namb += amb
            med_i = float(np.median(ident[slc])); med_e = float(np.median(100 - dhat[slc])); med_x = float(np.median(excess[slc]))
            # the cell's extra divergence beyond expectation, as corrected distance in units of m (a typical
            # locus's divergence): a foreign or paralogous sequence adds a constant here whatever the comparator
            exm = float((corr(100 - med_i) - corr(100 - med_e)) / rate_c)
            crow.append(f"{name}\t{names[g]}\t{kc}\t{names[ci[s0]]}\t{mph[s0]:.3f}\t{nres[g]}\t{int(o[s0])}"
                        f"\t{ident[s0]:.2f}\t{100 - dhat[s0]:.2f}\t{excess[s0]:.1f}\t{med_i:.2f}\t{med_e:.2f}\t{med_x:.1f}"
                        f"\t{nf}\t{nsim}\t{flag}\t{amb}\t{exm:.3f}\n")
            if (flag or amb) and mask_hi:
                mask_cell(name, g, sites, "divergent",
                          [nf, nsim, names[ci[s0]], f"{mph[s0]:.3f}", f"{med_i:.2f}", f"{med_e:.2f}", f"{med_x:.1f}", "", f"{exm:.3f}"])
        cells.writelines(crow)

        # ---- similar side: any failing pair; LR direction attribution on shared sites
        simpairs = {}
        for j in np.flatnonzero(fail_lo):
            key = (int(min(gi[j], ci[j])), int(max(gi[j], ci[j])))
            simpairs.setdefault(key, j)
        Xs = X[:, sites] if simpairs else None
        for (A, B) in sorted(simpairs):
            lam, donor, recip = attribute(A, B, Xs, m, ov, expect, a.fit_overlap)
            o_ab = float(ov[A, B]); id_ab = 100.0 * sm[A, B] / max(o_ab, 1.0)
            exp_ab = float(expect(m[A, B])); defi_ab = o_ab * exp_ab / 100 - (o_ab - sm[A, B])
            targets = ([recip] if recip is not None else [A, B]) if mask_lo else []
            side = ("both" if recip is None else names[recip]) if mask_lo else "none"
            for g in targets:
                mask_cell(name, g, sites, "similar",
                          ["", "", names[B if g == A else A], f"{m[A, B]:.3f}", f"{id_ab:.2f}", f"{100 - exp_ab:.2f}",
                           f"{-defi_ab:.1f}", f"{lam:.1f}", f"{float((corr(100 - id_ab) - corr(exp_ab)) / rate_c):.3f}"])
            simf.write(f"{name}\t{names[A]}\t{names[B]}\t{m[A, B]:.3f}\t{int(o_ab)}\t{id_ab:.2f}\t{100 - exp_ab:.2f}"
                       f"\t{defi_ab:.1f}\t{lam:.1f}\t{names[donor] if donor is not None else 'NA'}"
                       f"\t{names[recip] if recip is not None else 'NA'}\t{side}\n")

        gd = expect(grid_m)
        trends.writelines(f"{name}\t{gm:.4f}\t{100 - gdv:.4f}\t{100 - min(a.min_ratio * gdv, plateaus[li]):.4f}\t{100 - gdv / 3:.4f}\n"
                          for gm, gdv in zip(grid_m, gd))
        nflag = sum(fcount.values()); fstr = "/".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
        per.write(f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t{len(idx)}\t{len(gi)}\t"
                  f"{rate[li]:.3f}\t{plateaus[li]:.1f}\t{kq:.2f}\t{nflag}\t"
                  + "\t".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
                  + f"\t{namb}\t{len(simpairs)}\ttested\n")
        disp.write(f"{name}\t{(phi_hat[li] if okp[li] else float('nan')):.4f}\t{int(phi_n[li])}\t{phi:.4f}\n")
        print(f"{name}: {len(idx)} cells, {len(gi)} pairs, flagged {nflag} "
              f"(fails {M_eff}..{K_eff} of {K_eff}: {fstr}), ambiguous {namb}, similar {len(simpairs)}", file=sys.stderr)
    for fh in (pairs, cells, per, trends, simf, disp): fh.close()
    with open(f"{a.out}/masked.csv", "w", newline="") as fh:
        mw = csv.writer(fh)
        mw.writerow(["locus", "taxon", "sites_masked", "side", "n_fail", "n_similar", "partner",
                     "pair_div", "identity", "expected", "excess", "lambda", "excess_m"])
        mw.writerows(mrows.values())
    return len(mrows)


def attribute(A, B, Xs, m, ov, expect, fit_overlap, nthird=10):
    """direction of a similar pair: likelihood ratio over up to nthird third parties, with both
    hypotheses evaluated on the SAME sites (those shared by A, B and C), so the comparison is
    not driven by how much data each party has. The shared sequence should fit the DONOR's
    expectations. |2 ln LR| > LR_BAR attributes; returns (lambda, donor, recipient) or
    (lambda, None, None) when unidentifiable."""
    mA, mB = m[A], m[B]
    cand = np.flatnonzero(np.isfinite(mA) & np.isfinite(mB) & (ov[A] >= fit_overlap) & (ov[B] >= fit_overlap))
    cand = cand[(cand != A) & (cand != B)]
    if len(cand) < 3: return 0.0, None, None
    cand = cand[np.argsort(np.minimum(mA[cand], mB[cand]))[:nthird]]
    xa, xb = Xs[A], Xs[B]; both = (xa >= 0) & (xb >= 0)
    LL = {A: 0.0, B: 0.0}; used = 0
    for C in cand:
        xc = Xs[C]; ok = both & (xc >= 0); o_ = float(ok.sum())
        if o_ < fit_overlap: continue
        used += 1
        for Xi, x in ((A, xa), (B, xb)):
            LL[Xi] += lbinom(float((x[ok] != xc[ok]).sum()), o_, float(expect(m[Xi, C])) / 100.0)
    if used < 3: return 0.0, None, None
    lam = 2.0 * (LL[A] - LL[B])
    if abs(lam) <= LR_BAR: return lam, None, None
    return (lam, A, B) if lam > 0 else (lam, B, A)


if __name__ == "__main__":
    main()
