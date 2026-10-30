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
that are left untouched in the output and named in the log and the report. At least 5 taxa.

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
  identical a genome identical to the focal at every site the two share is never one of its
            comparators, never its closest match and never a rate voter, and no two comparators
            are identical to each other (0.13.1): copies of a sequence carry its errors and
            would vouch for each other, and five copies of one genome are one comparison, not
            five. A genome identical at
            most loci but different at some is a comparator like any other, tested at the
            pair's floor, which is where such differences stand out most.
  rates     a locus rate is the median, over a set of taxon pairs, one equal vote per pair,
            of c(d + half a mismatch) / c(m), where d is the pair's mismatch fraction at the
            locus, c(p) = -P ln(1 - p/P) the multiple-hit correction with P the locus's
            random-identity plateau (from its residue frequencies), and m the pair's
            divergence (or its floor); rates are rescaled so the median locus has rate 1.
            Rates are estimated twice. Before m exists, each taxon's --n-comp nearest
            RESOLVED neighbours by the raw site-weighted median of locus divergences vote
            (0.13.1: unresolved and identical relatives are skipped, so a dataset of
            near-identical sequences still gets its rates from the pairs that differ); these provisional rates
            normalise the loci when m is computed. In the test, each taxon's --n-comp closest
            comparators by m vote, resolved pairs only; that rate gives the expectation. Where
            a locus's comparator pairs cannot set a rate (fewer than 8 resolved: a clade of
            identical sequences) the provisional rate stands in, so the contaminant among
            identical relatives is still tested, at the pair's floor.
            Nearest-neighbour pairs are used because test 1 compares close relatives, and a
            median over all pairs is a median over the deep branches of the tree, where a
            locus's rate can differ several-fold from its rate near the tips. Test 2 uses the
            same expectation for every pair, deep ones included: at a fast locus (a
            mitochondrial gene several times the typical rate) the expectation for deep pairs
            approaches the plateau faster than real sequences do (no rate variation among
            sites), so deep pairs there can look too similar.
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
            among the o shared sites, the observed count is binomial(o, E/o), or beta-binomial
            with the same mean and variance E + phi E^2 when the locus is dispersed. Every
            sequence (taxon x locus) gets two tests.
            TEST 1, too divergent: each comparison with one of its comparators fails when
            P(X >= observed) is below the bar; the sequence fails when a majority of its
            comparators do (fewer comparators than requested failing a majority = "ambiguous").
            TEST 2, too similar: is the sequence closer to another sequence than it should be?
            Two arms, one decision. (a) Against each comparator, lower tail: the expectation is
            the locus curve at the pair's m (not the floor: a pair cannot be more similar than
            its loci allow), held against E and E x f, f <= 1 the row factor, the median over
            each taxon's other comparators of (observed + 0.5) / (expected + 0.5), the smaller of
            the pair's two: a lineage that is slow at the locus is close to everything there,
            and that is not a copy. (b) Against its CLOSEST MATCH at the locus (fewest mismatches
            per shared site, half a mismatch added, ties by m) when that match is not one of its
            comparators: is it closer to the match than its own relatives are? The expected
            mismatch rate is the median, over its comparators that pass test 1 here and share
            >= --fit-overlap sites with the match (at least three), of their rate to the match
            at this locus: no rate model and no m, so a conserved locus, a slow lineage or a
            deep pair's saturation move the relatives as much as the sequence. A copy from
            outside its five closest genomes shows here; one from a close relative in (a).
            A failing pair is too similar. Its COPY is the member that fits its own relatives
            worse: the larger |test-1 excess| (median over its comparators of observed minus
            expected divergence, in units of m; absolute, because a foreign sequence in a fast
            lineage is CLOSER to the lineage's relatives than expected). The copy is masked; the
            other member is its source and is not masked for it.
            The bar is alpha / N, N every test the run performs: one too-divergent and one
            too-similar test per distinct comparator pair, and one closest-match test per
            sequence tested. Bonferroni over the whole run. There is no effect-size floor: the
            locus's dispersion sets how far from expectation a long locus has to be.
  decide    a sequence is MASKED when it fails test 1, or when it is the copy of a too-similar
            pair (or a member of one whose direction is uncertain), or both. --mask-similar-only / --mask-different-only restrict masking to one
            test. Masking blanks the sequence's whole stretch of the locus.

Outputs (in --out; <stem> is the alignment file's name without extension)
  <stem>.masked.<ext>  the alignment in its input format with every masked sequence blanked
  masked.csv           ONE row per masked sequence and why: reason (divergent, similar or both),
                       residues_masked (its residues at the locus; the whole stretch is blanked), and
                       the evidence of both tests. Test 1: divergent_fails of divergent_comparators,
                       the observed and expected divergence (%; medians over the comparators),
                       divergent_excess_m (the extra divergence beyond expectation, multiple-hit
                       corrected and divided by the locus rate: % on the scale of m) and
                       divergent_log10p, the p-value of the comparison that decides the majority.
                       Test 2: similar_test (comparator or closest: the arm whose test fired lowest, or
                       the closest-match test the sequence got when none fired), similar_partner (the
                       other member of the pair, or the closest match) and similar_role (copy, source or
                       uncertain: the direction is called when the larger of the two members' test-1
                       excesses is more than --donor-threshold times the smaller, else both are masked;
                       empty when no pair fired), the pair's observed and expected divergence (%; the
                       expectation before the row factor), similar_row_factor (comparator arm),
                       similar_log10p, similar_relatives (closest-match arm: the comparators whose rate
                       to the match set the expectation) and partner_excess_m (the partner's test-1
                       excess, next to this sequence's divergent_excess_m: the margin of the copy call).
                       Then its comparators, closest first, and log10_bar, the bar every p-value was
                       held to. p-values are floored at 1e-300 (log10 -300).
  tests.tsv.gz         with --all-tests only: one row per test performed, the focal sequence first.
                       test (divergent, similar_comparator, similar_closest), locus, taxon, comparator,
                       overlap (shared sites), mismatches, expected (expected mismatches; before the row
                       factor), row_factor (similar_comparator rows), log10p (exact), fail (0 or 1).
                       Large (every comparison at every locus); tabs or line breaks in names are
                       written as spaces.
  report.html          the interactive report, one file with its data inside (report.npz holds the
                       same data for scripts). --no-report skips it.
Options
  --datatype dna|aa (required)  --n-comp 5 (odd >= 3)  --fit-overlap 20
  --alpha 0.05  --mask-similar-only  --mask-different-only  --all-tests  --no-report
  --threads 1 (0: every CPU this process may use; N: at most N; each stage takes as many as pay off)
  --memory GB (default: half of what is available when the run starts; 0: no budget; a run that
  cannot fit stops before it starts and says what it needs)
  (constants in the code: HALF = 0.5 mismatch for the resolution floor and the rate
  pseudocount; dispersion from pairs with >= 10 expected mismatches, shrinkage weight 30)
"""
import argparse, base64, collections, csv, gzip, hashlib, io, json, math, os, platform, re, subprocess, sys, threading, time, zipfile
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):   # BLAS stays single-threaded (its matrices are
    os.environ.setdefault(_v, "1")                                          # small); --threads parallelises over loci and pairs
import numpy as np
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from itertools import chain

PROTEIN = "ACDEFGHIKLMNPQRSTVWY"
DNA = "ACGT"
HALF = 0.5              # mismatches imputed to a locus that shows none (half the detection limit);
                        # also the pseudocount on the locus side of the rate ratio
PHI_MIN_E = 10.0        # a pair contributes to the locus dispersion only if it expects this many mismatches
PHI_N0 = 30.0           # shrinkage weight: a locus's dispersion is pulled toward the dataset's with this many pairs' worth
MIN_LOCI_FOR_M = 5      # pairs overlapping in fewer loci have no m and are never compared
MIN_LOCI_PER_TAXON = 10 # a taxon needs this many loci with >= fit-overlap residues to take part
MIN_FIT_PAIRS = 8       # a locus needs this many fit pairs to be tested (or to get a rate)
TEST_THREAD_PAIRS = 4000  # the dispersion pre-pass uses threads only when a locus has this many comparisons on average
TEST_BATCH_COMPARISONS = 50_000   # the test computes the exact tails of this many comparisons at a time (loci batched), so small
                                  # loci do not pay numpy's per-call overhead hundreds of times each; about 250 MB of temporaries per batch
TEST_BATCH_MIN = 5_000            # with worker processes the batches shrink to keep every worker busy, but never below this
PAIRS_PER_THREAD = 4000           # the pair stage uses a thread per this many pairs: beyond that the threads contend (Shen_2018,
                                  # 59K pairs x 2,408 loci: 16 threads 1.3 s, 32 threads 3.8 s, 96 threads 5.2 s; Ivan_2026, 246K
                                  # pairs: 96 threads 0.9 s)
LSET_ELEMENTS = 1 << 20           # a locus's closest matches are found in taxa x taxa slices of about this many elements (rows = this / taxa)
__version__ = "0.13.1"
REPORT_FORMAT = 7                 # bundle format version; the viewer refuses a newer one (3: tree branch lengths are F81-corrected; 4: the too-similar test over every pair; 7: every test of every cell, the closest match at the locus, the uncertain-direction pairs)
DENSE_CELLS_MAX = 8_000_000       # per-cell numeric report columns are dense (taxa x loci) up to this many cells, sparse beyond
BEST_MATCH_ALL_MAX_TAXA = 2000    # best matches for every tested cell up to this many taxa; beyond, only for cells with a fail
BEST_MATCH_MIN_SITES = 30         # a best match needs this many sites shared with the focal (as plot_panels.best_matches)
REPORT_SEQ_ALL_MAX = 32 << 20     # --report-sequences auto carries the whole alignment up to this many bytes, else the cells' rows
CODE_SIMILAR, CODE_TICK, CODE_ABSENT = 128, 254, 255
PASS1_FLOPS_PER_THREAD = 6e7    # pass one adds a thread per this many floating-point operations per locus: below it the
                                # interpreter lock around the small matrix products costs more than the thread gains
                                # (measured: 90 taxa x 590 sites, 4e7 flops per locus with the per-letter products, is
                                # fastest at 1-2 threads and slower at 4 or more)


# ------------------------------------------------------------------ input
def read_fasta(txt):
    """names and sequences, one per header line, in order: a header with nothing under it gives an empty sequence
    (read_input then rejects it by name); text before the first header is an error, not a nameless sequence"""
    names, seqs, cur = [], [], None
    for line in txt.splitlines():
        if line.startswith(">"):
            if cur is not None: seqs.append("".join(cur))
            name = line[1:].split()
            if not name: sys.exit("FASTA: a header line ('>') has no name")
            names.append(name[0]); cur = []
        elif cur is None:
            if line.strip(): sys.exit(f"FASTA: text before the first header line ('>'): {line.strip()[:40]!r} (a sequence whose header line is missing?)")
        else:
            cur.append("".join(line.split()))
    if cur is not None: seqs.append("".join(cur))
    return names, seqs


def read_phylip(txt):
    """relaxed PHYLIP: after the counts line, 'name sequence' lines (a sequence may continue on the following lines),
    or interleaved blocks with the names in the first block only or repeated in every block"""
    lines = [l for l in txt.splitlines() if l.strip()]
    hdr = lines[0].split() if lines else []
    if len(hdr) < 2 or not hdr[0].isdigit() or not hdr[1].isdigit():
        sys.exit("PHYLIP: the file is empty" if not lines else
                 f"cannot parse PHYLIP: the first line must give the number of sequences and the number of sites, not {lines[0].strip()[:40]!r}")
    n, L = int(hdr[0]), int(hdr[1]); body = lines[1:]; squash = lambda x: "".join(x.split())
    def sequential():
        names, seqs = [], []
        for line in body:
            p = line.split(None, 1)
            if len(p) == 2 and len(names) < n: names.append(p[0]); seqs.append(squash(p[1]))
            elif names: seqs[-1] += squash(line)                 # a sequence continued on the next line
        return names, seqs
    def interleaved():
        names, seqs = [], []
        for i, line in enumerate(body):
            p = line.split(None, 1)
            if i < n: names.append(p[0]); seqs.append(squash(p[1]) if len(p) == 2 else "")
            else: seqs[i % n] += squash(p[1] if len(p) == 2 and p[0] == names[i % n] else line)   # the names may be repeated in every block
        return names, seqs
    names, seqs = sequential()
    if len(names) != n or any(len(sq) != L for sq in seqs): names, seqs = interleaved()
    if len(names) == n and any(len(sq) != L for sq in seqs):
        t = next(t for t, sq in zip(names, seqs) if len(sq) != L); l = len(next(sq for sq in seqs if len(sq) != L))
        sys.exit(f"PHYLIP: the header says {L:,} sites but {t} has {l:,}: the sequences are not aligned (every sequence must have {L:,} characters), "
                 "or the file is laid out in a form the reader does not handle (relaxed one-line and interleaved layouts are supported; try FASTA plus a partition file)")
    if len(names) != n or any(len(sq) != L for sq in seqs):
        sys.exit(f"cannot parse PHYLIP: the header says {n} sequences of {L:,} sites, {len(names)} were found; relaxed one-line or interleaved layouts are supported; "
                 "try FASTA plus a RAxML-style partition file")
    return names, seqs


def parse_raxml_partitions(path, L):
    """name -> 0-based site index array; ranges 'a-b', 'a-b\\s' (step) and comma lists;
    every site may belong to one partition only"""
    loci, seen = [], set()
    owner = np.full(L, -1, dtype=np.int64)                           # the partition that holds each site, for the overlap message
    try:
        with open(path, encoding="utf-8") as fh: lines = fh.read().splitlines()
    except UnicodeDecodeError: sys.exit(f"the partition file {path} is not a text file")
    except OSError as e: sys.exit(f"cannot read the partition file {path}: {e.strerror or e}")
    for line in lines:
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
                if b > L: sys.exit(f"partition {name} runs past the alignment ({b:,} > {L:,} sites)")
                sites.extend(range(a - 1, b, s))
            elif re.match(r"^\d+$", part):
                if int(part) < 1: sys.exit(f"bad site {part!r} in partition {name}: sites are 1-based")
                if int(part) > L: sys.exit(f"partition {name} runs past the alignment ({int(part):,} > {L:,} sites)")
                sites.append(int(part) - 1)
            else:
                sys.exit(f"cannot parse partition range {part!r} in {line!r}")
        if not sites: sys.exit(f"partition {name} has no sites")
        arr = np.array(sorted(set(sites)), dtype=int)
        if arr.max() >= L: sys.exit(f"partition {name} runs past the alignment ({arr.max() + 1} > {L})")
        dup = arr[owner[arr] >= 0]
        if len(dup):
            first = int(dup[0]); others = sorted({loci[int(owner[x])][0] for x in dup})
            sys.exit(f"site {first + 1} is in partition {name} and in partition {others[0]}: a site may be in one partition only "
                     f"({len(dup)} site{'s' if len(dup) > 1 else ''} of {name} {'are' if len(dup) > 1 else 'is'} already in {', '.join(others)}; "
                     f"the alignment has {L} sites, a site may appear only once)")
        owner[arr] = len(loci)
        loci.append((name, arr))
    return loci


class _Unsupported(Exception):
    """the fast NEXUS reader met something it does not handle; Bio.Nexus takes over"""


# the characters an alignment may hold, as Bio.Nexus accepts them: the residues, the IUPAC ambiguity codes (missing data to
# the test), '*' for protein, and the gap and missing symbols ('-' and '?', or what a NEXUS file declares); U is read as T
ALPHABET_CHARS = {"dna": "GATCURYWSMKHBVDN", "aa": "BZXACDEFGHIKLMNPQRSTVWY*"}
AMBIGUITY_CHARS = {"dna": "RYWSMKHBVDN", "aa": "BZX*"}


def check_characters(path, names, rows, dtype, gap=None, missing=None):
    """exit with every illegal character, how often it occurs and where it first occurs; warn when the data look like the
    other datatype"""
    legal = ALPHABET_CHARS[dtype]; legal = legal + legal.lower() + "-?" + (gap or "") + (missing or "")
    ok = np.zeros(256, bool); ok[[ord(c) for c in legal if ord(c) < 256]] = True
    counts = np.zeros(256, np.int64)
    for r in rows: counts += np.bincount(np.frombuffer(r, dtype=np.uint8), minlength=256)
    bad = [int(c) for c in np.flatnonzero((counts > 0) & ~ok)]
    if bad:
        first = {}
        for t, r in zip(names, rows):
            b = np.frombuffer(r, dtype=np.uint8)
            for c in bad:
                if c not in first:
                    pos = np.flatnonzero(b == c)
                    if len(pos): first[c] = (t, int(pos[0]) + 1)
            if len(first) == len(bad): break
        shown = "; ".join(f"{(chr(c) if 32 <= c < 127 else '0x%02x' % c)!r} x{counts[c]:,} (first in {first[c][0]} at site {first[c][1]:,})" for c in bad[:8]) + (f"; and {len(bad) - 8} more" if len(bad) > 8 else "")
        letters = sum(counts[ord(c)] for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"); acgt = sum(counts[ord(c)] for c in "ACGTUNacgtun")
        hint = (" Fewer than half the letters are A, C, G, T or N: is this a protein alignment (--datatype aa)?" if dtype == "dna" and letters and acgt < 0.5 * letters else "")
        sys.exit(f"illegal characters in the alignment {os.path.basename(path)}: {shown}. The {dtype} alphabet is {' '.join(ALPHABET_CHARS[dtype].replace('*', ''))}"
                 f"{' (U is read as T)' if dtype == 'dna' else ' and * for a stop'}, with the ambiguity codes {' '.join(AMBIGUITY_CHARS[dtype].replace('*', ''))} as missing data, '-' for a gap and '?' for missing data"
                 f"{' (this NEXUS file declares gap=' + repr(gap) + ' and missing=' + repr(missing) + ')' if (gap and gap != '-') or (missing and missing != '?') else ''}.{hint}")
    if dtype == "aa":
        acgt = sum(counts[ord(c)] for c in "ACGTUNacgtun"); letters = sum(counts[ord(c)] for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
        if letters and acgt > 0.95 * letters: LOG(f"WARNING: {100.0 * acgt / letters:.1f}% of the letters are A, C, G, T, U or N: this looks like a nucleotide alignment, but --datatype aa was given")


class NexusInfo:
    """what the writer needs from a NEXUS file read by the fast path"""
    def __init__(self, datatype, missing, gap, charsets, taxsets):
        self.datatype, self.missing, self.gap, self.charsets, self.taxsets = datatype, missing, gap, charsets, taxsets


_NEX_PUNCT = "()[]{}\\,;:=*\\'\"`+-<>"
_NEX_VALID = {"dna": "GATCRYWSMKHBVDN", "nucleotide": "GATCRYWSMKHBVDN", "rna": "GAUCRYWSMKHBVDN", "protein": "BZXACDEFGHIKLMNPQRSTVWY*"}


_RE_NEXUS_HEAD = re.compile(rb"\s*#NEXUS[ \t]*\r?\n", re.I)
_RE_BEGIN = re.compile(rb"\s*begin\s+(\w+)\s*;", re.I)
_RE_END = re.compile(rb"(?im)^[ \t]*end(?:block)?[ \t]*;")
_RE_MATRIX = re.compile(rb"(?i)(?<![\w])matrix(?![\w])")


def _safename(name):
    """a taxon or set name as Bio.Nexus writes it: quoted when it holds whitespace or punctuation"""
    safe = name.replace("'", "''")
    return "'" + safe + "'" if set(safe) & set(" \t\n" + _NEX_PUNCT) else safe


def _compact4nexus(orig_list):
    """site list (0-based) -> NEXUS ranges (1-based), exactly as Bio.Nexus writes them"""
    if not orig_list: return ""
    orig_list = sorted(set(orig_list)); shortlist = []
    clist = orig_list[:]; clist.append(clist[-1] + 0.5)
    while len(clist) > 1:
        step = 1
        for i, x in enumerate(clist):
            if x == clist[0] + i * step: continue
            elif i == 1 and len(clist) > 3 and clist[i + 1] - x == x - clist[0]: step = x - clist[0]
            else:
                sub = clist[:i]
                if len(sub) == 1: shortlist.append(str(sub[0] + 1))
                elif step == 1: shortlist.append("%d-%d" % (sub[0] + 1, sub[-1] + 1))
                else: shortlist.append("%d-%d\\%d" % (sub[0] + 1, sub[-1] + 1, step))
                clist = clist[i:]; break
    return " ".join(shortlist)


def read_nexus_fast(path):
    """a NEXUS alignment read straight into byte rows: one DATA or CHARACTERS block (plain names, no comments, no
    quoting, format options limited to datatype / missing / gap / interleave) and at most one SETS block of charset
    and taxset statements with plain site ranges. Anything else raises _Unsupported and Bio.Nexus is used instead.
    Works on the file's bytes with offsets, so the only copies made are the sequence rows themselves.
    Returns names, rows (bytearrays), charsets (name -> sorted 0-based sites), NexusInfo."""
    with open(path, "rb") as fh: data = fh.read()
    for ch in (b"[", b"'", b'"', b"(", b"{"):
        if ch in data: raise _Unsupported("comments, quotes or ambiguity sets")
    m = _RE_NEXUS_HEAD.match(data)
    if not m: raise _Unsupported("header")
    pos = m.end(); datab = setsb = None
    while True:
        m = _RE_BEGIN.match(data, pos)
        if not m:
            if data[pos:].strip(): raise _Unsupported("text between blocks")
            break
        kind = m.group(1).lower().decode("ascii"); body_start = m.end()
        from_ = body_start
        if kind in ("data", "characters"):                             # skip the matrix itself: its closing ';' is a plain find
            mi = _RE_MATRIX.search(data, body_start); semi = data.find(b";", mi.end()) if mi else -1
            if semi >= 0: from_ = semi
        e = _RE_END.search(data, from_)
        if not e: raise _Unsupported("block end")
        span = (body_start, e.start())
        if kind in ("data", "characters"):
            if datab is not None: raise _Unsupported("two data blocks")
            datab = span
        elif kind == "sets":
            if setsb is not None: raise _Unsupported("two sets blocks")
            setsb = data[span[0]:span[1]]
        else: raise _Unsupported("block " + kind)
        pos = e.end()
    if datab is None: raise _Unsupported("no data block")
    # the data block: dimensions, format, matrix
    mi = _RE_MATRIX.search(data, datab[0], datab[1])
    if not mi: raise _Unsupported("no matrix")
    header = data[datab[0]:mi.start()].decode("ascii"); mstart = mi.end()
    cmds = [c.strip() for c in header.split(";") if c.strip()]
    ntax = nchar = None; opts = {}
    for c in cmds:
        w = c.split(None, 1); key = w[0].lower(); rest = w[1] if len(w) > 1 else ""
        if key == "dimensions":
            for k, v in re.findall(r"(\w+)\s*=\s*(\d+)", rest):
                if k.lower() == "ntax": ntax = int(v)
                elif k.lower() == "nchar": nchar = int(v)
                else: raise _Unsupported("dimensions " + k)
            if len(re.findall(r"\w+\s*=", rest)) != len(rest.split()): raise _Unsupported("dimensions")
        elif key == "format":
            for tok in re.findall(r"\w+(?:\s*=\s*\S+)?", rest):
                kv = [t.strip() for t in tok.split("=", 1)]
                k = kv[0].lower(); v = kv[1] if len(kv) > 1 else None
                if k in ("datatype", "missing", "gap") and v is not None: opts[k] = v
                elif k == "interleave": opts[k] = (v or "yes").lower()
                else: raise _Unsupported("format option " + k)
        else: raise _Unsupported("command " + key)
    if ntax is None or nchar is None: raise _Unsupported("dimensions")
    dt = opts.get("datatype", "").lower()
    if dt not in _NEX_VALID: raise _Unsupported("datatype")
    missing = opts.get("missing", "?")[0] if "missing" in opts else None
    gap = opts.get("gap", "-")[0] if "gap" in opts else None
    interleave = opts.get("interleave", "no") in ("yes", "")
    end = data.find(b";", mstart, datab[1])
    if end < 0: raise _Unsupported("matrix end")
    if data[end + 1:datab[1]].strip(): raise _Unsupported("text after the matrix")
    names = []; rows = {}
    def lines_between(a_, b_):                                     # non-empty lines of data[a_:b_], one at a time
        p_ = a_
        while p_ < b_:
            q = data.find(b"\n", p_, b_); q = b_ if q < 0 else q
            ln = data[p_:q].strip()
            if ln: yield ln
            p_ = q + 1
    if interleave:
        for ln in lines_between(mstart, end):
            w = ln.split(None, 1)
            if len(w) != 2: raise _Unsupported("interleaved line")
            nm = w[0].decode("ascii")
            if nm not in rows: names.append(nm); rows[nm] = bytearray()
            rows[nm] += b"".join(w[1].split())
    else:
        pending = None
        for ln in lines_between(mstart, end):
            if pending is None:
                w = ln.split(None, 1); nm = w[0].decode("ascii")
                if nm in rows: raise _Unsupported("duplicate taxon")
                names.append(nm); rows[nm] = bytearray(b"".join(w[1].split())) if len(w) > 1 else bytearray()
                pending = nm if len(rows[nm]) < nchar else None
            else:
                rows[pending] += b"".join(ln.split())
                if len(rows[pending]) >= nchar: pending = None
    if len(names) != ntax or any(len(rows[t]) != nchar for t in names):
        if not interleave:                                             # one taxon per line, no continuation: an unaligned matrix reads cleanly this way
            n2, r2 = [], {}
            for ln in lines_between(mstart, end):
                w = ln.split(None, 1); nm = w[0].decode("ascii", "replace")
                if nm in r2: break
                n2.append(nm); r2[nm] = b"".join(w[1].split()) if len(w) > 1 else b""
            if len(n2) == ntax:
                lens = collections.Counter(len(r2[t]) for t in n2); L2 = lens.most_common(1)[0][0]; off = [t for t in n2 if len(r2[t]) != nchar]
                if off and len(off) == ntax and len(lens) == 1:
                    sys.exit(f"NEXUS {os.path.basename(path)}: every sequence has {L2:,} characters but the dimensions line says nchar={nchar:,}: fix nchar")
                if off and (L2 == nchar or lens[L2] >= 2):
                    t = off[0]; sys.exit(f"NEXUS {os.path.basename(path)}: taxon {t} has {len(r2[t]):,} characters but nchar is {nchar:,} ({len(off)} of {ntax} taxa differ): the sequences are not aligned "
                                         f"(every sequence must have nchar characters), or the matrix is interleaved without 'interleave' in its format line")
        raise _Unsupported("matrix dimensions")
    del data
    check_characters(path, names, [rows[t] for t in names], "aa" if dt == "protein" else "dna", gap, missing)
    # the sets block
    charsets = {}; taxsets = {}
    for c in (setsb.decode("ascii") if setsb else "").split(";"):
        c = c.strip()
        if not c: continue
        m = re.match(r"(?is)(charset|taxset)\s+(\S+)\s*=\s*(.*)$", c)
        if not m: raise _Unsupported("sets statement")
        kind, nm, body = m.group(1).lower(), m.group(2), m.group(3).split()
        if kind == "charset":
            sites = []
            for tok in body:
                r = re.match(r"^(\d+)(?:-(\d+))?(?:\\(\d+))?$", tok)
                if not r: raise _Unsupported("charset range " + tok)
                a_, b_, st = int(r.group(1)), int(r.group(2) or r.group(1)), int(r.group(3) or 1)
                if a_ < 1 or b_ < a_ or b_ > nchar: raise _Unsupported("charset range")
                sites.extend(range(a_ - 1, b_, st))
            if nm in charsets: raise _Unsupported("charset twice")
            charsets[nm] = sorted(set(sites))
        else:
            if any(t not in rows for t in body): raise _Unsupported("taxset member")
            taxsets[nm] = sorted(set(body))                           # as Bio.Nexus stores them
    return names, [rows[t] for t in names], charsets, NexusInfo(dt, missing, gap, charsets, taxsets)


def _check_rows(path, names, rows):
    """every sequence ASCII and all of one length, or the message says which is not"""
    if not rows: sys.exit(f"no sequences in {path}")
    for t, r in zip(names, rows):
        b = np.frombuffer(r, dtype=np.uint8); hi = np.flatnonzero(b >= 128)
        if len(hi):
            ch = bytes(b[hi[0]:hi[0] + 4]).decode("utf-8", "replace")[0]
            sys.exit(f"non-ASCII character {ch!r} in the sequence of taxon {t} at site {int(hi[0]) + 1:,} of {os.path.basename(path)}: sequences may hold only ASCII letters, '-' and '?'")
    lens = collections.Counter(len(r) for r in rows); L = lens.most_common(1)[0][0]
    bad = [(t, len(r)) for t, r in zip(names, rows) if len(r) != L]
    if bad:
        shown = ", ".join(f"{t} has {l:,}" for t, l in bad[:5]) + (f", and {len(bad) - 5} more" if len(bad) > 5 else "")
        sys.exit(f"the sequences in {os.path.basename(path)} are not all the same length: {shown}, while {lens[L]} of the {len(rows)} sequences have {L:,} sites. "
                 "The input must be an alignment: every sequence padded with gaps to the same length")


def read_input(path, partitions):
    """returns names, rows (one bytearray per taxon), loci, format, datatype declared by a NEXUS file (or None), and
    the NEXUS record the writer needs (a Bio.Nexus object, or NexusInfo from the fast reader)"""
    try:
        with open(path, "rb") as fh: raw = fh.read(64)
    except OSError as e:
        sys.exit(f"cannot read the alignment {path}: {e.strerror or e}")
    if not raw.strip(): sys.exit(f"the alignment {path} is empty")
    try: head = raw.decode("utf-8").lstrip()[:7].upper()
    except UnicodeDecodeError: sys.exit(f"the alignment {path} is not a text file (compressed? binary?): give a plain NEXUS, FASTA or PHYLIP file")
    nexus = None; declared = None
    if head.startswith("#NEXUS"):
        try:
            names, rows, charsets, nexus = read_nexus_fast(path); dt = nexus.datatype
        except _Unsupported:
            try: from Bio.Nexus import Nexus
            except ImportError: sys.exit(f"{path} is a NEXUS file the fast reader does not handle (comments, quoted names, duplicates); parsing it needs Biopython, which is not installed: install biopython, or supply plain NEXUS or FASTA/PHYLIP plus a partition file")
            try:
                nexus = Nexus.Nexus(path)
            except Exception as e:
                msg = str(e) or type(e).__name__
                if isinstance(e, StopIteration) or "does not match data length" in msg: msg += " (a sequence shorter or longer than nchar: the sequences are not aligned?)"
                sys.exit(f"Bio.Nexus could not parse {path}:\n  {msg}\n"
                         "This tool reads NEXUS with Biopython and accepts only what Bio.Nexus accepts. Common causes: "
                         "commas inside CHARSET site lists (NEXUS wants whitespace: 'CHARSET x = 1-10 21-30;'), "
                         "taxon names with '-' or other punctuation that are not quoted, or unusual block formatting. "
                         "Check for commas in the SETS block, and/or supply the alignment as FASTA with a RAxML-style "
                         "partition file (--partitions).")
            names = list(nexus.taxlabels); rows = [bytearray(str(nexus.matrix[t]).encode()) for t in names]
            charsets = nexus.charsets; dt = (nexus.datatype or "").lower()
        fmt = "nexus"
        if charsets and partitions:
            sys.exit("the NEXUS file has charsets and a partition file was given: use one or the other, not both")
        if charsets:
            loci = [(k, np.array(sorted(v), dtype=int)) for k, v in charsets.items()]
            owner = np.full(len(rows[0]), -1, dtype=np.int64)
            for j, (k, arr) in enumerate(loci):
                dup = arr[owner[arr] >= 0]
                if len(dup):
                    first = int(dup[0]); others = sorted({loci[int(owner[x])][0] for x in dup})
                    sys.exit(f"site {first + 1} is in charset {k} and in charset {others[0]}: a site may be in one charset only "
                             f"({len(dup)} site{'s' if len(dup) > 1 else ''} of {k} {'are' if len(dup) > 1 else 'is'} already in {', '.join(others)})")
                owner[arr] = j
        elif partitions:
            loci = parse_raxml_partitions(partitions, len(rows[0]))
        else:
            sys.exit("the NEXUS file has no charsets: add a SETS block or give --partitions")
        if dt in ("dna", "rna", "nucleotide"): declared = "dna"
        elif dt == "protein": declared = "aa"
    else:
        if not partitions: sys.exit("a FASTA or PHYLIP alignment needs --partitions")
        try:
            with open(path, encoding="utf-8") as fh: txt = fh.read()
        except UnicodeDecodeError: sys.exit(f"the alignment {path} is not a UTF-8 text file (compressed? binary?)")
        if head.startswith(">"):
            names, seqs = read_fasta(txt); fmt = "fasta"
            empty = [t for t, sq in zip(names, seqs) if not sq]
            if empty: sys.exit(f"FASTA: {len(empty)} record{'s' if len(empty) > 1 else ''} with no sequence under the header in {path}: {', '.join(empty[:5])}{' ...' if len(empty) > 5 else ''}")
        else:
            if not head[:1].isdigit(): sys.exit(f"cannot tell the format of {path}: a NEXUS file starts with #NEXUS, a FASTA file with '>', a PHYLIP file with its counts")
            names, seqs = read_phylip(txt); fmt = "phylip"
        rows = [bytearray(s.encode()) for s in seqs]; del seqs, txt
        _check_rows(path, names, rows)
        loci = parse_raxml_partitions(partitions, len(rows[0]))
    _check_rows(path, names, rows)
    dup = sorted(t for t, c in collections.Counter(names).items() if c > 1)
    if dup: sys.exit(f"duplicate taxon names: {dup[:5]}{' ...' if len(dup) > 5 else ''}")
    if not loci: sys.exit("no loci defined")
    return names, rows, loci, fmt, declared, nexus


# ------------------------------------------------------------------ output
def write_masked(path, out_dir, fmt, names, rows, nexus):
    stem, ext = os.path.splitext(os.path.basename(path)); out = os.path.join(out_dir, f"{stem}.masked{ext or '.' + fmt}")
    if fmt == "nexus" and isinstance(nexus, NexusInfo):
        with open(out, "wb") as fh:
            fh.write(b"#NEXUS\nbegin data;\n" + f"dimensions ntax={len(names)} nchar={len(rows[0])};\n".encode())
            fh.write(("format datatype=" + nexus.datatype + (" missing=" + nexus.missing if nexus.missing else "")
                      + (" gap=" + nexus.gap if nexus.gap else "") + ";\nmatrix\n").encode())
            pad = max(len(_safename(t)) for t in names)
            for t, r in zip(names, rows): fh.write(_safename(t).ljust(pad + 1).encode() + r + b"\n")
            fh.write(b";\nend;\n")
            if nexus.charsets or nexus.taxsets:
                lines = ["\nbegin sets"] + [f"charset {_safename(k)} = {_compact4nexus(v)}" for k, v in nexus.charsets.items() if v] + \
                        [f"taxset {_safename(k)} = {' '.join(_safename(t) for t in v)}" for k, v in nexus.taxsets.items() if v] + ["end;\n"]
                fh.write(";\n".join(lines).encode())
    elif fmt == "nexus":
        try: from Bio.Seq import Seq
        except ImportError: sys.exit("writing the masked alignment back into this NEXUS file needs Biopython, which is not installed")
        for t, r in zip(names, rows): nexus.matrix[t] = Seq(r.decode())
        nexus.write_nexus_data(filename=out, interleave=False, append_sets=True)
    elif fmt == "fasta":
        with open(out, "wb") as fh:
            for t, r in zip(names, rows): fh.write(f">{t}\n".encode() + r + b"\n")
    else:
        with open(out, "wb") as fh:
            fh.write(f"{len(names)} {len(rows[0])}\n".encode())
            for t, r in zip(names, rows): fh.write(f"{t} ".encode() + r + b"\n")
    return out


# ------------------------------------------------------------------ shared machinery
def wmedian_stack(vals, wts):
    """weighted median of vals over axis 0, nan-aware; (L, P) -> (P,). The work is done with the loci as
    the last (contiguous) axis, one row per pair, which sorts and sums far faster than along axis 0."""
    ok = np.isfinite(vals) & (wts > 0)
    v = np.where(ok, vals, np.inf).T.copy(); w = np.where(ok, wts, 0.0).T.copy()      # (P, L)
    order = np.argsort(v, axis=1)
    v_s = np.take_along_axis(v, order, 1); w_s = np.take_along_axis(w, order, 1)
    cw = np.cumsum(w_s, axis=1, dtype=np.float64); tot = cw[:, -1]
    idx = np.minimum((cw < tot[:, None] / 2).sum(1), v.shape[1] - 1)
    out = v_s[np.arange(len(idx)), idx]
    return np.where(tot > 0, out, np.nan)


_LANCZOS = np.array([0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
                     12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7])
def lgamma_v(x):
    """log Gamma for arrays of positive x: Lanczos (g = 7, 9 terms), arguments below 1 shifted up by one; within
    about 3e-15 relative of math.lgamma"""
    x = np.asarray(x, dtype=np.float64); small = x < 1.0
    xs = np.where(small, x + 1.0, x) - 1.0; t = xs + 7.5
    s = _LANCZOS[0] + np.sum(_LANCZOS[1:] / (xs[..., None] + np.arange(1, 9)), axis=-1)
    out = 0.5 * math.log(2 * math.pi) + (xs + 0.5) * np.log(t) - t + np.log(s)
    return np.where(small, out - np.log(x), out)
_LF = np.zeros(1)
_LOG_FLOOR = math.log(1e-300)
def log_fact(upto):
    """log k! for k = 0..upto, cached and grown on demand"""
    global _LF
    if len(_LF) <= upto: _LF = np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, upto + 1, dtype=float)))])
    return _LF
def _block_terms(k0, cnt, n, d, a_, b_, cpair, binom, B):
    """log pmf of the mismatch counts k0[i] .. k0[i] + cnt[i] - 1 for each pair i, as a (pairs, B) array with -inf
    beyond cnt. Binomial coefficients from the log-factorial table; the beta-binomial's two gamma terms from one
    log-gamma per pair at the block's first term and the recurrence lgamma(x + 1) = lgamma(x) + log x along the
    block, which costs a log per term instead of a log-gamma."""
    nint = n.astype(np.int64); lf = log_fact(int(n.max()) + 1)
    j = np.arange(B); kk = k0[:, None] + j[None, :]; valid = j[None, :] < cnt[:, None]
    ksafe = np.where(valid, kk, 0); nk = np.where(valid, nint[:, None] - kk, 0)
    lchoose = lf[nint][:, None] - lf[ksafe] - lf[nk]
    kf = ksafe.astype(float); nf_ = n[:, None]
    if binom: lp = lchoose + kf * np.log(d)[:, None] + (nf_ - kf) * np.log(1 - d)[:, None]
    else:
        k0f = k0.astype(float)
        base_a = lgamma_v(k0f + a_)                                            # lgamma(k0 + a), then + log(k0 + a + t)
        inc_a = np.log(k0f[:, None] + a_[:, None] + j[None, :])
        lg_a = base_a[:, None] + np.cumsum(inc_a, axis=1) - inc_a
        x = n - k0f + b_                                                        # lgamma(n - k0 + b), then - log(x - t)
        base_b = lgamma_v(x)
        inc_b = np.log(np.maximum(x[:, None] - (j[None, :] + 1), 1e-300))
        cs_b = np.cumsum(inc_b, axis=1)
        lg_b = base_b[:, None] - np.concatenate([np.zeros((len(k0), 1)), cs_b[:, :-1]], axis=1)
        lp = lchoose + lg_a + lg_b + cpair[:, None]
        small = n < 2                                                           # pairs with fewer than two sites stay binomial
        if small.any():
            lbin = lchoose + kf * np.log(d)[:, None] + (nf_ - kf) * np.log(1 - d)[:, None]
            lp = np.where(small[:, None], lbin, lp)
    return np.where(valid, lp, -np.inf)


def _logsumexp_rows(lp):
    mx = np.max(lp, axis=1); mx = np.where(np.isfinite(mx), mx, 0.0)
    return mx + np.log(np.sum(np.exp(lp - mx[:, None]), axis=1))


def bb_tails(o, n, E, phi, upper):
    """log tail probability of the mismatch count among n shared sites: binomial with mean E, or, when the locus
    is overdispersed, beta-binomial with mean E and variance E + phi E^2. upper -> log P(X >= o), else
    log P(X <= o). All pairs at once, in blocks of 64 terms per pair. Each tail is summed from the observed count
    away from the mean (upward for the upper tail, downward for the lower) until a block adds less than 1e-20 of the
    running sum while the pmf is falling (or the end of the range is reached), so far tails come out to full
    precision instead of as 1 - P(the other side). A zero expectation gives log 1 for the lower tail and the floor
    for an upper tail above zero. Where the overdispersion is negligible (Beta concentration a + b above 1e9) the
    binomial limit is used: there the beta-binomial's large gamma differences would cancel in float64 and leave
    nats of rounding error. Results are floored at log(1e-300) and capped at log(1). phi may be one value or one per
    pair, so the pairs of many loci can be computed in one call."""
    o = np.asarray(o, float); n = np.asarray(n, float); E = np.asarray(E, float); m = len(o)
    if m == 0: return np.empty(0)
    phi = np.broadcast_to(np.asarray(phi, float), (m,))
    zero = E <= 1e-12
    if zero.any():
        out = np.empty(m)
        out[zero] = np.where(o[zero] <= 0, 0.0, _LOG_FLOOR) if upper else 0.0
        if (~zero).any(): out[~zero] = bb_tails(o[~zero], n[~zero], E[~zero], phi[~zero], upper)
        return out
    E = np.minimum(E, n * (1 - 1e-9)); d = E / n
    bbm = phi >= 1e-8
    if bbm.any() and not bbm.all():                                   # binomial and beta-binomial pairs in one call: each its own way
        out = np.empty(m)
        out[~bbm] = bb_tails(o[~bbm], n[~bbm], E[~bbm], 0.0, upper); out[bbm] = bb_tails(o[bbm], n[bbm], E[bbm], phi[bbm], upper)
        return out
    a_ = b_ = cpair = None; binom = not bbm.any()
    if not binom:
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = np.where(n >= 2, np.minimum((d + phi * n * d) / ((1 - d) * np.maximum(n - 1, 1)), 1 - 1e-6), 0.5)
        a_ = d * (1 / rho - 1); b_ = (1 - d) * (1 / rho - 1)
        lim = (a_ + b_) > 1e9
        if lim.any():
            out = np.empty(m)
            out[lim] = bb_tails(o[lim], n[lim], E[lim], 0.0, upper)
            if (~lim).any(): out[~lim] = bb_tails(o[~lim], n[~lim], E[~lim], phi[~lim], upper)
            return out
        cpair = -lgamma_v(n + a_ + b_) - (lgamma_v(a_) + lgamma_v(b_) - lgamma_v(a_ + b_))
    oi = o.astype(np.int64); nint = n.astype(np.int64); B = 64
    out = np.full(m, -np.inf)
    if upper: active = np.flatnonzero(oi <= nint); kstart = oi.copy(); last = nint
    else: active = np.flatnonzero(oi >= 0); top_k = np.minimum(oi, nint)
    while len(active):
        if upper:
            cnt = np.minimum(B, last[active] - kstart[active] + 1); k0 = kstart[active]
        else:
            cnt = np.minimum(B, top_k[active] + 1); k0 = top_k[active] - cnt + 1
        has = cnt > 0; active, cnt, k0 = active[has], cnt[has], k0[has]
        if not len(active): break
        lp = _block_terms(k0, cnt, n[active], d[active], None if a_ is None else a_[active],
                          None if b_ is None else b_[active], None if cpair is None else cpair[active], binom, B)
        blk = _logsumexp_rows(lp)
        out[active] = np.logaddexp(out[active], blk)
        rows = np.arange(len(active))
        if upper:
            kstart[active] += cnt
            done = kstart[active] > last[active]
            falling = lp[rows, cnt - 1] < lp[rows, 0]                          # moving up, the pmf falls
        else:
            top_k[active] -= cnt
            done = top_k[active] < 0
            falling = lp[rows, 0] < lp[rows, cnt - 1]                          # moving down, the pmf falls
        done |= falling & (blk < out[active] - 46.0)
        active = active[~done]
    return np.maximum(np.minimum(out, 0.0), _LOG_FLOOR)


def bb_logpmf(k, n, E, phi):
    """log P(X = k) under the same model as bb_tails (beta-binomial with mean E and variance E + phi E^2, the binomial
    limit where the overdispersion is negligible or phi is 0); used to screen pairs: P(X <= k) >= P(X = k)"""
    k = np.asarray(k, float); n = np.asarray(n, float); E = np.asarray(E, float)
    E = np.minimum(np.maximum(E, 1e-300), n * (1 - 1e-9)); d = E / n
    lchoose = lgamma_v(n + 1) - lgamma_v(k + 1) - lgamma_v(n - k + 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        lbin = lchoose + k * np.log(d) + (n - k) * np.log1p(-d)
        if phi < 1e-8: return lbin
        rho = np.where(n >= 2, np.minimum((d + phi * n * d) / ((1 - d) * np.maximum(n - 1, 1)), 1 - 1e-6), 0.5)
        a_ = d * (1 / rho - 1); b_ = (1 - d) * (1 / rho - 1)
        lim = (a_ + b_) > 1e9
        a2 = np.where(lim, 1.0, a_); b2 = np.where(lim, 1.0, b_)
        lbb = (lchoose + lgamma_v(k + a2) + lgamma_v(n - k + b2) - lgamma_v(n + a2 + b2)
               - lgamma_v(a2) - lgamma_v(b2) + lgamma_v(a2 + b2))
    return np.where(lim, lbin, lbb)


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


def _expand(tri, n, iu):
    """an upper-triangle vector (diagonal included, row-major) as the full symmetric n x n matrix"""
    M = np.empty((n, n), tri.dtype); M[iu] = tri; M[iu[1], iu[0]] = tri; return M


def _cpus():
    """CPUs this process may run on (its affinity mask on Linux), else the machine's count"""
    try: return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError): return max(1, os.cpu_count() or 1)


def _mem_available():
    """bytes of memory the run can take without swapping: Linux MemAvailable, else the free physical pages, else 2 GB"""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"): return int(line.split()[1]) * 1024
    except OSError: pass
    try: return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (ValueError, OSError, AttributeError): return 2 << 30


MEM_BUDGET = None                 # bytes the run may use (--memory; main sets it, half of what is available by default); None: no budget


def _rss():
    """bytes resident now (Linux: /proc; elsewhere the peak so far)"""
    try:
        with open("/proc/self/statm") as fh: return int(fh.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError, IndexError): pass
    try:
        import resource; r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(r if sys.platform == "darwin" else r * 1024)
    except Exception: return 0


def _mem_room():
    """bytes a stage may take for its threads and chunks: the budget less what is resident now (never below a sixteenth
    of the budget, so a stage can always run small), or what the machine has available when there is no budget"""
    if MEM_BUDGET is None: return _mem_available()
    return max(MEM_BUDGET - _rss(), MEM_BUDGET // 16)


def fmem(b):
    """bytes as GB or MB"""
    return f"{b / float(1 << 30):.1f} GB" if b >= (1 << 30) else f"{b / float(1 << 20):.0f} MB"


def memory_estimate(n, L, loci):
    """bytes the run holds resident, from the alignment's shape: the residue codes (taxa x sites), the per-locus pair
    counts (two arrays of loci x pairs), and the taxa x taxa matrices of the pair stage and the test (about 72 bytes per
    taxon pair, counting the pair stage's provisional copies and the pair index). The stages' own temporaries are sized
    to the room left and are not counted. Returns the total and the three parts."""
    nloc = len(loci); ntri = n * (n + 1) // 2; longest = max((len(st) for _, st in loci), default=0)
    x = n * L; sw = 2 * nloc * ntri * (2 if longest < 65536 else 4); nn = 72 * n * n + 16 * ntri
    return x + sw + nn, x, sw, nn



class _Log:
    """stderr prose (a file too, when one is opened); per-locus detail lines reach the terminal only when no progress bar is drawn"""
    def __init__(self): self.real = sys.stderr; self.fh = None; self.bar = False
    def open(self, path): self.fh = open(path, "w", encoding="utf-8")
    def __call__(self, msg):
        if self.bar: self.real.write("\r\x1b[K")
        self.real.write(msg + "\n"); self.real.flush()
        if self.fh: self.fh.write(msg + "\n"); self.fh.flush()
    def detail(self, msg):
        if self.fh: self.fh.write(msg + "\n")
        if not self.bar: self.real.write(msg + "\n")
    def close(self):
        if self.fh: self.fh.close(); self.fh = None
LOG = _Log()


class Progress:
    """stage progress on stderr. `bar` draws one updating line per stage; `auto` is a bar on a terminal and, elsewhere (a
    log file), one progress<TAB>stage<TAB>done<TAB>total<TAB>elapsed_s<TAB>extra line at 2% steps; `none` is silent.
    Ticks may come from worker threads."""
    def __init__(self, mode="none", t0=None):
        self.mode = ("bar" if sys.stderr.isatty() else "sparse") if mode == "auto" else mode
        self.t0 = time.time() if t0 is None else t0; self.lock = threading.Lock(); self.stage = None; self.out = sys.stderr
        self.durations = {}
        enc = (getattr(self.out, "encoding", "") or "").lower().replace("-", "")
        self.fill, self.empty = ("\u2588", "\u2591") if enc == "utf8" else ("#", ".")
    def start(self, stage, total, extra=""):
        with self.lock:
            if self.stage is not None: self._finish()
            self.stage, self.total, self.done, self.extra = stage, max(int(total), 1), 0, extra
            self.ts = time.time(); self.last = 0.0; self.next_frac = 0.02
            self._emit(self.ts)
    def tick(self, k=1, extra=None):
        if self.stage is None: return
        with self.lock:
            self.done += k
            if extra is not None: self.extra = extra
            now = time.time()
            if self.mode == "sparse":
                if self.done >= self.total or self.done / self.total >= self.next_frac:
                    self.next_frac = math.floor(self.done / self.total * 50 + 1) / 50; self._emit(now)
            elif now - self.last >= 0.5 or self.done >= self.total: self._emit(now)
    def finish(self):
        with self.lock: self._finish()
    def _finish(self):
        if self.stage is None: return
        self.done = self.total; now = time.time(); self.durations[self.stage] = round(now - self.ts, 2); self._emit(now, final=True); self.stage = None
    def _emit(self, now, final=False):
        if self.mode == "none": return
        el = now - self.t0; st = now - self.ts
        if self.mode == "sparse":
            self.out.write(f"progress\t{self.stage}\t{self.done}\t{self.total}\t{el:.1f}\t{self.extra}\n"); self.out.flush(); self.last = now; return
        frac = self.done / self.total; w = 28; filled = int(round(frac * w))
        eta = (st / frac - st) if 0 < frac < 1 else 0.0
        left = "" if final else (f"{eta:4.0f} s left" if frac > 0 else "")
        line = f"\r{self.stage:<12s} {self.fill * filled}{self.empty * (w - filled)} {self.done:>8,}/{self.total:<8,} {st:5.0f} s {left:>11s}  {self.extra}"
        self.out.write(line[:220] + "\x1b[K" + ("\n" if final else "")); self.out.flush(); self.last = now
PROG = Progress("none")


def _git_commit():
    """the short commit of the script's checkout, when it is one; "+dirty" when the script differs from it"""
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        c = subprocess.run(["git", "-C", here, "rev-parse", "--short", "HEAD"], capture_output=True, text=True, timeout=5).stdout.strip() or "unknown"
        if c != "unknown" and subprocess.run(["git", "-C", here, "status", "--porcelain", "--", os.path.abspath(__file__)], capture_output=True, text=True, timeout=5).stdout.strip():
            c += "+dirty"
        return c
    except Exception:
        return "unknown"


def _ranges(sites):
    """a locus's site indices as [start, end, step] runs (0-based, inclusive): one run for a contiguous locus"""
    sites = np.asarray(sites); runs = []; i = 0
    while i < len(sites):
        j = i + 1; step = int(sites[j] - sites[i]) if j < len(sites) else 1
        while j + 1 < len(sites) and sites[j + 1] - sites[j] == step: j += 1
        if j == i + 1 and step != 1 and j < len(sites):           # two sites with an odd gap: emit the first alone
            runs.append([int(sites[i]), int(sites[i]), 1]); i += 1; continue
        if j >= len(sites): j = len(sites) - 1
        runs.append([int(sites[i]), int(sites[j]), step if j > i else 1]); i = j + 1
    return runs


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""): h.update(blk)
    return h.hexdigest()



def _nj(D):
    """neighbour joining on a distance matrix (numpy, n x n): an unrooted tree as an adjacency dict
    node -> [(neighbour, length)], tips numbered 0..n-1; negative branch lengths are clamped to zero"""
    n = len(D); adj = {i: [] for i in range(n)}; nxt = n
    dist = {(i, j): float(D[i, j]) for i in range(n) for j in range(i + 1, n)}
    d = lambda i, j: dist[(i, j)] if i < j else dist[(j, i)]
    active = list(range(n))
    while len(active) > 2:
        m = len(active); r = {i: sum(d(i, j) for j in active if j != i) for i in active}
        best = None
        for a in range(m):
            for b in range(a + 1, m):
                i, j = active[a], active[b]; q = (m - 2) * d(i, j) - r[i] - r[j]
                if best is None or q < best[0]: best = (q, i, j)
        _, i, j = best; dij = d(i, j)
        li = 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2)); lj = dij - li
        li, lj = max(li, 0.0), max(lj, 0.0)
        u = nxt; nxt += 1; adj[u] = [(i, li), (j, lj)]; adj[i].append((u, li)); adj[j].append((u, lj))
        for k in active:
            if k not in (i, j): dist[(min(u, k), max(u, k))] = max(0.5 * (d(i, k) + d(j, k) - dij), 0.0)
        active = [k for k in active if k not in (i, j)] + [u]
    if len(active) == 2:
        i, j = active; l = max(d(i, j), 0.0); adj[i].append((j, l)); adj[j].append((i, l))
    return adj


def _root_on_edge(adj, u, v, f):
    """a new root node on edge (u, v), the fraction f of the edge's length from u"""
    l = next(ln for y, ln in adj[u] if y == v); r = max(adj) + 1
    adj[u] = [(y, ln) for y, ln in adj[u] if y != v] + [(r, l * f)]
    adj[v] = [(y, ln) for y, ln in adj[v] if y != u] + [(r, l * (1 - f))]
    adj[r] = [(u, l * f), (v, l * (1 - f))]
    return r


def _midpoint_edge(adj, ntips):
    """the edge and fraction where the midpoint of the longest tip-to-tip path lies: (u, v, f)"""
    best = None
    for a in range(ntips):                                       # paths from every tip: n is tiny
        stack = [(a, None, 0.0, [a])]
        while stack:
            x, parent, dist_, path = stack.pop()
            if x < ntips and x != a and (best is None or dist_ > best[0]): best = (dist_, path)
            for y, ln in adj[x]:
                if y != parent: stack.append((y, x, dist_ + ln, path + [y]))
    total, path = best; half = total / 2.0; cum = 0.0
    for u, v in zip(path, path[1:]):
        l = next(ln for y, ln in adj[u] if y == v)
        if cum + l >= half or (u, v) == (path[-2], path[-1]): return u, v, 0.0 if l == 0 else min(max((half - cum) / l, 0.0), 1.0)
        cum += l
    return path[0], path[1], 0.5


def _newick(adj, root, labels):
    def rec(x, parent):
        kids = [(y, ln) for y, ln in adj[x] if y != parent]
        if not kids: return labels[x]
        return "(" + ",".join(f"{rec(y, x)}:{ln:.3f}" for y, ln in kids) + ")"
    return rec(root, None) + ";"


def cell_trees(D_exp, D_obs, labels, og=None):
    """the face-to-face trees of one cell: neighbour joining on the expected and on the observed distances (both
    F81-corrected with the locus plateau, in percent), each rooted at the node where the outgroup tip `og` attaches
    (midpoint-rooted when there is none). Returns (expected newick, observed newick)"""
    n = len(labels)
    if n < 2: return None
    out = []
    for D in (D_exp, D_obs):
        t = _nj(D)
        if og is None: u, v, f = _midpoint_edge(t, n)
        else: u, v, f = og, t[og][0][0], 1.0
        r = _root_on_edge(t, u, v, f); out.append(_newick(t, r, labels))
    return tuple(out)

def write_report(a, names, loci, rep, dtype, fmt, K, K_eff, M_eff, t0):
    """report.npz, a ZIP of .npy members and JSON that np.load opens like an .npz, and report.html, the viewer with
    the bundle inline as base64: one file whatever its size"""
    n, nloc = len(names), len(loci)
    dense = rep["dense"]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED, compresslevel=4) as zf:          # level 4: half the time of 6 for 5% more bytes
        def put(name, arr):
            b = io.BytesIO(); np.lib.format.write_array(b, np.ascontiguousarray(arr), version=(1, 0)); zf.writestr(name, b.getvalue())
        def put_json(name, obj, compact=False): zf.writestr(name, json.dumps(obj, separators=(",", ":")) if compact else json.dumps(obj, indent=1))
        for k in ("n_comp", "n_fail", "n_similar", "side", "focal_sites"): put(f"cells/{k}.npy", rep[k])
        numcols = ("med_identity", "med_expected", "excess_m", "closest", "closest_overlap", "closest_identity", "closest_expected",
                   "best_match", "best_match_identity", "partner", "div_log10p", "sim_log10p", "locus_closest", "locus_closest_log10p")
        if not dense: put("cells/sparse_locus.npy", rep["sparse_locus"]); put("cells/sparse_taxon.npy", rep["sparse_taxon"])
        for k in numcols: put(f"cells/{k}.npy", rep[k])
        for k in ("comparators", "comp_div_log10p", "comp_sim_log10p"): put(f"cells/{k}.npy", rep[k])
        for li, (cells_, blk) in rep["codes"].items():
            put(f"codes/{li}.npy", blk); put(f"codes/{li}_cells.npy", cells_)
        for li, tr in rep["trees"].items(): put_json(f"trees/{li}.json", tr, compact=True)         # machine-read only: no indentation
        for li, ev in rep["similar"].items(): put_json(f"similar/{li}.json", ev, compact=True)
        for li, (rows_, blk) in rep["seqs"].items():
            put(f"seq/{li}.npy", blk); put(f"seq/{li}_taxa.npy", rows_)
        put("taxa/closest.npy", rep["taxa_closest"])
        side = rep["side"]; n_comp = rep["n_comp"]; n_fail = rep["n_fail"]
        tested = n_comp > 0
        counts = dict(taxa=n, loci=nloc, cells=int(tested.sum()), present=int((rep["focal_sites"] > 0).sum()),
                      masked=int(((side & 3) > 0).sum()), masked_divergent=int((side & 1 > 0).sum()), masked_similar=int((side & 2 > 0).sum()),
                      masked_both=int(((side & 3) == 3).sum()), ambiguous=int((side & 4 > 0).sum()),
                      similar_tested=int((side & 64 > 0).sum()), similar_pairs=int(sum(sum(1 for e in v if e["fail"]) for v in rep["similar"].values())),
                      similar_sources=int((side & 128 > 0).sum()), similar_uncertain=int((side & 8 > 0).sum()),
                      near_miss={str(j): int(((n_fail == j) & tested & ((side & 1) == 0)).sum()) for j in range(1, M_eff)},
                      excluded_taxa=int(len(rep["thin"])), skipped_loci=int(sum(1 for st in rep["loc_status"] if st != "tested")),
                      untested_sequences=int(rep["untested"]["sequences"]), untested=rep["untested"], uncovered_sites=int(rep.get("uncovered_sites", 0)))
        thin = set(int(t) for t in rep["thin"])
        put_json("taxa.json", dict(
            names=list(names),
            status=[("excluded: %d loci with at least %d residues, %d needed" % (int(rep["present"][i]), a.fit_overlap, MIN_LOCI_PER_TAXON)) if i in thin else "tested" for i in range(n)],
            loci_present=[int(v) for v in (rep["focal_sites"] > 0).sum(1)]))
        put_json("loci.json", dict(
            names=[nm for nm, _ in loci], length=[int(len(st)) for _, st in loci], ranges=[_ranges(st) for _, st in loci],
            status=list(rep["loc_status"]), rate=[None if not np.isfinite(v) else round(float(v), 5) for v in rep["rate"]],
            plateau=[round(float(v), 3) for v in rep["plateaus"]], phi=[None if not np.isfinite(v) else float(v) for v in rep["phi"]],
            test_rate=[None if not np.isfinite(v) else round(float(v), 5) for v in rep["test_rate"]], similar_tests=[int(v) for v in rep["pair_tests"]]))
        put_json("manifest.json", dict(
            format_version=REPORT_FORMAT, tool="LocusMasker", version=__version__, commit=_git_commit(),
            python=platform.python_version(), numpy=np.__version__,
            system=platform.platform(), date=time.strftime("%Y-%m-%d %H:%M:%S"), command=" ".join(sys.argv),
            options=dict(datatype=dtype, n_comp=K, K_eff=K_eff, M_eff=M_eff, fit_overlap=a.fit_overlap, alpha=a.alpha, donor_threshold=a.donor_threshold,
                         mask_similar_only=a.mask_similar_only, mask_different_only=a.mask_different_only, threads=a.threads, all_tests=a.all_tests,
                         report_codes=a.report_codes, report_sequences=rep["seq_policy"]),
            input=dict(file=os.path.basename(a.aln), bytes=os.path.getsize(a.aln), sha256=_sha256(a.aln), format=fmt,
                       partitions=os.path.basename(a.partitions) if a.partitions else None),
            dataset=a.name or (os.path.basename(os.path.dirname(os.path.abspath(a.aln))) if os.path.splitext(os.path.basename(a.aln))[0].lower() == "alignment" and os.path.basename(os.path.dirname(os.path.abspath(a.aln))) not in ("", ".", "work") else os.path.splitext(os.path.basename(a.aln))[0]),
            counts=counts, test=dict(tests=int(rep["tests"]["total"]), divergent_tests=int(rep["tests"]["divergent"]),
                                     similar_comparator_tests=int(rep["tests"]["similar_comparator"]), similar_closest_tests=int(rep["tests"]["similar_closest"]),
                                     similar_candidates=int(rep["tests"]["locus_sets"]), comparisons=int(rep["tests"]["divergent"]),
                                     bonferroni_k=round(float(rep["kq"]), 4), log_bar=round(float(rep["logbar"]), 4), log10_bar=round(float(rep["logbar"]) / LN10, 4),
                                     phi_dataset=round(float(rep["phi_bar"]), 6), half=HALF),
            carried=dict(numeric_columns="dense" if dense else "sparse", codes=a.report_codes, sequences=rep["seq_policy"], trees="masked",
                         tree_partners="the pair's partner", similar="every too-similar test and every source, per locus (similar/<locus>.json)",
                         tests="every cell's comparators and both of each comparator's tests (cells/comp_div_log10p, cells/comp_sim_log10p, log10 p by comparator rank), the closest match at the locus and its too-similar test (cells/locus_closest, cells/locus_closest_log10p)",
                         best_match_cells="all" if n <= BEST_MATCH_ALL_MAX_TAXA else "fails"),
            code_map=dict(mismatch_count_max=K_eff, similar_base=CODE_SIMILAR, tick=CODE_TICK, absent=CODE_ABSENT),
            side_bits=dict(masked_divergent=1, masked_similar=2, ambiguous=4, similar_uncertain=8, flagged_not_masked=16, similar_not_masked=32, similar_tested=64, similar_source=128),
            timings=dict(PROG.durations, total=round(time.time() - t0, 2))))
    data = buf.getvalue()
    with open(f"{a.out}/report.npz", "wb") as fh: fh.write(data)
    tpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viewer"); tpl = os.path.join(tpl_dir, "viewer.html")
    if not os.path.exists(tpl):
        LOG(f"report.html not written: the viewer template is not beside the script ({tpl}); report.npz holds the data"); return
    def text(name):
        with open(os.path.join(tpl_dir, name), encoding="utf-8") as fh: return fh.read()
    html = text("viewer.html")
    def inline_script(mt): return "<script>\n" + text(mt.group(1)) + "\n</script>"
    def inline_style(mt): return "<style>\n" + text(mt.group(1)) + "\n</style>"
    def inline_icon(mt):                                           # the logo and favicons, as data: URIs
        src = os.path.join(tpl_dir, mt.group(2)); kind = "image/svg+xml" if src.endswith(".svg") else "image/png"
        with open(src, "rb") as fh: return f'{mt.group(1)}="data:{kind};base64,' + base64.b64encode(fh.read()).decode("ascii") + '"'
    html = html.replace("@@BUNDLE@@", base64.b64encode(data).decode("ascii"))   # the marker first, so nothing inlined afterwards can carry it
    html = re.sub(r'<script src="([^"]+)"></script>', inline_script, html)
    html = re.sub(r'<link rel="stylesheet" href="([^"]+)">', inline_style, html)
    html = re.sub(r'(href|src)="(icons/[^"]+)"', inline_icon, html)
    with open(f"{a.out}/report.html", "w", encoding="utf-8") as fh: fh.write(html)
    LOG(f"report: report.html ({os.path.getsize(f'{a.out}/report.html') >> 20} MB), report.npz ({len(data) >> 20} MB)")


def pair_divergences(S, W, n, iu, rr=None, ok_r=None, threads=1, progress=None, stage="pairs"):
    """for every taxon pair (flattened n*n): the site-weighted median over loci of the pair's divergence
    (rate-normalised by rr when given), the same median of the resolution floor (half a mismatch per locus),
    and the number of loci contributing; pairs with fewer than MIN_LOCI_FOR_M loci get NaN. Self pairs carry
    no vote; loci without a rate (ok_r False) carry none either. Computed in chunks of pairs straight from
    the per-locus counts, so the peak memory is a few chunk-sized temporaries rather than several copies of
    the whole loci x pairs matrix; the per-pair arithmetic is exactly that of pair_median."""
    nloc, npair = S.shape                                               # the pairs are the upper triangle, diagonal included
    m = np.full(npair, np.nan); m_fl = np.full(npair, np.nan); nl = np.zeros(npair, dtype=np.int64)
    diag = np.flatnonzero(iu[0] == iu[1])
    # chunk size from the machine: each chunk's temporaries are about 64 bytes per locus x pair, the in-flight chunks
    # (one per thread) may take a quarter of the memory available now, and no chunk is below 8 MB (the per-chunk
    # Python overhead must stay negligible against the sorts) or above 256 MB (nothing is gained beyond it)
    threads = max(1, min(threads, npair // PAIRS_PER_THREAD))                   # no more threads than the pairs warrant
    per_chunk = min(256 << 20, max(32 << 20, _mem_room() // (4 * threads)))
    chunk = max(int(per_chunk // (64 * max(nloc, 1))), 1)
    if threads > 1: chunk = max(min(chunk, -(-npair // (4 * threads))), int((32 << 20) // (64 * max(nloc, 1))), 1)   # four chunks per thread, but never below 32 MB
    chunk = min(chunk, npair)
    threads = max(1, min(threads, -(-npair // chunk)))                       # no more workers than chunks
    LOG(f"pair stage: {-(-npair // chunk)} chunk{'s' if npair > chunk else ''} of {chunk} pairs, about {max(chunk * 64 * nloc >> 20, 1)} MB each, {threads} thread{'s' if threads > 1 else ''}")
    def do_chunk(c0):
        c1 = min(c0 + chunk, npair)
        Sc = S[:, c0:c1].astype(np.float32); Wc = W[:, c0:c1].astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            Dv = np.where(Wc > 0, 100.0 * (1.0 - Sc / np.maximum(Wc, 1)), np.nan)
        Wfl = 100.0 * HALF / np.maximum(Wc, 1.0)
        Wm = Wc.copy()
        d = diag[(diag >= c0) & (diag < c1)] - c0
        Wm[:, d] = 0.0
        if rr is not None:
            Dv = Dv / rr[:, None]; Wfl = Wfl / rr[:, None]
            Wm = np.where(ok_r[:, None], Wm, 0.0)
        nan_fb = np.full(c1 - c0, np.nan)
        m[c0:c1], nl[c0:c1] = pair_median(Dv, Wm, MIN_LOCI_FOR_M, nan_fb)
        m_fl[c0:c1], _ = pair_median(Wfl, Wm, MIN_LOCI_FOR_M, nan_fb)
        if progress: progress.tick()
    starts = range(0, npair, chunk)
    if progress: progress.start(stage, len(starts))
    if threads > 1:
        with ThreadPoolExecutor(threads) as ex: list(ex.map(do_chunk, starts))
    else:
        for c0 in starts: do_chunk(c0)
    if progress: progress.finish()
    return _expand(m, n, iu), _expand(m_fl, n, iu), _expand(nl, n, iu)


def flatten_pairs(idx, top, valid, kk):
    gi = np.repeat(idx, kk); ci = top[idx].ravel(); vi = valid[idx].ravel(); rk = np.tile(np.arange(kk), len(idx))
    return gi[vi], ci[vi], rk[vi]


def norm_quantile(p):
    """inverse standard normal CDF, via bisection on erfc (no scipy)"""
    if p < 0.5: return -norm_quantile(1 - p)
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


def prepare(a, names, rows, loci, dtype, K_eff, t0):
    """everything the test needs before it runs: the alignment as residue codes, the per-locus shared and same counts
    of every pair, the locus plateaus and rates, and the genome-wide pair divergences m (with the resolution floor
    m_hi and the ranking key mj). a needs fit_overlap and threads."""
    n, L = len(names), len(rows[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = np.empty((n, L), dtype=np.int8)                               # one row at a time: no joined copy of the alignment
    letters = np.zeros(256, np.int64)
    for i, r in enumerate(rows): b = np.frombuffer(r, dtype=np.uint8); X[i] = lut[b]; letters += np.bincount(b, minlength=256)
    # letters outside the alphabet are missing data; when they are most of the letters the datatype (or the file) is wrong
    is_letter = np.zeros(256, bool); is_letter[65:91] = True; is_letter[97:123] = True
    n_letters = int(letters[is_letter].sum()); n_in = int(letters[is_letter & (lut >= 0)].sum())
    if n_letters and n_in < 0.95 * n_letters:
        LOG(f"WARNING: {n_letters - n_in:,} of the {n_letters:,} letters ({100.0 * (n_letters - n_in) / n_letters:.1f}%) are ambiguity codes and count as missing data")

    # ---- pass 1: per-locus shared/same counts; alignment-wide identity falls out of the sums
    nloc = len(loci)
    # the same-residue and shared-site counts of every pair, per locus: symmetric, so the upper triangle (diagonal
    # included, row-major) is kept and a locus's full matrix is rebuilt when the test needs it
    iu = np.triu_indices(n); ntri = len(iu[0]); diag_pos = np.flatnonzero(iu[0] == iu[1])
    tri_idx = np.empty((n, n), np.int32); tri_idx[iu] = np.arange(ntri); tri_idx[iu[1], iu[0]] = np.arange(ntri)
    cdt = np.uint16 if max(len(sites) for _, sites in loci) < 65536 else np.uint32   # counts never exceed a locus's length
    S = np.empty((nloc, ntri), cdt); W = np.empty((nloc, ntri), cdt)
    plateaus = np.empty(nloc)
    def locus_counts(li):
        sites = loci[li][1]
        Xm = X[:, sites]; Pm = (Xm >= 0).astype(np.float32)
        sh = Pm @ Pm.T; sm = np.zeros((n, n), np.float32)
        # each letter's product runs over the sites where the letter occurs at all: the other columns are zero and
        # add nothing, and a site holds only a few of the letters (3 of 20 in a typical protein alignment)
        if len(alphabet) > 4:
            occ = np.zeros((len(alphabet) + 1, Xm.shape[1]), bool); occ[Xm + 1, np.arange(Xm.shape[1])] = True   # row 0: missing
            where = [np.flatnonzero(occ[i + 1]) for i in range(len(alphabet))]
        else:
            where = [np.flatnonzero((Xm == i).any(0)) for i in range(len(alphabet))]
        for i, cols in enumerate(where):
            if len(cols):
                Ai = (Xm[:, cols] == i).astype(np.float32); sm += Ai @ Ai.T
        S[li] = sm[iu]; W[li] = sh[iu]                                  # exact: float32 counts below 2^24
        cnt = np.bincount(Xm[Xm >= 0], minlength=len(alphabet)).astype(float)
        f = cnt / max(cnt.sum(), 1.0)
        plateaus[li] = max(100.0 * (1.0 - float((f ** 2).sum())), 1.0)
        PROG.tick()
    # loci in flight at once: the threads, but no more than a quarter of the memory available now allows, each
    # locus needing about three taxa x taxa float32 matrices plus two one-hot copies of its sites
    per_locus = 3 * n * n * 4 + 2 * n * max(len(sites) for _, sites in loci) * 4
    # the matrix-product work of a typical locus: the presence product plus one product per letter over the sites where
    # the letter occurs, so the mean number of distinct letters per site (from a sample of loci) sets the width
    sample = [loci[i][1] for i in np.linspace(0, nloc - 1, min(nloc, 64)).astype(int)]
    kbar = float(np.mean([np.mean([(X[:, st] == i).any(0).mean() for i in range(len(alphabet))]) * len(alphabet) for st in sample]))
    flops = 2.0 * n * n * float(np.mean([len(sites) for _, sites in loci])) * (1.0 + kbar)
    workers = max(1, min(a.threads, nloc, int((_mem_room() // 2) // per_locus), int(round(flops / PASS1_FLOPS_PER_THREAD))))
    PROG.start("counts", nloc)
    if workers > 1:
        with ThreadPoolExecutor(workers) as ex: list(ex.map(locus_counts, range(nloc), chunksize=max(1, nloc // (8 * workers))))
    else:
        for li in range(nloc): locus_counts(li)
    PROG.finish()
    if workers < a.threads: LOG(f"pass one: {workers} of {a.threads} threads (the loci are small, few, or large in memory)")
    focal_sites = np.ascontiguousarray(W[:, diag_pos].T).astype(np.uint32)   # each taxon's residues at each locus, before the exclusion zeroes anything
    present = (W[:, diag_pos] >= a.fit_overlap).sum(0)
    thin = np.flatnonzero(present < MIN_LOCI_PER_TAXON)
    if len(thin):
        drop = np.isin(iu[0], thin) | np.isin(iu[1], thin)
        S[:, drop] = 0; W[:, drop] = 0                                     # out of the test entirely
        LOG(f"{len(thin)} taxa with fewer than {MIN_LOCI_PER_TAXON} loci present are excluded from the test and left unmasked: "
            + ", ".join(names[i] for i in thin[:20]) + (f" and {len(thin) - 20} more (the report lists them)" if len(thin) > 20 else ""))
    tot_same = _expand(S.sum(0, dtype=np.float64), n, iu); tot_shared = _expand(W.sum(0, dtype=np.float64), n, iu)
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(tot_shared > 0, 100.0 * tot_same / tot_shared, np.nan)
    np.fill_diagonal(ani, np.nan)
    ident = (tot_shared > 0) & (tot_same == tot_shared); np.fill_diagonal(ident, False)   # genomes identical at every shared site: six copies of a bad sequence must not vouch for each other
    del tot_same, tot_shared
    cls = np.arange(n)                                                 # identical genomes form a class (its lowest index); one class, one comparator
    while True:
        new = np.minimum(cls, np.where(ident, cls[None, :], n).min(1))
        if np.array_equal(new, cls): break
        cls = new
    has_dup = bool(ident.any())
    if has_dup: LOG(f"{int((cls != np.arange(n)).sum())} genomes are identical to another at every shared site: an identical genome is never a comparator, and no two comparators are identical to each other")
    # deterministic tie-break for comparator ranking: pooled divergence, then taxon order
    tie = np.nan_to_num(100.0 - ani, nan=100.0) * 1e-9 + np.arange(n)[None, :] * 1e-12
    del ani
    LOG(f"per-locus identities done, {time.time() - t0:.0f}s")

    # ---- pair divergence m: capped site-weighted median across loci of rate-normalised divergence
    # ONE rate per locus, estimated once and used everywhere (to normalise loci inside m and as the
    # expectation in the test). It is estimated where the test looks: each taxon's nearest neighbours
    # by the provisional pair divergence m0 (the raw site-weighted median over loci, available before rates exist), one equal
    # vote per neighbour pair resolved beyond its floor, of the multiple-hit-corrected ratio
    # c(d + half a mismatch) / c(m0), m0 the provisional pair divergence (raw site-weighted median).
    m0, m0_fl, _ = pair_divergences(S, W, n, iu, threads=a.threads, progress=PROG, stage="pairs0")   # the floor: the imputed half mismatch per locus, as a percentage
    m0_hi = np.maximum(m0, m0_fl); res0 = m0 > m0_fl
    dist = np.where(np.isfinite(m0) & res0 & ~ident & (cls == np.arange(n))[None, :], m0, np.inf); np.fill_diagonal(dist, np.inf)   # neighbours by the provisional pair divergence (a median over loci, not swayed by one long fast locus); an unresolved pair casts no vote, an identical genome is no neighbour, one genome stands for its copies
    Kn = min(K_eff, n - 1)
    nb = np.argsort(dist, axis=1, kind="stable")[:, :Kn]      # stable: a tie goes to the lower index on every platform (numpy's default sort is unstable and its tie order differs with the CPU's SIMD level)
    gi0 = np.repeat(np.arange(n), Kn); ci0 = nb.ravel()
    okn = np.isfinite(dist[gi0, ci0]) & res0[gi0, ci0]
    gi0, ci0 = gi0[okn], ci0[okn]
    m0n = m0_hi[gi0, ci0]
    rate = np.full(nloc, np.nan)
    t0i = tri_idx[gi0, ci0]
    PROG.start("rates", nloc)
    for li in range(nloc):
        PROG.tick()
        o_, s_ = W[li][t0i].astype(np.float32), S[li][t0i].astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            d_ = np.where(o_ > 0, 100.0 * (1.0 - s_ / np.maximum(o_, 1)), np.nan)
        val = np.isfinite(d_) & (o_ >= a.fit_overlap)
        if val.sum() < MIN_FIT_PAIRS: continue
        P_ = plateaus[li]
        corr = lambda p: -P_ * np.log(1.0 - np.minimum(np.asarray(p, np.float64), 0.98 * P_) / P_)   # float64 under every numpy (1.x kept float32 here when p was float32, 2.x promotes)
        rate[li] = float(np.median(corr(d_[val] + 100.0 * HALF / o_[val]) / corr(m0n[val])))
    PROG.finish()
    del m0, m0_fl, m0_hi, res0, dist, nb, m0n
    lsites = np.array([len(s) for _, s in loci], dtype=float)
    ok_r = np.isfinite(rate) & (rate > 0)
    if ok_r.sum() >= 3:
        scale = wmedian_stack(rate[ok_r, None], lsites[ok_r, None])[0]
        rate = rate / max(scale, 1e-6)
    rr = np.where(ok_r, rate, np.nan)
    m, m_fl, nl = pair_divergences(S, W, n, iu, rr, ok_r, threads=a.threads, progress=PROG, stage="pairs")
    np.fill_diagonal(m, np.nan)
    # resolution floor of the pair: the median over its loci of half a mismatch per locus (rate-normalised);
    # the divergence used on the too-divergent side is the larger of the measured median and this floor
    m_hi = np.maximum(m, m_fl); np.fill_diagonal(m_hi, np.nan)
    m_xs = np.maximum(m - m_fl, 0.0)                                  # divergence resolved beyond the floor: the rate-vote weight
    del m_fl
    keep = np.ones(n, bool); keep[thin] = False
    nlk = nl[np.ix_(keep, keep)][np.triu_indices(int(keep.sum()), 1)]
    fb = float((nlk < MIN_LOCI_FOR_M).mean()) if len(nlk) else 0.0
    del nl, nlk
    LOG(f"pair divergences done, {time.time() - t0:.0f}s"
          + (f"; note: no locus has a rate (no nearest-neighbour pairs resolved beyond their resolution floor), so no pair has a divergence" if ok_r.sum() == 0 else
             f"; note: {fb:.0%} of taxon pairs overlap in <{MIN_LOCI_FOR_M} loci with a rate and are never compared" if fb > 0.05 else ""))

    mj = np.where(np.isfinite(m) & ~ident, m, np.inf) + tie           # the comparator ranking key: m with the tie-break; an identical genome is never a comparator (Rob, 30 Oct 2026)
    del tie
    return argparse.Namespace(X=X, S=S, W=W, iu=iu, diag_pos=diag_pos, tri_idx=tri_idx, m=m, m_hi=m_hi, m_xs=m_xs, mj=mj,
                              plateaus=plateaus, rate=rate, focal_sites=focal_sites, thin=thin, present=present, cls=cls, has_dup=has_dup)


def divergence_matrix(aln, partitions=None, datatype="dna", fit_overlap=20, n_comp=5, threads=1):
    """the genome-wide pair divergences exactly as a run computes them, for scripts that need them (a run writes no
    matrix): returns (names, m, m_hi), m in percent with NaN where a pair overlaps in fewer than 5 loci or a taxon is
    excluded, m_hi = max(m, the pair's resolution floor)"""
    datatype = str(datatype).lower()
    if datatype not in ("dna", "aa"): raise ValueError(f"datatype must be 'dna' or 'aa', not {datatype!r}")
    names, rows, loci, fmt, declared, nexus = read_input(aln, partitions)
    if declared and declared != datatype: raise ValueError(f"datatype {datatype!r} was given but the NEXUS file declares {declared!r}")
    if not isinstance(nexus, NexusInfo): check_characters(aln, names, rows, datatype, getattr(nexus, "gap", None), getattr(nexus, "missing", None))
    n = len(names); K_eff = min(n_comp, n - 1); K_eff -= 1 - K_eff % 2
    P = prepare(argparse.Namespace(fit_overlap=fit_overlap, threads=max(1, min(threads, _cpus()))), names, rows, loci, datatype, K_eff, time.time())
    return names, P.m, P.m_hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True); ap.add_argument("--partitions"); ap.add_argument("--out", required=True)
    ap.add_argument("--datatype", required=True, choices=("dna", "aa"))
    ap.add_argument("--n-comp", type=int, default=5); ap.add_argument("--fit-overlap", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--donor-threshold", type=float, default=2.0, help="the direction of a too-similar pair is called (the copy masked, its source not) when the larger of the "
                    "two members' test-1 excesses, in absolute value, is more than this many times the smaller; otherwise the direction is uncertain and both members are masked (default 2)")
    ap.add_argument("--mask-similar-only", action="store_true"); ap.add_argument("--mask-different-only", action="store_true")
    ap.add_argument("--threads", type=int, default=1, help="at most this many threads or worker processes (0: every CPU this process may use); each stage takes as many as pay off")
    ap.add_argument("--memory", type=float, help="GB of memory the run may use (default: half of what is available when it starts; 0: no budget); the stages size their threads and chunks to stay within it, and a run that cannot fit stops before it starts")
    ap.add_argument("--progress", choices=("auto", "bar", "none"), default="auto", help="auto: a progress bar on a terminal, progress lines at 2%% steps otherwise")
    ap.add_argument("--report-codes", choices=("fails", "masked", "none"), default="fails", help="site codes in the report: for every cell with a failing comparator (default), masked cells only, or none")
    ap.add_argument("--report-sequences", choices=("auto", "all", "fails", "masked", "none"), default="auto",
                    help="sequences in the report: the whole alignment, the rows of every cell with a failing comparator, of masked cells only, or none; auto: all when the alignment is small, else fails")
    ap.add_argument("--no-report", action="store_true", help="write no report.html / report.npz")
    ap.add_argument("--all-tests", action="store_true", help="also write tests.tsv.gz: one row per test the run performed (every comparison, every pair of taxa at every locus, every direction test) with its statistic, exact p-value and result; large")
    ap.add_argument("--name", help="dataset name shown in the report (default: the alignment file's stem, or its folder when the file is called alignment.*)")
    a = ap.parse_args(); K = a.n_comp
    global PROG
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    stale = [f for f in ("tests.tsv.gz",) if not a.all_tests and os.path.exists(os.path.join(a.out, f))] + \
            [f for f in ("report.html", "report.npz") if a.no_report and os.path.exists(os.path.join(a.out, f))] + \
            [f for f in ("cells.tsv", "pairs.npy", "similar.tsv", "per_locus.tsv", "run.log", "report_data.js") if os.path.exists(os.path.join(a.out, f))]
    PROG = Progress(a.progress, t0); LOG.bar = PROG.mode == "bar"
    LOG(f"LocusMasker {__version__} ({_git_commit()}) on {time.strftime('%Y-%m-%d %H:%M:%S')}; python {platform.python_version()}, numpy {np.__version__}, {platform.platform()}")
    LOG("command: " + " ".join(sys.argv))
    if stale: LOG(f"note: {a.out} holds files this run does not write ({', '.join(stale)}): they are from an earlier run")
    ncpu = _cpus(); asked = a.threads; a.threads = ncpu if asked <= 0 else max(1, min(asked, ncpu))
    if asked > ncpu: LOG(f"threads: {asked} asked, {ncpu} CPUs available to this process: using {a.threads}")
    global MEM_BUDGET
    avail = _mem_available(); MEM_BUDGET = None if a.memory is not None and a.memory <= 0 else int(a.memory * (1 << 30)) if a.memory else avail // 2
    if MEM_BUDGET is not None and MEM_BUDGET > avail: LOG(f"WARNING: --memory {a.memory:g} GB is more than the {avail / (1 << 30):.1f} GB available now")
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    if a.mask_similar_only and a.mask_different_only:
        sys.exit("--mask-similar-only and --mask-different-only are mutually exclusive")
    if a.partitions and not os.path.isfile(a.partitions): sys.exit(f"the partition file {a.partitions} does not exist")
    PROG.start("read", 1)
    names, rows, loci, fmt, declared, nexus = read_input(a.aln, a.partitions)
    PROG.finish()
    covered = np.zeros(len(rows[0]), bool)
    for _, st in loci: covered[st] = True
    nuncov = int((~covered).sum())
    if nuncov:
        runs = _ranges(np.flatnonzero(~covered)); shown = ", ".join(f"{r[0] + 1}-{r[1] + 1}" if r[1] > r[0] else str(r[0] + 1) for r in runs[:6])
        LOG(f"WARNING: {nuncov:,} of the {len(covered):,} sites are in no {'charset' if fmt == 'nexus' and not a.partitions else 'partition'} and are ignored by every test and never masked"
            f" (sites {shown}{', ...' if len(runs) > 6 else ''})")
    if declared and declared != a.datatype:
        sys.exit(f"--datatype {a.datatype} was given but the NEXUS file declares {declared}")
    dtype = a.datatype
    if not isinstance(nexus, NexusInfo): check_characters(a.aln, names, rows, dtype, getattr(nexus, "gap", None), getattr(nexus, "missing", None))
    if len(loci) < 10:
        sys.exit(f"only {len(loci)} loci defined: this is a phylogenomic test and needs at least 10")
    n, L = len(names), len(rows[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    if n < 5: sys.exit(f"only {n} taxa: the test needs at least 5")
    K_eff = min(K, n - 1); K_eff -= 1 - K_eff % 2                    # largest odd <= min(K, n-1)
    M_eff = (K_eff + 1) // 2                                          # fails that make a majority
    if K_eff < K: LOG(f"n-comp shrunk to {K_eff} ({n} taxa)")
    LOG(f"{n} taxa, {L} sites, {len(loci)} loci, {dtype}, {fmt}")
    need, mx, msw, mnn = memory_estimate(n, L, loci)
    parts = f"{fmem(msw)} of pair counts for {n:,} taxa x {len(loci):,} loci, {fmem(mnn)} of taxa x taxa matrices, {fmem(mx)} of residue codes"
    if MEM_BUDGET is None: LOG(f"memory: about {fmem(need)} resident ({parts}); no budget (--memory 0)")
    else:
        LOG(f"memory: about {fmem(need)} resident ({parts}); budget {fmem(MEM_BUDGET)}" + ("" if a.memory else " (half of the memory available now; --memory sets it)"))
        if need > MEM_BUDGET:
            sys.exit(f"this run needs about {fmem(need)} of memory ({parts}), more than the budget of {fmem(MEM_BUDGET)}"
                     + ("" if a.memory else ", half of what is available now") + ": --memory raises the budget (--memory 0 removes it), or give fewer taxa or loci")

    P = prepare(a, names, rows, loci, dtype, K_eff, t0)
    X, S, W, iu, diag_pos, tri_idx = P.X, P.S, P.W, P.iu, P.diag_pos, P.tri_idx
    m, m_hi, m_xs, mj, plateaus, rate = P.m, P.m_hi, P.m_xs, P.mj, P.plateaus, P.rate
    run_z.cls, run_z.has_dup = P.cls, P.has_dup
    focal_sites, thin, present = P.focal_sites, P.thin, P.present
    masked = rows                                                      # blanked in place; X is its own copy

    seq_policy = a.report_sequences if a.report_sequences != "auto" else ("all" if n * L <= REPORT_SEQ_ALL_MAX else "fails")
    rep = dict(focal_sites=focal_sites, thin=thin, present=present, seq_policy=seq_policy if not a.no_report else "none",
               codes_policy=a.report_codes if not a.no_report else "none", no_report=a.no_report, uncovered_sites=nuncov)
    nmasked = run_z(a, names, loci, X, S, W, iu, diag_pos, tri_idx, m, m_hi, m_xs, mj, plateaus, rate, masked, K_eff, M_eff, rep)

    PROG.start("write", 2)
    out = write_masked(a.aln, a.out, fmt, names, masked, nexus)
    PROG.tick()
    if not a.no_report: write_report(a, names, loci, rep, dtype, fmt, K, K_eff, M_eff, t0)
    PROG.finish()
    LOG(f"outputs: {out}, {a.out}/masked.csv" + (f", {a.out}/tests.tsv.gz" if a.all_tests else "") + ("" if a.no_report else f", {a.out}/report.html")
        + f"\ndone {time.time() - t0:.0f}s")
    LOG.close()


# ------------------------------------------------------------------ the counting test
LSET_CHUNK = 256        # the fewest taxa per slice when a locus's closest matches are found (LSET_ELEMENTS sets the rest)
LN10 = math.log(10.0)
_TEST_BATCH = None
def _test_batch(lis):
    """a batch of loci in a worker process: the test's per-locus function, inherited at the fork"""
    return _TEST_BATCH(lis)


MASKED_COLUMNS = ["locus", "taxon", "reason", "residues_masked",
                  "divergent_fails", "divergent_comparators", "divergent_observed", "divergent_expected", "divergent_excess_m",
                  "divergent_log10p", "similar_test", "similar_partner", "similar_role", "similar_observed", "similar_expected",
                  "similar_row_factor", "similar_log10p", "similar_relatives", "partner_excess_m", "comparators", "log10_bar"]


def fnum(v, nd):
    """a number for the text outputs: fixed decimals, no negative zero, empty when not finite"""
    v = float(v)
    return "" if not math.isfinite(v) else f"{round(v, nd) + 0.0:.{nd}f}"
TESTS_COLUMNS = ["test", "locus", "taxon", "comparator", "overlap", "mismatches", "expected", "row_factor", "log10p", "fail"]


def run_z(a, names, loci, X, S, W, iu, diag_pos, tri_idx, m, m_hi, m_xs, mj, plateaus, rate, masked, K_eff, M_eff, rep):
    n = len(names); nloc = len(loci); ntri = len(iu[0])
    mask_hi = not a.mask_similar_only
    mask_lo = not a.mask_different_only
    views = [np.frombuffer(mb, dtype=np.uint8) for mb in masked]        # writable views for blanking

    # each taxon's candidates in order of pair divergence (the tie-break makes every row a strict order), computed
    # once: a locus's comparators are the first K candidates that share enough sites with the taxon at that locus
    big = n >= 1000
    prepool = ThreadPoolExecutor(min(a.threads, nloc)) if a.threads > 1 and big else None
    if a.threads > 1: LOG(f"comparators: {min(a.threads, nloc) if prepool else 1} thread{'s' if prepool and min(a.threads, nloc) > 1 else ''} ({n:,} taxa; threads pay off from 1,000)")
    order = np.empty((n, n), np.int32)
    def sort_rows(b0):
        b1 = min(b0 + 256, n); order[b0:b1] = np.argsort(mj[b0:b1], axis=1, kind="stable")
    if prepool: list(prepool.map(sort_rows, range(0, n, 256)))
    else:
        for b0 in range(0, n, 256): sort_rows(b0)

    def choose_comparators(li):
        """for every taxon at this locus: its first K_eff candidates (in order of m) that have the overlap the taxon
        needs here, walking the sorted candidates in windows that grow until every taxon has K or runs out"""
        Wl = W[li]; nres = Wl[diag_pos].astype(int); need = np.minimum(a.fit_overlap, nres)
        top = np.full((n, K_eff), -1, np.int64); cnt = np.zeros(n, np.int64)
        rows = np.flatnonzero(nres > 0); start = 0; width = 4 * K_eff
        if getattr(run_z, "has_dup", False):                               # identical genomes in the data: one comparator per class of copies, so the walk is per taxon
            cls = run_z.cls
            for g in rows.tolist():
                seen = set(); k = 0; tg = tri_idx[g]
                for cnd in order[g].tolist():
                    if not np.isfinite(mj[g, cnd]): break
                    if Wl[tg[cnd]] < need[g] or cls[cnd] in seen: continue
                    seen.add(cls[cnd]); top[g, k] = cnd; k += 1
                    if k == K_eff: break
                cnt[g] = k
            rows = rows[:0]
        while len(rows) and start < n:
            stop = min(start + width, n)
            cand = order[rows, start:stop]
            mjc = mj[rows[:, None], cand]
            ok = np.isfinite(mjc) & (cand != rows[:, None]) & (Wl[tri_idx[rows[:, None], cand]] >= need[rows][:, None])
            pos = cnt[rows][:, None] + np.cumsum(ok, axis=1) - ok                 # the comparator slot each accepted candidate takes
            r_, j_ = np.nonzero(ok & (pos < K_eff))
            top[rows[r_], pos[r_, j_]] = cand[r_, j_]
            cnt[rows] = np.minimum(cnt[rows] + ok.sum(1), K_eff)
            exhausted = ~np.isfinite(mjc[:, -1])                                  # non-finite candidates sort last: nothing beyond
            rows = rows[(cnt[rows] < K_eff) & ~exhausted]
            start = stop; width *= 4
        # the closest match AT THE LOCUS: the taxon in the test with the fewest mismatches per shared site (half a mismatch
        # added; ties by m) among those sharing >= --fit-overlap sites with the sequence. When it is not one of the K_eff
        # closest genomes (a tie for closest that includes one counts as in-set) the sequence is tested against it (test 2).
        # Taxa x taxa slices are taken LSET_CHUNK rows at a time.
        Sl = S[li]; closest = np.full(n, -1, np.int64); fork = np.zeros(n, bool)
        rows = np.flatnonzero(cnt > 0)
        lset = max(LSET_CHUNK, LSET_ELEMENTS // max(n, 1))
        for b0 in range(0, len(rows), lset):
            rr = rows[b0:b0 + lset]; t_ = tri_idx[rr]
            ov = Wl[t_].astype(np.float64); sm = Sl[t_].astype(np.float64); mjr = mj[rr]
            okc = (ov >= a.fit_overlap) & np.isfinite(mjr)                       # mj is non-finite on the diagonal, for excluded taxa and for identical genomes
            with np.errstate(invalid="ignore", divide="ignore"):
                d = np.where(okc, (ov - sm + HALF) / np.maximum(ov, 1.0), np.inf)
            dmin = d.min(1); has = np.isfinite(dmin)
            tops = np.zeros((len(rr), n), bool); r_, c_ = np.nonzero(top[rr] >= 0); tops[r_, top[rr][r_, c_]] = True
            inset = ((d <= dmin[:, None]) & tops).any(1)
            key = np.where(okc, d + 1e-9 * np.where(np.isfinite(mjr), mjr, 0.0) / 100.0, np.inf)   # ties by m, well below the count resolution
            closest[rr] = np.where(has, key.argmin(1), -1)
            fork[rr] = has & ~inset
        PROG.tick()
        return top, top >= 0, cnt, nres, closest, fork

    PROG.start("comparators", nloc)
    sel = list(prepool.map(choose_comparators, range(nloc))) if prepool else [choose_comparators(li) for li in range(nloc)]
    if prepool: prepool.shutdown()
    PROG.finish()
    Npairs = 0                                                        # comparisons over all loci: sizes the thread decision only
    for sel_ in sel: Npairs += int(sel_[1][sel_[2] > 0].sum())

    def locus_pairs(li):
        """the comparator pairs of a locus and the locus's rate; None when the locus is skipped"""
        Wl, Sl = W[li], S[li]
        top, valid, ncomp, nres, closest, fork = sel[li]
        idx = np.flatnonzero(ncomp > 0)
        if len(idx) == 0: return "skipped-no-cells"
        gi, ci, rk = flatten_pairs(idx, top, valid, K_eff)
        t = tri_idx[gi, ci]; o = Wl[t].astype(np.float32); same = Sl[t].astype(np.float32); mis = o - same
        ident = 100.0 * same / o; dob = 100.0 - ident
        mp = m[gi, ci]; mph = m_hi[gi, ci]; mpx = m_xs[gi, ci]
        fit = (o >= a.fit_overlap) & (mpx > 0)
        if fit.sum() < MIN_FIT_PAIRS and not np.isfinite(rate[li]): return f"skipped-few-resolved-pairs({int(fit.sum())})"
        P_ = plateaus[li]
        corr = lambda p: -P_ * np.log(1.0 - np.minimum(np.asarray(p, float), 0.98 * P_) / P_)
        pc = 100.0 * HALF / o
        if fit.sum() >= MIN_FIT_PAIRS: rate_c = max(float(np.median((corr(dob + pc) / corr(mph))[fit])), 1e-3)
        else: rate_c = max(float(rate[li]), 1e-3)            # the comparator pairs cannot set a rate (identical relatives): the provisional rate, from the nearest resolved pairs
        def expect(mv):
            mvc = np.maximum(np.asarray(mv, float), 1e-4)
            return np.maximum(P_ * (1.0 - np.exp(-rate_c * corr(mvc) / P_)), 1e-4)
        return dict(idx=idx, gi=gi, ci=ci, rk=rk, o=o, same=same, mis=mis, ident=ident, dob=dob, mp=mp, mph=mph, fit=fit,
                    Wl=Wl, Sl=Sl, ncomp=ncomp, nres=nres, expect=expect, corr=corr, rate_c=rate_c)

    # worker threads pay off in the test only when a locus's arrays are large (thousands of comparisons, or the pairs of
    # a thousand taxa and more): with small loci the per-call interpreter-lock traffic of many small numpy operations
    # makes threads slower than one
    ntested = sum(1 for t in sel if (t[2] > 0).any())
    pool = ThreadPoolExecutor(min(a.threads, nloc)) if a.threads > 1 and (Npairs / max(ntested, 1) >= TEST_THREAD_PAIRS or big) else None
    # dispersion pre-pass: the robust spread of observed/expected among each locus's resolved comparator pairs with
    # enough expected mismatches, with the Poisson part removed: Var(O/E) = 1/E + phi. The same pass counts the tests
    # each locus will run (its distinct comparisons and its testable pairs), so the Bonferroni bar is known before any
    # test runs. The comparator arrays are kept for the test when they fit in a quarter of the memory available (about
    # 120 bytes per comparison).
    nrows_all = sum(int(t[2].sum()) for t in sel)
    keep_pairs = nrows_all * 120 < _mem_room() // 2
    cache = {}
    def locus_phi(li):
        L = locus_pairs(li); PROG.tick()
        if keep_pairs: cache[li] = L
        if isinstance(L, str): return np.nan, 0.0, 0, 0
        gmin = np.minimum(L["gi"], L["ci"]).astype(np.int64); gmax = np.maximum(L["gi"], L["ci"]).astype(np.int64)
        ndist = len(np.unique(gmin * n + gmax))
        fk = np.flatnonzero(sel[li][5]); Y_ = sel[li][4][fk]; top_ = sel[li][0][fk]
        okr = (top_ >= 0) & (top_ != Y_[:, None]) & (L["Wl"][tri_idx[np.where(top_ >= 0, top_, 0), Y_[:, None]]] >= a.fit_overlap)
        npt = int((okr.sum(1) >= 3).sum())
        E = L["o"] * L["expect"](L["mph"]) / 100.0
        use = L["fit"] & (E >= PHI_MIN_E)
        if use.sum() < 10: return np.nan, 0.0, ndist, npt
        r = np.log((L["mis"][use] + HALF) / E[use])
        s_ = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        return max(s_ ** 2 - float(np.median(1.0 / E[use])), 0.0), float(use.sum()), ndist, npt
    phi_hat = np.full(nloc, np.nan); phi_n = np.zeros(nloc); ndist_l = np.zeros(nloc, np.int64); npt_l = np.zeros(nloc, np.int64)
    PROG.start("dispersion", nloc)
    for li, (h, c, nd_, np_) in enumerate(pool.map(locus_phi, range(nloc)) if pool else map(locus_phi, range(nloc))):
        phi_hat[li] = h; phi_n[li] = c; ndist_l[li] = nd_; npt_l[li] = np_
    PROG.finish()
    if int(ndist_l.sum()) == 0 and int(npt_l.sum()) == 0:
        # no test can run: say why, in the words of the thresholds, and stop before writing outputs that would mean nothing
        reasons = collections.Counter(locus_pairs(li) if isinstance(cache.get(li, None), str) or not keep_pairs else cache[li] for li in range(nloc))
        reasons = collections.Counter(v if isinstance(v, str) else "tested" for v in reasons.elements())
        thin = rep["thin"]; keep = np.ones(n, bool); keep[thin] = False
        nk = int(keep.sum()); fin = int(np.isfinite(m[np.ix_(keep, keep)]).sum() // 2); npairs = nk * (nk - 1) // 2; nrated = int(np.isfinite(rate).sum())
        why = [f"{len(thin)} of {n} taxa excluded (fewer than {MIN_LOCI_PER_TAXON} loci with at least {a.fit_overlap} residues)"]
        if nrated == 0: why.append(f"no locus has a rate: at no locus are {MIN_FIT_PAIRS} nearest-neighbour pairs resolved beyond their resolution floor (the sequences are too similar, or share too few sites, for their divergence to show), so no pair has a divergence m")
        elif nk >= 2: why.append(f"{npairs - fin:,} of the {npairs:,} pairs of the remaining taxa overlap in fewer than {MIN_LOCI_FOR_M} loci with a rate and are never compared")
        for st, k in sorted(reasons.items()):
            if st == "skipped-no-cells": why.append(f"{k} loci with no sequence that has a comparator (no taxon pair shares enough sites, or the taxa are excluded)")
            elif st.startswith("skipped-few-resolved-pairs"): why.append(f"{k} loci with fewer than {MIN_FIT_PAIRS} comparator pairs resolved beyond their resolution floor (the sequences are too similar, or too short, to set a locus rate)")
        sys.exit(f"nothing could be tested: no sequence has a comparator at any locus. {'; '.join(why)}. "
                 f"Check --datatype, the partitions and the missing data; --fit-overlap {a.fit_overlap} sets the sites a pair must share.")
    okp = phi_n > 0
    phi_bar = wmedian1(phi_hat[okp], phi_n[okp]) if okp.any() else 0.0
    LOG(f"dispersion: dataset-wide phi {phi_bar:.4f} (robust SD of O/E beyond Poisson {math.sqrt(phi_bar):.2f}); "
          f"per-locus phi 10/50/90%: {np.nanpercentile(phi_hat[okp], [10, 50, 90]).round(4).tolist() if okp.any() else 'n/a'}")
    # the Bonferroni family is every test the run performs: one too-divergent and one too-similar test per distinct
    # comparator pair, and one too-similar test per sequence tested against its closest match
    D_tests = int(ndist_l.sum()); S_tests = int(npt_l.sum()); N_tests = 2 * D_tests + S_tests
    logbar = math.log(a.alpha / max(N_tests, 1))
    kq = norm_quantile(1 - a.alpha / max(N_tests, 1))
    lbar10 = logbar / LN10
    nfork = int(sum(int(sel_[5].sum()) for sel_ in sel))
    LOG(f"tests: {D_tests:,} too-divergent comparisons + {D_tests:,} too-similar comparisons + {S_tests:,} closest-match tests "
        f"({nfork:,} sequences whose closest match at the locus is not among their comparators, {S_tests:,} of them with three relatives "
        f"to set the expectation) = N {N_tests:,}; every test is held to alpha / N = {a.alpha / max(N_tests, 1):.3g} (log10 {lbar10:.2f})")

    # ---- the report's per-cell columns: dense taxa x loci while the grid is small enough, else rows for cells with a fail
    dense = n * nloc <= DENSE_CELLS_MAX
    NUMCOLS = (("med_identity", np.float32, np.nan), ("med_expected", np.float32, np.nan), ("excess_m", np.float32, np.nan),
               ("closest", np.int32, -1), ("closest_overlap", np.uint32, 0), ("closest_identity", np.float32, np.nan), ("closest_expected", np.float32, np.nan),
               ("best_match", np.int32, -1), ("best_match_identity", np.float32, np.nan), ("partner", np.int32, -1),
               ("div_log10p", np.float32, np.nan), ("sim_log10p", np.float32, np.nan), ("locus_closest", np.int32, -1), ("locus_closest_log10p", np.float32, np.nan))
    for k in ("n_comp", "n_fail", "n_similar", "side"): rep[k] = np.zeros((n, nloc), np.uint8)
    for k, dt, fill in NUMCOLS: rep[k] = np.full((n, nloc), fill, dt) if dense else []
    rep["comparators"] = np.full((n, nloc, K_eff), -1, np.int32) if dense else []
    for k in ("comp_div_log10p", "comp_sim_log10p"): rep[k] = np.full((n, nloc, K_eff), np.nan, np.float32) if dense else []
    rep["dense"] = dense; sparse_loc, sparse_tax = [], []
    rep["codes"] = {}; rep["seqs"] = {}; rep["trees"] = {}; rep["similar"] = {}; rep["loc_status"] = ["tested"] * nloc
    rep["phi"] = np.full(nloc, np.nan); rep["test_rate"] = np.full(nloc, np.nan); rep["pair_tests"] = npt_l

    codes_policy, seq_policy = rep["codes_policy"], rep["seq_policy"]
    bm_all = n <= BEST_MATCH_ALL_MAX_TAXA

    tests_fh = None
    if a.all_tests:
        tests_fh = gzip.open(f"{a.out}/tests.tsv.gz", "wt", compresslevel=4, encoding="utf-8", newline="\n")
        tests_fh.write("\t".join(TESTS_COLUMNS) + "\n")
    every = tests_fh is not None

    phi_l = (phi_n * np.where(okp, phi_hat, 0.0) + PHI_N0 * phi_bar) / (phi_n + PHI_N0)   # each locus's dispersion, shrunk toward the dataset's
    def process_locus(li):
        """everything the test decides and writes for one locus, as rows, text blocks and mask actions; touches no
        shared state, so loci can be processed in parallel and written in order. A generator: it yields its requests
        for exact tails (kind, observed, shared sites, expected) and receives the tails; run_batch drives a batch of
        loci so that one call serves them all"""
        name, sites = loci[li]
        L = cache.pop(li, None) if keep_pairs else None
        if L is None: L = locus_pairs(li)
        if isinstance(L, str): return dict(skipped=True, log=f"{name}: {L}", status=L)
        idx, gi, ci, rk, o, same, mis, ident, dob, mp, mph, fit, Wl, Sl, ncomp, nres, expect, corr, rate_c = \
            (L[k] for k in ("idx", "gi", "ci", "rk", "o", "same", "mis", "ident", "dob", "mp", "mph", "fit", "Wl", "Sl", "ncomp", "nres", "expect", "corr", "rate_c"))
        phi = phi_l[li]
        top = sel[li][0]

        # the count test with a quasi-Poisson variance: Var(O) = E + phi E^2, phi the locus's dispersion (shrunk toward
        # the dataset's); phi = 0 is the plain binomial. Exact tails, no normal approximation, no effect-size floor.
        dhat = expect(mph)                                          # too-divergent side: the divergence we cannot exclude
        E = o * dhat / 100.0
        excess = mis - E
        up = excess > 0                                              # test 1's upper tails are asked for these comparisons
        kc = ncomp[idx].astype(np.intp); ncell = len(idx); ar = np.arange(ncell)
        starts = np.r_[0, np.cumsum(kc)[:-1]].astype(np.intp)                   # the platform index type: reduceat, repeat and bincount refuse int64 on 32-bit builds (WebAssembly)
        row = np.repeat(ar, kc); col = np.arange(len(gi)) - np.repeat(starts, kc)
        def cell_matrix(v):                                         # the comparisons laid out one row per cell, NaN beyond its count
            Pm = np.full((ncell, K_eff), np.nan); Pm[row, col] = v; return Pm
        def cell_median(v):
            Pm = cell_matrix(v); Pm.sort(axis=1)                    # NaN sorts last
            return (Pm[ar, (kc - 1) // 2] + Pm[ar, kc // 2]) / 2    # np.median's mean of the middle two
        med_i = cell_median(ident); med_e = cell_median(100 - dhat)
        # the cell's extra divergence beyond expectation, as corrected distance in units of m (a typical locus's
        # divergence): a foreign or paralogous sequence adds a constant here whatever the comparator
        exm = (corr(100 - med_i) - corr(100 - med_e)) / rate_c

        # ---- test 2, too similar, two arms that feed one pair decision.
        # COMPARATOR ARM: every comparison of test 1 again, lower tail: closer to a comparator than the locus rate at
        # their genome-wide divergence predicts? The expectation is evaluated at m itself (a pair cannot be more similar
        # than its loci allow), and held against E and E x f, f <= 1 the row factor: the median over each taxon's OTHER
        # comparators of (observed + 0.5) / (expected + 0.5), the smaller of the pair's two, so a lineage that is slow
        # at the locus (close to everything here) is not a copy of each of its relatives; the half mismatch keeps
        # comparators that expect about none from saying anything. A pair fails only against both expectations.
        # CLOSEST-MATCH ARM: a sequence whose closest match Y at the locus is not one of its comparators is tested
        # against Y, lower tail: is it closer to Y than its own relatives are? The expected mismatch rate is the median,
        # over its comparators that pass test 1 here and share >= --fit-overlap sites with Y (at least three), of their
        # mismatch rate to Y at this locus: no rate model, no m; a conserved locus, a slow lineage or a deep pair's
        # saturation move the relatives as much as the sequence.
        # A failing pair is too similar; the COPY is the member that fits its own relatives worse (the larger test-1
        # excess in absolute value: a foreign sequence in a fast lineage is CLOSER to the lineage's relatives than
        # expected), and it is masked; the other is its source.
        cpos = {int(g): q for q, g in enumerate(idx.tolist())}; Xs = None
        qpos = np.full(n, -1, np.int64); qpos[idx] = ar                     # taxon -> its cell here, -1 when it has none
        pairs = {}                                                   # (min, max) -> the record of the test that fired lowest
        def record(key, rec):
            if key not in pairs or rec["lnp"] < pairs[key]["lnp"]: pairs[key] = rec
        # comparator arm
        E_lo = o * expect(mp) / 100.0
        oe = cell_matrix((mis + HALF) / (E_lo + HALF))
        def med_capped(Q, k):
            Qs = np.sort(Q, axis=1); kk = np.maximum(k, 1)
            v = (Qs[ar, (kk - 1) // 2] + Qs[ar, kk // 2]) / 2
            return np.where(k >= 2, np.minimum(v, 1.0), 1.0)
        f_all = np.ones(n); f_all[idx] = med_capped(oe, kc)                        # the partner is not a comparator
        f_loo = np.ones((n, K_eff))                                                # the partner is comparator r: leave it out
        for r_ in range(K_eff):
            Q = oe.copy(); Q[:, r_] = np.nan
            f_loo[idx, r_] = med_capped(Q, kc - (kc > r_))
        def row_factor(A_, B_):
            out = []
            for x, y in ((A_, B_), (B_, A_)):
                eq = top[x] == y[:, None]; has = eq.any(1)
                out.append(np.where(has, f_loo[x, eq.argmax(1)], f_all[x]))
            return np.minimum(out[0], out[1])
        ff = row_factor(gi, ci)
        dn = mis < E_lo; r_ff = np.flatnonzero(dn & (ff < 1.0))
        # the exact tails: asked for here, computed by the driver for a batch of loci at once, answered in order
        t_up, t_dn, t_ff = yield [("up", mis[up], o[up], E[up]), ("dn", mis[dn], o[dn], E_lo[dn]), ("dn", mis[r_ff], o[r_ff], E_lo[r_ff] * ff[r_ff])]

        # ---- test 1, too divergent: each comparison's upper tail P(X >= observed); a cell fails when a majority of its
        # comparators do (with fewer comparators than requested, a majority of those it has: "ambiguous")
        lp_hi = np.full(len(gi), np.nan)                            # natural-log upper tails where computed
        lp_hi[up] = t_up
        fail_hi = up & (lp_hi < logbar)
        nf = np.add.reduceat(fail_hi.astype(np.intp), starts)
        flag = (kc >= K_eff) & (nf >= M_eff)
        amb = (kc < K_eff) & (nf >= kc // 2 + 1)
        dflag = flag | amb
        lp_in = np.zeros(len(gi))                                   # natural-log lower tails; 0 (p = 1) above expectation
        lp_in[dn] = t_dn
        if len(r_ff): lp_in[r_ff] = np.maximum(lp_in[r_ff], t_ff)
        fail_in = lp_in < logbar
        for j in np.flatnonzero(fail_in).tolist():
            g, c = int(gi[j]), int(ci[j])
            record((min(g, c), max(g, c)), dict(arm="comparator", g=g, y=c, o=float(o[j]), mis=float(mis[j]), E=float(E_lo[j]), f=float(ff[j]),
                                                 lnp=float(lp_in[j]), rel=[]))
        # closest-match arm, every forked cell at once: its comparators' mismatch rates to Y where they pass test 1 here,
        # are not Y and share enough sites with Y; three or more of them set the expectation
        closest, fork = sel[li][4], sel[li][5]
        fk = fork[idx]; Y = np.where(fk, closest[idx], -1)
        tested2 = np.zeros(ncell, bool); o2 = np.zeros(ncell); mis2 = np.zeros(ncell); E2 = np.full(ncell, np.nan)
        lp_lo = np.zeros(ncell); rel2 = [[] for _ in range(ncell)]
        fq = np.flatnonzero(fk)
        if len(fq):
            G = idx[fq]; Yv = Y[fq]; C = top[G]; C0 = np.maximum(C, 0); qc = qpos[C0]
            valid = (C >= 0) & (C != Yv[:, None]) & ~((qc >= 0) & dflag[np.maximum(qc, 0)])
            t = tri_idx[C0, Yv[:, None]]; w = Wl[t].astype(np.float64); valid &= w >= a.fit_overlap
            with np.errstate(invalid="ignore", divide="ignore"):
                rates = np.where(valid, (w - Sl[t].astype(np.float64)) / np.maximum(w, 1.0), np.nan)
            ok3 = valid.sum(1) >= 3
            if ok3.any():
                q3 = fq[ok3]; tgy = tri_idx[G[ok3], Yv[ok3]]
                o2[q3] = Wl[tgy].astype(np.float64); mis2[q3] = o2[q3] - Sl[tgy].astype(np.float64)
                with np.errstate(all="ignore"): E2[q3] = o2[q3] * np.nanmedian(rates[ok3], axis=1)
                tested2[q3] = True
                for i in np.flatnonzero(ok3).tolist(): rel2[fq[i]] = C[i][valid[i]].tolist()
        tq = np.flatnonzero(tested2 & (mis2 < E2))
        (t_lo,) = yield [("dn", mis2[tq], o2[tq], E2[tq])]
        lp_lo[tq] = t_lo
        fail2 = tested2 & (lp_lo < logbar)
        for q in np.flatnonzero(fail2).tolist():
            g = int(idx[q]); y = int(Y[q])
            record((min(g, y), max(g, y)), dict(arm="closest", g=g, y=y, o=o2[q], mis=mis2[q], E=float(E2[q]), f=1.0, lnp=float(lp_lo[q]), rel=list(rel2[q])))
        # the pairs that failed, each once; the copy is the member with the larger |test-1 excess| (median over its
        # comparators of the divergence beyond expectation, in units of m), the source the other
        copy_q = np.zeros(ncell, bool); src_q = np.zeros(ncell, bool); partner_ = np.full(ncell, -1, np.int32); role_ = {}
        npair = np.zeros(ncell, np.intp)
        pair_of = {}                                                 # cell -> (record, partner, is_copy) of its lowest pair
        for (g, y), rec in pairs.items():
            qg, qy = cpos.get(g), cpos.get(y)
            xg = abs(exm[qg]) if qg is not None and np.isfinite(exm[qg]) else -np.inf
            xy = abs(exm[qy]) if qy is not None and np.isfinite(exm[qy]) else -np.inf
            cp, sr = (g, y) if xg >= xy else (y, g)
            # the direction is uncertain, and both members are masked, when the larger excess is not more than --donor-threshold
            # times the smaller (a member outside the test, with no excess, cannot be masked: the other is the copy)
            uncertain = bool(np.isfinite(xg) and np.isfinite(xy) and max(xg, xy) <= a.donor_threshold * min(xg, xy)); rec["uncertain"] = uncertain
            for member, other, is_copy in ((cp, sr, True), (sr, cp, False)):
                qm = cpos.get(member)
                if qm is None: continue
                npair[qm] += 1
                role = "uncertain" if uncertain else "copy" if is_copy else "source"
                if role == "source": src_q[qm] = True
                else: copy_q[qm] = True
                if qm not in pair_of or rec["lnp"] < pair_of[qm][0]["lnp"]: pair_of[qm] = (rec, other, role)
        for qm, (rec, other, role) in pair_of.items(): partner_[qm] = other; role_[qm] = role
        unc_q = np.zeros(ncell, bool)
        for qm, role in role_.items(): unc_q[qm] = role == "uncertain"
        for q in np.flatnonzero(tested2 & (partner_ < 0)).tolist(): partner_[q] = Y[q]
        sflag = copy_q

        # ---- the decision: masked when too divergent (test 1) or the copy of a too-similar pair (test 2)
        masked_q = (dflag & mask_hi) | (sflag & mask_lo)
        fail_q = (nf >= 1) | fail2 | copy_q | src_q | masked_q
        # exact upper tails for every comparison of the cells the outputs describe (and of every comparison under --all-tests),
        # so a cell's decisive comparison (the one that makes or misses the majority) has its exact p-value
        rowsel = np.isnan(lp_hi) & (np.repeat(fail_q, kc) | every)
        (t_rows,) = yield [("up", mis[rowsel], o[rowsel], E[rowsel])]
        lp_hi[rowsel] = t_rows
        need = np.where(kc >= K_eff, M_eff, kc // 2 + 1)
        LP = cell_matrix(lp_hi); LP.sort(axis=1)
        div_lp = np.where(fail_q, LP[ar, np.maximum(need, 1) - 1], np.nan)
        # every comparator's two tests for every cell: log10 p of the too-divergent upper tail and of the too-similar lower
        # tail, by comparator rank as in the report's comparators (NaN: no such comparator). A comparison at or below
        # expectation keeps the upper tail's convention of 0 (p = 1) unless the cell failed, where every tail is exact
        def slot_matrix(v):
            Pm = np.full((ncell, K_eff), np.nan, np.float32); Pm[row, rk] = v / LN10; return Pm
        comp_div_ = slot_matrix(np.where(np.isnan(lp_hi), 0.0, lp_hi)); comp_sim_ = slot_matrix(lp_in)
        # the closest match at the locus for every cell (-1: no taxon shares enough sites) and its too-similar test: the
        # comparator arm when it is one of the comparators, the closest-match arm when it is not, NaN when not tested
        lclos_ = closest[idx].astype(np.int32); lclos_p = np.full(ncell, np.nan, np.float32)
        hasy = lclos_ >= 0; viaarm = hasy & tested2 & (Y == lclos_)
        lclos_p[viaarm] = lp_lo[viaarm] / LN10
        eq = top[idx] == lclos_[:, None]; has = hasy & ~viaarm & eq.any(1)
        lclos_p[has] = comp_sim_[np.flatnonzero(has), eq.argmax(1)[has]]

        def sim_cols(q):
            """the too-similar columns of a cell: its pair's numbers from the test that fired, else its own closest-match test"""
            if q in pair_of:
                rec, other, _ = pair_of[q]
            elif tested2[q]: rec, other = dict(arm="closest", o=o2[q], mis=mis2[q], E=float(E2[q]), f=1.0, lnp=float(lp_lo[q]), rel=rel2[q]), int(Y[q])
            else: return ["", "", "", "", "", "", "", "", ""]
            qo = cpos.get(other); px = fnum(exm[qo], 3) if qo is not None else ""
            return [rec["arm"], names[other], role_.get(q, ""), fnum(100 * rec["mis"] / rec["o"], 2), fnum(100 * rec["E"] / rec["o"], 2),
                    fnum(rec["f"], 3) if rec["arm"] == "comparator" else "", fnum(rec["lnp"] / LN10, 3), "; ".join(names[c] for c in rec["rel"]), px]

        mrows_out = []
        for q in np.flatnonzero(masked_q).tolist():
            g = int(idx[q]); dv = bool(dflag[q] and mask_hi); sm = bool(sflag[q] and mask_lo)
            reason = "both" if dv and sm else ("divergent" if dv else "similar")
            s = int(starts[q]); comps = [names[int(c)] for c in ci[s:s + kc[q]]]
            row_ = [name, names[g], reason, int(nres[g]), int(nf[q]), int(kc[q]), fnum(100 - med_i[q], 2), fnum(100 - med_e[q], 2),
                    fnum(exm[q], 3), fnum(div_lp[q] / LN10, 3)] + sim_cols(q) + ["; ".join(comps), fnum(lbar10, 3)]
            mrows_out.append((g, reason, row_))

        # ---- the report's view of this locus: cell columns, best matches, site codes, sequence rows, trees, pairs
        side_ = np.zeros(ncell, np.uint8)
        side_ |= np.where(dflag, 1 if mask_hi else 16, 0).astype(np.uint8); side_ |= (amb.astype(np.uint8) << 2)
        side_ |= np.where(sflag, 2 if mask_lo else 32, 0).astype(np.uint8); side_ |= np.where(tested2, 64, 0).astype(np.uint8)
        side_ |= np.where(src_q, 128, 0).astype(np.uint8); side_ |= np.where(unc_q, 8, 0).astype(np.uint8)
        simlp_ = np.full(ncell, np.nan, np.float32); simlp_[tested2] = lp_lo[tested2] / LN10
        for q, (rec, _, _) in pair_of.items(): simlp_[q] = rec["lnp"] / LN10
        rep_q = fail_q if codes_policy == "fails" else masked_q if codes_policy == "masked" else np.zeros(ncell, bool)
        # best match: per-site divergence compared on the SAME sites. The closest taxon among those sharing the most sites
        # with the focal is the champion; a taxon closer over its own shared sites (>= BEST_MATCH_MIN_SITES) challenges
        # it on the sites all three have and wins only if it is closer there (a fragment covering a locus's conserved
        # part is not "closer"). Ties go to the closest by m. Identity is over the sites the pair shares.
        rows_bm = idx if bm_all else idx[fail_q]
        best_ = np.full(ncell, -1, np.int32); bident_ = np.full(ncell, np.nan, np.float32)
        if len(rows_bm):
            srow = Sl[tri_idx[rows_bm]].astype(np.float64); wrow = Wl[tri_idx[rows_bm]].astype(np.float64)
            okc = wrow >= BEST_MATCH_MIN_SITES; okc[np.arange(len(rows_bm)), rows_bm] = False
            with np.errstate(invalid="ignore", divide="ignore"):
                dvg = np.where(okc, (wrow - srow) / np.maximum(wrow, 1.0), np.inf)
            mmb = np.where(np.isfinite(m[rows_bm]), m[rows_bm], np.inf)
            wmax = np.where(okc, wrow, -1.0).max(1)
            dtier = np.where(okc & (wrow >= wmax[:, None]), dvg, np.inf); lo_t = dtier.min(1)
            tied = dtier == lo_t[:, None]; tkey = np.where(tied, mmb, np.inf); champ = tkey.argmin(1)
            champ = np.where(np.isfinite(tkey[np.arange(len(rows_bm)), champ]), champ, tied.argmax(1))   # none of the tied has an m: the first tied
            okb = np.isfinite(lo_t)
            Xs_ = Xs if Xs is not None else X[:, sites]
            for r_ in np.flatnonzero(okb & (dvg < lo_t[:, None]).any(1)).tolist():
                g = int(rows_bm[r_]); fx = Xs_[g]; fp = fx >= 0; c = int(champ[r_])
                T = np.array(sorted(np.flatnonzero(dvg[r_] < lo_t[r_]).tolist(), key=lambda t: (-wrow[r_, t], dvg[r_, t], mmb[r_, t])))   # broadest coverage first
                while len(T):                                                                  # every remaining challenger against the champion at once;
                    XT = Xs_[T]; common = fp[None, :] & (Xs_[c] >= 0)[None, :] & (XT >= 0)      # the first that wins takes over, and the rest face it
                    dt_ = ((XT != fx[None, :]) & common).sum(1); dc_ = ((Xs_[c] != fx)[None, :] & common).sum(1)
                    wins = (common.sum(1) >= BEST_MATCH_MIN_SITES) & ((dt_ < dc_) | ((dt_ == dc_) & (mmb[r_, T] < mmb[r_, c])))
                    if not wins.any(): break
                    k = int(wins.argmax()); c = int(T[k]); T = T[k + 1:]
                champ[r_] = c
            qb = np.array([cpos[int(g)] for g in rows_bm.tolist()], dtype=np.intp)
            best_[qb] = np.where(okb, champ, -1); bident_[qb] = np.where(okb, 100.0 * (1.0 - dvg[np.arange(len(rows_bm)), champ]), np.nan)
        codes_out = None
        if rep_q.any():
            # k = the comparators present that differ from the focal residue; where the cell has a too-similar partner,
            # the partner is left out of k and the sites where it agrees with the focal while at least two comparators
            # differ carry CODE_SIMILAR + k (the blue of a copy), so one row holds both tests' evidence
            Xs_ = Xs if Xs is not None else X[:, sites]
            qs = np.flatnonzero(rep_q); blk_c = np.empty((len(qs), len(sites)), np.uint8)
            step = max(1, (16 << 20) // max(K_eff * len(sites), 1))                    # cells per slice: the slice's block stays under 16 MB
            for c0 in range(0, len(qs), step):
                qq = qs[c0:c0 + step]; G = idx[qq]; FX = Xs_[G]; pres = FX >= 0
                Bq = np.where(copy_q[qq] | src_q[qq], partner_[qq], -1)
                C = top[G]; cv = (C >= 0) & (C != Bq[:, None]); XC = Xs_[np.maximum(C, 0)]
                k_ = ((XC != FX[:, None, :]) & (XC >= 0) & cv[:, :, None]).sum(1)
                rows_c = np.where(pres, np.minimum(k_, 127), CODE_ABSENT).astype(np.uint8)
                hb = Bq >= 0
                if hb.any():
                    XB = Xs_[np.maximum(Bq, 0)]; shared = hb[:, None] & pres & (XB >= 0) & (FX == XB) & (k_ >= 2)
                    rows_c[shared] = (CODE_SIMILAR + np.minimum(k_[shared], 125)).astype(np.uint8)
                blk_c[c0:c0 + step] = rows_c
            codes_out = (idx[qs].astype(np.int32), blk_c)
        # face-to-face trees for the masked cells: the focal, its comparators and its locus set, plus an outgroup to root
        # both trees on: the taxon at the locus sharing >= fit-overlap sites with every tree taxon that is the most distant
        # from them genome-wide (largest mean m)
        trees_out = None
        if masked_q.any():
            Xs_ = Xs if Xs is not None else X[:, sites]
            tcells, texp, tobs, ttips, togs, tclos = [], [], [], [], [], []
            for q in np.flatnonzero(masked_q).tolist():
                g = int(idx[q]); tips = [g] + [int(c) for c in top[g] if c >= 0]
                if partner_[q] >= 0 and int(partner_[q]) not in tips: tips.append(int(partner_[q]))
                if len(tips) < 2: continue
                og = -1; cand = np.ones(n, bool); cand[tips] = False; cidx = np.flatnonzero(cand)
                if len(cidx):
                    tarr0 = np.array(tips); Wsub = Wl[tri_idx[np.ix_(cidx, tarr0)]]; msub = m[np.ix_(cidx, tarr0)]
                    okog = (Wsub >= a.fit_overlap).all(1) & np.isfinite(msub).all(1)
                    if okog.any(): og = int(cidx[okog][np.argmax(msub[okog].mean(1))]); tips.append(og)
                tarr = np.array(tips); nt = len(tips); rows_t = Xs_[tarr]; pres_t = rows_t >= 0
                both = pres_t[:, None, :] & pres_t[None, :, :]; nb = both.sum(2)
                with np.errstate(invalid="ignore", divide="ignore"):
                    Do = np.where(nb > 0, 100.0 * ((rows_t[:, None, :] != rows_t[None, :, :]) & both).sum(2) / np.maximum(nb, 1), np.nan)
                    mv = m[np.ix_(tarr, tarr)]; De = np.where(np.isfinite(mv), expect(np.where(np.isfinite(mv), mv, 1.0)), np.nan)
                np.fill_diagonal(De, 0.0); np.fill_diagonal(Do, 0.0)
                if np.isnan(De).any(): De = np.where(np.isnan(De), np.nanmax(De) if np.isfinite(np.nanmax(De)) else 50.0, De)
                if np.isnan(Do).any(): Do = np.where(np.isnan(Do), De, Do)
                # both matrices on the method's own scale: the F81 correction with the locus plateau, so the expected
                # tree is the genome-wide tree scaled by the locus rate (r c(m)) and the observed tree is c(d) over shared
                # sites, capped at c(0.98 P) like every other corrected distance in the tool
                with np.errstate(invalid="ignore", divide="ignore"):
                    De = corr(De); Do = corr(Do)
                np.fill_diagonal(De, 0.0); np.fill_diagonal(Do, 0.0)
                t = cell_trees(De, Do, [str(x) for x in tips], tips.index(og) if og >= 0 else None)
                if t:
                    # the tip closest to the focal in each tree, by that tree's own distances (the outgroup excluded)
                    cand = np.ones(nt, bool); cand[0] = False
                    if og >= 0: cand[tips.index(og)] = False
                    ce = int(np.argmin(np.where(cand, De[0], np.inf))) if cand.any() else -1; co = int(np.argmin(np.where(cand, Do[0], np.inf))) if cand.any() else -1
                    tcells.append(g); texp.append(t[0]); tobs.append(t[1]); ttips.append(tips); togs.append(tips.index(og) if og >= 0 else -1); tclos.append([ce, co])
            if tcells: trees_out = dict(cells=tcells, expected=texp, observed=tobs, tips=ttips, outgroup=togs, closest=tclos)
        seq_out = None
        if seq_policy != "none":
            if seq_policy == "all": rows_s = np.flatnonzero(nres > 0)
            else:
                want = fail_q if seq_policy == "fails" else masked_q
                rs = set()
                for q in np.flatnonzero(want).tolist():
                    g = int(idx[q]); rs.add(g); rs.update(int(c) for c in top[g] if c >= 0)
                    if partner_[q] >= 0: rs.add(int(partner_[q]))
                rows_s = np.array(sorted(rs), dtype=np.int64)
            if len(rows_s):
                seq_out = (rows_s.astype(np.int32), np.stack([views[g][sites] for g in rows_s.tolist()]))
        sim_json = ([dict(test=rec["arm"], taxon=g, partner=y, overlap=int(rec["o"]), mismatches=int(rec["mis"]), expected=round(rec["E"], 3),
                          row_factor=round(rec["f"], 4), relatives=rec["rel"], log10p=round(rec["lnp"] / LN10, 3), fail=True,
                          copy=(-1 if rec.get("uncertain") else int(g if role_.get(cpos.get(g)) == "copy" else y)), uncertain=bool(rec.get("uncertain", False)),
                          excess_m=[(round(float(exm[cpos[x]]), 4) if x in cpos and np.isfinite(exm[cpos[x]]) else None) for x in (g, y)])
                     for (g, y), rec in pairs.items()]
                    + [dict(test="closest", taxon=int(idx[q]), partner=int(Y[q]), overlap=int(o2[q]), mismatches=int(mis2[q]), expected=round(float(E2[q]), 3),
                            row_factor=1.0, relatives=[int(c) for c in rel2[q]], log10p=round(float(lp_lo[q]) / LN10, 3), fail=False, copy=-1,
                            excess_m=[(round(float(exm[q]), 4) if np.isfinite(exm[q]) else None), None])
                       for q in np.flatnonzero(tested2 & ~fail2).tolist()])
        repblk = dict(idx=idx, n_comp=kc.astype(np.uint8), n_fail=nf.astype(np.uint8), n_similar=np.minimum(npair, 255).astype(np.uint8), side=side_,
                      med_identity=med_i, med_expected=med_e, excess_m=exm, closest=ci[starts].astype(np.int32), closest_overlap=o[starts].astype(np.uint32),
                      closest_identity=ident[starts], closest_expected=100.0 - dhat[starts], best_match=best_, best_match_identity=bident_,
                      partner=partner_, div_log10p=div_lp / LN10, sim_log10p=simlp_,
                      comparators=top[idx].astype(np.int32), fail_q=fail_q, codes=codes_out, seq=seq_out, trees=trees_out, similar=sim_json,
                      locus_closest=lclos_, locus_closest_log10p=lclos_p, comp_div_log10p=comp_div_, comp_sim_log10p=comp_sim_)

        # ---- every test of the locus, for --all-tests: one row per comparison, the focal sequence first
        tests_block = ""
        if every:
            clean = lambda x: x.replace("\t", " ").replace("\n", " ").replace("\r", " ")
            lname = clean(name).replace("%", "%%"); tnames = [clean(x) for x in names]
            r3 = lambda v: (np.round(np.asarray(v, float), 3) + 0.0).tolist()
            tests_block = ("divergent\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t\t%.3f\t%d\n") * len(gi) % tuple(chain.from_iterable(zip(
                [tnames[x] for x in gi.tolist()], [tnames[x] for x in ci.tolist()], o.astype(np.int64).tolist(), mis.astype(np.int64).tolist(),
                r3(E), r3(lp_hi / LN10), fail_hi.astype(np.int64).tolist())))
            tests_block += ("similar_comparator\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t%.4f\t%.3f\t%d\n") * len(gi) % tuple(chain.from_iterable(zip(
                [tnames[x] for x in gi.tolist()], [tnames[x] for x in ci.tolist()], o.astype(np.int64).tolist(), mis.astype(np.int64).tolist(),
                r3(E_lo), ff.tolist(), r3(lp_in / LN10), fail_in.astype(np.int64).tolist())))
            tq_ = np.flatnonzero(tested2)
            if len(tq_):
                tests_block += ("similar_closest\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t\t%.3f\t%d\n") * len(tq_) % tuple(chain.from_iterable(zip(
                    [tnames[int(idx[q])] for q in tq_.tolist()], [tnames[int(Y[q])] for q in tq_.tolist()], o2[tq_].astype(np.int64).tolist(),
                    mis2[tq_].astype(np.int64).tolist(), r3(E2[tq_]), r3(lp_lo[tq_] / LN10), fail2[tq_].astype(np.int64).tolist())))

        log = (f"{name}: {ncell} cells, {len(gi)} comparator slots, {int(tested2.sum())} closest-match tests; too divergent {int(flag.sum())}"
               f" (+{int(amb.sum())} with fewer comparators), too-similar pairs {len(pairs)} ({sum(1 for r in pairs.values() if r['arm'] == 'closest')} by the closest match)")
        return dict(skipped=False, mrows=mrows_out, tests=tests_block, log=log, rep=repblk, phi=phi, rate_c=rate_c)

    mrows = []; tally = collections.Counter()
    def consume(li, R):
        name, sites = loci[li]
        if R["skipped"]:
            LOG.detail(R["log"]); rep["loc_status"][li] = R["status"]; PROG.tick(); return
        for g, reason, row_ in R["mrows"]:
            views[g][sites] = 45                                    # 45 = '-': the whole stretch of the locus
            mrows.append(row_); tally[reason] += 1
        if tests_fh is not None: tests_fh.write(R["tests"])
        LOG.detail(R["log"])
        r = R["rep"]; idx_ = r["idx"]; rep["phi"][li] = R["phi"]; rep["test_rate"][li] = R["rate_c"]
        for k in ("n_comp", "n_fail", "n_similar", "side"): rep[k][idx_, li] = r[k]
        if dense:
            for k, _, _ in NUMCOLS: rep[k][idx_, li] = r[k]
            rep["comparators"][idx_, li] = r["comparators"]
            for k in ("comp_div_log10p", "comp_sim_log10p"): rep[k][idx_, li] = r[k]
        else:
            sel_ = r["fail_q"]; sparse_loc.append(np.full(int(sel_.sum()), li, np.int32)); sparse_tax.append(idx_[sel_].astype(np.int32))
            for k, _, _ in NUMCOLS: rep[k].append(np.asarray(r[k])[sel_])
            rep["comparators"].append(r["comparators"][sel_])
            for k in ("comp_div_log10p", "comp_sim_log10p"): rep[k].append(r[k][sel_])
        if r["codes"] is not None: rep["codes"][li] = r["codes"]
        if r["seq"] is not None: rep["seqs"][li] = r["seq"]
        if r["trees"] is not None: rep["trees"][li] = r["trees"]
        if r["similar"]: rep["similar"][li] = r["similar"]
        PROG.tick(extra=f"{tally['divergent'] + tally['both']} too divergent, {tally['similar'] + tally['both']} copied")

    # loci in batches of about TEST_BATCH_COMPARISONS comparisons: the exact tails of a batch are computed in one call
    # per kind (upper, lower), the loci answered in turn, their results consumed in locus order. The stage is Python
    # work per locus around those calls, so it runs in worker PROCESSES (forked: they inherit every array without a
    # copy; threads only queue on the interpreter lock here), one batch per task, the batches sized to keep them busy
    if pool: pool.shutdown(wait=True); pool = None                     # no live threads when the workers fork
    ncmp_l = [int(t[1][t[2] > 0].sum()) for t in sel]
    workers = min(a.threads, nloc) if a.threads > 1 and nloc > 1 and "fork" in multiprocessing.get_all_start_methods() else 1
    target = max(TEST_BATCH_MIN, min(TEST_BATCH_COMPARISONS, sum(ncmp_l) // (4 * workers))) if workers > 1 else TEST_BATCH_COMPARISONS
    batches = []; cur = []; acc = 0
    for li in range(nloc):
        cur.append(li); acc += ncmp_l[li]
        if acc >= target: batches.append(cur); cur = []; acc = 0
    if cur: batches.append(cur)
    workers = min(workers, len(batches), max(1, _mem_room() // (TEST_BATCH_COMPARISONS * 64 * 80)))   # a batch's tails: about 256 MB of temporaries
    if a.threads > 1: LOG(f"test stage: {workers} worker process{'es' if workers > 1 else ''} ({len(batches)} batch{'es' if len(batches) > 1 else ''} of loci, about {target:,} comparisons each)")
    def run_batch(lis):
        gens, reqs, done = {}, {}, {}
        for li in lis:
            g = process_locus(li)
            try: reqs[li] = next(g); gens[li] = g
            except StopIteration as e: done[li] = e.value
        while reqs:
            ans = {li: [None] * len(r) for li, r in reqs.items()}
            for kind, upper in (("up", True), ("dn", False)):
                slots = [(li, k) for li, r in reqs.items() for k, x in enumerate(r) if x[0] == kind]
                if not slots: continue
                parts = [reqs[li][k] for li, k in slots]
                res = bb_tails(np.concatenate([x[1] for x in parts]), np.concatenate([x[2] for x in parts]), np.concatenate([x[3] for x in parts]),
                               np.concatenate([np.full(len(x[1]), phi_l[li]) for (li, _), x in zip(slots, parts)]), upper)
                pos = 0
                for (li, k), x in zip(slots, parts): ans[li][k] = res[pos:pos + len(x[1])]; pos += len(x[1])
            nxt = {}
            for li in list(reqs):
                try: nxt[li] = gens[li].send(ans[li])
                except StopIteration as e: done[li] = e.value
            reqs = nxt
        return [done[li] for li in lis]
    PROG.start("test", nloc, "0 too divergent, 0 copied")
    if workers > 1:
        global _TEST_BATCH; _TEST_BATCH = run_batch                       # the forked workers find it in their copy of this module
        with multiprocessing.get_context("fork").Pool(workers) as mp:
            for lis, Rs in zip(batches, mp.imap(_test_batch, batches)):   # results in locus order, consumed as they arrive
                for li, R in zip(lis, Rs): consume(li, R)
        _TEST_BATCH = None
    else:
        for lis in batches:
            for li, R in zip(lis, run_batch(lis)): consume(li, R)
    PROG.finish()
    if tests_fh is not None: tests_fh.close()
    # the report's remaining pieces: sparse rows stacked, the closest taxa by m, the test's constants
    if not dense:
        rep["sparse_locus"] = np.concatenate(sparse_loc) if sparse_loc else np.zeros(0, np.int32)
        rep["sparse_taxon"] = np.concatenate(sparse_tax) if sparse_tax else np.zeros(0, np.int32)
        for k, dt, _ in NUMCOLS: rep[k] = np.concatenate(rep[k]).astype(dt) if rep[k] else np.zeros(0, dt)
        rep["comparators"] = np.concatenate(rep["comparators"]) if rep["comparators"] else np.zeros((0, K_eff), np.int32)
        for k in ("comp_div_log10p", "comp_sim_log10p"): rep[k] = np.concatenate(rep[k]) if rep[k] else np.zeros((0, K_eff), np.float32)
    cl = order[:, :min(5, n)].astype(np.int32); cl[~np.isfinite(mj[np.arange(n)[:, None], cl])] = -1   # the closest non-identical relatives
    rep.update(taxa_closest=cl, rate=rate, plateaus=plateaus, kq=kq, logbar=logbar, phi_bar=phi_bar,
               tests=dict(divergent=D_tests, similar_comparator=D_tests, similar_closest=S_tests, total=N_tests, locus_sets=nfork))
    with open(f"{a.out}/masked.csv", "w", newline="") as fh:
        mw = csv.writer(fh); mw.writerow(MASKED_COLUMNS); mw.writerows(mrows)
    skipped = [(loci[i][0], st) for i, st in enumerate(rep["loc_status"]) if st != "tested"]
    if skipped:
        LOG(f"{len(skipped)} loci not tested: " + "; ".join(f"{nm} ({st})" for nm, st in skipped[:20]) + (f"; and {len(skipped) - 20} more" if len(skipped) > 20 else ""))
    present_c = rep["focal_sites"] > 0; tested_c = rep["n_comp"] > 0; untested = present_c & ~tested_c
    thin_t = np.zeros(n, bool); thin_t[rep["thin"]] = True; skipped_l = np.array([st != "tested" for st in rep["loc_status"]])
    u_thin = int((untested & thin_t[:, None]).sum()); u_loc = int((untested & ~thin_t[:, None] & skipped_l[None, :]).sum()); u_comp = int(untested.sum()) - u_thin - u_loc
    rep["untested"] = dict(sequences=int(untested.sum()), of_excluded_taxa=u_thin, at_untested_loci=u_loc, without_a_comparator=u_comp)
    LOG(f"not tested: {len(rep['thin'])} of {n} taxa (excluded: fewer than {MIN_LOCI_PER_TAXON} loci with at least {a.fit_overlap} residues"
        f"{': ' + ', '.join(names[i] for i in rep['thin'][:10]) + (' ...' if len(rep['thin']) > 10 else '') if len(rep['thin']) else ''}); "
        f"{int(skipped_l.sum())} of {nloc} loci ({', '.join(f'{loci[i][0]} ({st})' for i, st in enumerate(rep['loc_status']) if st != 'tested')[:300] or 'none'}); "
        f"{int(untested.sum()):,} of {int(present_c.sum()):,} present sequences ({u_thin:,} of excluded taxa, {u_loc:,} at untested loci, {u_comp:,} with no comparator: no resolved relative sharing enough sites)")
    npair = sum(sum(1 for e in v if e['fail']) for v in rep['similar'].values()); ncl = sum(sum(1 for e in v if e['fail'] and e['test'] == 'closest') for v in rep['similar'].values())
    LOG(f"too-similar pairs: {npair:,} ({npair - ncl:,} against a comparator, {ncl:,} against the closest match; {S_tests:,} closest-match tests run)")
    nunc = int(((rep["side"] & 8) > 0).sum())
    LOG(f"masked {len(mrows)} sequences: {tally['divergent']} too divergent, {tally['similar']} too similar (copied from another taxon"
        f"{f'; {nunc} of a pair whose direction is uncertain, both members masked' if nunc else ''}), {tally['both']} both")
    return len(mrows)


if __name__ == "__main__":
    main()
