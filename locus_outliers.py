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
                       mask_test and mask_p - the sequence's one significance: the most significant
                       BAR-HELD p over the tests that masked it (the close-relative decisive p, the fired
                       closest-match p, the broad-sample stage-1 variance p) and the test it came from;
                       a rough level of belief, not a combined p. Every p in the text outputs is plain,
                       at six decimals (1e-7 and below print 0.000000; report.npz keeps exact values) -
                       residues_masked (its residues at the locus; the whole stretch is blanked), then
                       one column family per test, each with its effect and its p where the test ran.
                       CLOSE-RELATIVE test: close_relative_fails of close_relative_comparators, the
                       observed and expected divergence (%; medians over the comparators),
                       close_relative_excess_m (the extra divergence beyond expectation, multiple-hit
                       corrected and divided by the locus rate: % on the scale of m) and
                       close_relative_p, the p of the comparison that decides the majority.
                       CLOSEST-MATCH test: the fired pair when it fired, else the strongest comparison
                       it made; closest_match_arm names the expectation used - "relative" (the partner
                       is one of the focal's relatives: the pair machinery, with closest_match_row_factor)
                       or "outside" (the empirical expectation from the relatives' own distances to the
                       partner, with closest_match_relatives). closest_match_role (copy, source or
                       uncertain - the direction is called when the larger of the two members'
                       divergence excesses is more than --donor-threshold times the smaller, else both
                       are masked; empty when the test did not fire), observed and expected divergence
                       (%; before the row factor), closest_match_p and
                       closest_match_partner_excess_m (the partner's divergence excess, next to this
                       sequence's close_relative_excess_m: the margin of the copy call). BROAD-SAMPLE
                       test: broad_sample_partner (the most extreme aligned member of the sample, by
                       standardised residual), broad_sample_side (close or divergent - the test calls
                       one direction at most, the stronger side), observed and expected (%) against the
                       partner and broad_sample_p (on a divergent row with no broad mask these are
                       the donor-side label). Then its comparators, closest first, and p_bar, the
                       bar every p was held to.
  tests.tsv.gz         with --all-tests only: one row per test performed, the focal sequence first.
                       test (close_relative, closest_match_relative, closest_match_outside; one row
                       per comparison), locus, taxon, comparator,
                       overlap (shared sites), mismatches, expected (expected mismatches; before the row
                       factor), row_factor (closest_match_relative rows), p (six decimals), fail (0 or 1).
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
  --no-numba (the pair stage in numpy even when numba is installed: the compiled kernel gives the
  same numbers, 3-8x faster at scale, so this switch is for checking that)
  --candidates 200 (the closest genomes each taxon keeps as its candidate comparators; a taxon whose
  list runs out at a locus gets its full row computed)  --block-taxa N (the taxon block of the pair sweep;
  default from the taxa count, mainly for testing: any block size gives the same results)
  (constants in the code: HALF = 0.5 mismatch for the resolution floor and the rate
  pseudocount; dispersion from pairs with >= 10 expected mismatches, shrinkage weight 30)
"""
import argparse, base64, collections, csv, gzip, hashlib, io, json, math, os, platform, re, shutil, subprocess, sys, threading, time, types, zipfile
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
BS_PER_COMP = 10        # broad sample per cell: this many per requested comparator (50 at the default 5)
BS_POOL = 100           # each taxon's pool of candidate broad sample: this many per requested comparator (500)
BS_MIN = 5              # a cell with fewer valid broad sample than this is not tested
BS_VAR_SCREEN = 4.0     # a comparator whose own row variance exceeds this times the locus median is dropped
BS_SCREEN_LOGP = -6.907755278982137   # log 1e-3: a comparator whose own side p is below this is dropped from that side
BS_SEED = 0x5eed1517    # the constant behind the pools and the per-cell draws (our own hashes: no library generator)
PASS1_LOCI_PER_WORKER = 64             # pass one with small loci: a worker process per this many loci
PASS1_RETURN_MAX = 512 << 20           # ... when the pair counts the workers send back are at most this many bytes
PAIR_CHUNK_MIN = 32 << 20              # the pair stage's chunks of pairs hold between these many bytes of temporaries
PAIR_CHUNK_MAX = 256 << 20
DISPERSION_PAIRS_PER_WORKER = 200_000  # the dispersion pre-pass takes a worker process per this many comparisons (Ivan_2026, 4.2 M
                                       # comparisons in a 10 GB process: 1 worker 0.8 s, 8 workers 0.7 s, 42 workers 1.8 s: forks cost)
COMPARATOR_SECONDS_PER_WORKER = 0.3    # the comparator stage takes a worker process per this many seconds of estimated walk (measured 31 Oct 2026: about
                                       # 0.5 ms + 4 us per taxon per locus from the kept counts, 20 x taxa x sites / 1.5e8 s per locus from the residues)
TEST_BATCH_COMPARISONS = 50_000   # the test computes the exact tails of this many comparisons at a time (loci batched), so small
                                  # loci do not pay numpy's per-call overhead hundreds of times each; about 250 MB of temporaries per batch
TEST_BATCH_MIN = 5_000            # with worker processes the batches shrink to keep every worker busy, but never below this
PAIRS_PER_THREAD = 4000           # the pair stage uses a thread per this many pairs: beyond that the threads contend (Shen_2018,
                                  # 59K pairs x 2,408 loci: 16 threads 1.3 s, 32 threads 3.8 s, 96 threads 5.2 s; Ivan_2026, 246K
                                  # pairs: 96 threads 0.9 s)
NUMBA_PAIR_LOCI_MIN = 2_000_000   # the pair stage runs its compiled kernel (numba, when installed) from this many loci x pairs:
                                  # below it the half second numba takes to import is more than numpy needs for the stage
NUMBA_PAIR_SLAB = 1 << 22         # pairs per kernel call (the progress bar moves once per slab)
PAIR_BLOCK = 256                  # pairs a kernel block gathers at a time (its scratch is 32 bytes per locus x pair)
USE_NUMBA = True                  # --no-numba clears it: the numpy pair stage whatever the size
SWEEP_BLOCK_MIN = 1024            # up to this many taxa the pair sweep is one block pair whose per-locus work runs over the loci in threads or
                                  # processes, as pass one always did; beyond it the taxa are cut into blocks, a block pair per worker at a time
SWEEP_BLOCKS = 2.0                # ... about this many times the square root of the workers blocks, so each worker gets a few block pairs
SWEEP_BLOCK_FLOOR = 32            # ... never smaller than this (a block pair is then too little work for its overheads)
SWEEP_FLOPS_PER_WORKER = 2e10      # ... and the workers the blocks are cut for are at most one per this many floating-point operations of products in the
                                  # sweep (about a third of a second each): a small sweep over many taxa is not cut into more block pairs than its work is worth
CLOSEST_THREADS_MIN = 200_000     # the one block pair's closest reduce runs in threads from this many positions (taxa^2) per locus; below it the
                                  # per-locus arrays are too small for threads to pay (they queue on the interpreter lock) and the loci go to worker processes
RESTART_ARRAYS = dict(rates=("ident", "rate", "ok_r"), pairs=("ident", "rate", "ok_r", "cand", "key", "m", "mfl", "nfin", "cd", "cv", "cc", "counters"))
STOP_AFTER = None                 # test hook (LOCUSMASKER_STOP_AFTER=rates|pairs): the run stops right after writing that restart point
SWEEP_RESULTS_DIR = "/dev/shm"    # the sweep's worker processes leave each block pair's results in a file here (a memory-backed filesystem) when
                                  # it exists and has the room, else they pickle them through the pool's one pipe, which is far slower
CANDIDATES = 200                  # each taxon keeps this many closest candidates by the ranking key (--candidates); a taxon whose list runs
                                  # out in the comparator walk gets its full row computed afterwards
FULL_M_MAX_TAXA = 20000           # up to this many taxa the sweep also keeps the full m matrix (a few GB) for the consumers that read whole
                                  # rows (best-match ties, trees); above it those rows are computed on demand and the trees are rooted at
                                  # their midpoint (no outgroup search). A rule on the taxa count, never on memory: outputs never depend
                                  # on the machine
BAND_SLOTS = 4                    # per (taxon, locus) the closest-match band keeps one entry per distinct mismatch rate within 1e-9 of the
                                  # minimum, up to this many (one is the rule; more needs two overlaps with 2 o1 o2 >= 1e9)
LSET_ELEMENTS = 1 << 20           # a locus's closest matches are found in taxa x taxa slices of about this many elements (rows = this / taxa)
__version__ = "0.13.1"
REPORT_FORMAT = 9                 # bundle format version; the viewer refuses a newer one (3: tree branch lengths are F81-corrected; 4: the too-similar test over every pair; 7: every test of every cell, the closest match at the locus, the uncertain-direction pairs; 8: the broad-sample test, cells/side2 and cells/bs_*; 9: the broad-sample call - one direction at most, both sides' statistics and the standardised-residual representative for every gate-passer: cells/bs_hi_log10p, cells/bs_call; from 0.13.1 also cells/bs_rep, the stage-1 representative of every tested cell, and its sequence among the masked cells' carried rows)
DENSE_CELLS_MAX = 8_000_000       # per-cell numeric report columns are dense (taxa x loci) up to this many cells, sparse beyond
BEST_MATCH_ALL_MAX_TAXA = 2000    # best matches for every tested cell up to this many taxa; beyond, only for cells with a fail
BEST_MATCH_MIN_SITES = 30         # a best match needs this many sites shared with the focal (as plot_panels.best_matches)
REPORT_SEQ_ALL_MAX = 32 << 20     # --report-sequences auto carries the whole alignment up to this many bytes, else the cells' rows
REPORT_ONE_FILE_MAX = 256 << 20   # report.html carries the bundle (report.npz) inline up to this many bytes of bundle: a browser holds about four times
                                  # the bundle (the base64 text, its bytes, the decoded arrays), 1 GB here, which a laptop opens in seconds; above it the
                                  # run writes report.npz and the tables (taxa.tsv, loci.tsv, summary.json) and says so
CODE_SIMILAR, CODE_TICK, CODE_ABSENT = 128, 254, 255
PASS1_FLOPS_PER_THREAD = 6e7    # pass one adds a thread per this many floating-point operations per locus: below it the
                                # interpreter lock around the small matrix products costs more than the thread gains
                                # (measured: 90 taxa x 590 sites, 4e7 flops per locus with the per-letter products, is
                                # fastest at 1-2 threads and slower at 4 or more)


# ------------------------------------------------------------------ input
def read_fasta(txt):
    """names and sequences, one per header line, in order: a header with nothing under it gives an empty sequence
    (read_input then rejects it by name); text before the first header is an error, not a nameless sequence. txt may
    be str or bytes (the sequences come back in the same type, the names as str). Works record by record (split at the
    line-start '>') with whole-string operations, not line by line: a wrapped 1.4 GB file has 18 million lines"""
    isb = isinstance(txt, (bytes, bytearray)); gt, nl, cr, sp, tab, empty = (b">", b"\n", b"\r", b" ", b"\t", b"") if isb else (">", "\n", "\r", " ", "\t", "")
    body = txt.lstrip()
    if not body: return [], []
    if not body.startswith(gt):
        first = body.splitlines()[0].strip()[:40]
        sys.exit(f"FASTA: text before the first header line ('>'): {first.decode('utf-8', 'replace') if isb else first!r} (a sequence whose header line is missing?)")
    names, seqs = [], []
    for rec in body[1:].split(nl + gt):
        header, _, seq = rec.partition(nl)
        name = header.split()
        if not name: sys.exit("FASTA: a header line ('>') has no name")
        nm = name[0]
        if isb:
            try: nm = nm.decode("utf-8")
            except UnicodeDecodeError: sys.exit(f"FASTA: the name {nm[:40]!r} is not UTF-8 text")
        seq = seq.replace(nl, empty).replace(cr, empty)
        if seq.count(sp) or seq.count(tab): seq = empty.join(seq.split())
        names.append(nm); seqs.append(seq)
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
            what = f"non-ASCII character {ch!r}" if ch != "\ufffd" else f"byte 0x{int(b[hi[0]]):02X} (not UTF-8 text)"
            sys.exit(f"{what} in the sequence of taxon {t} at site {int(hi[0]) + 1:,} of {os.path.basename(path)}: sequences may hold only ASCII letters, '-' and '?'")
    lens = collections.Counter(len(r) for r in rows); L = lens.most_common(1)[0][0]
    bad = [(t, len(r)) for t, r in zip(names, rows) if len(r) != L]
    if bad:
        shown = ", ".join(f"{t} has {l:,}" for t, l in bad[:5]) + (f", and {len(bad) - 5} more" if len(bad) > 5 else "")
        sys.exit(f"the sequences in {os.path.basename(path)} are not all the same length: {shown}, while {lens[L]} of the {len(rows)} sequences have {L:,} sites. "
                 "The input must be an alignment: every sequence padded with gaps to the same length")


def read_input(path, partitions):
    """returns names, rows (the alignment as one (taxa, sites) uint8 array: a row per taxon), loci, format, datatype
    declared by a NEXUS file (or None), and the NEXUS record the writer needs (a Bio.Nexus object, or NexusInfo from the
    fast reader)"""
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
        if head.startswith(">"):                                        # FASTA is read as bytes: no decoding and re-encoding of gigabytes
            with open(path, "rb") as fh: raw = fh.read()
            names, seqs = read_fasta(raw); fmt = "fasta"; del raw
            empty = [t for t, sq in zip(names, seqs) if not sq]
            if empty: sys.exit(f"FASTA: {len(empty)} record{'s' if len(empty) > 1 else ''} with no sequence under the header in {path}: {', '.join(empty[:5])}{' ...' if len(empty) > 5 else ''}")
            rows = [bytearray(sq) for sq in seqs]; del seqs
        else:
            if not head[:1].isdigit(): sys.exit(f"cannot tell the format of {path}: a NEXUS file starts with #NEXUS, a FASTA file with '>', a PHYLIP file with its counts")
            try:
                with open(path, encoding="utf-8") as fh: txt = fh.read()
            except UnicodeDecodeError: sys.exit(f"the alignment {path} is not a UTF-8 text file (compressed? binary?)")
            names, seqs = read_phylip(txt); fmt = "phylip"
            rows = [bytearray(sq.encode()) for sq in seqs]; del seqs, txt
        _check_rows(path, names, rows)
        loci = parse_raxml_partitions(partitions, len(rows[0]))
    _check_rows(path, names, rows)
    dup = sorted(t for t, c in collections.Counter(names).items() if c > 1)
    if dup: sys.exit(f"duplicate taxon names: {dup[:5]}{' ...' if len(dup) > 5 else ''}")
    if not loci: sys.exit("no loci defined")
    R = np.empty((len(names), len(rows[0])), np.uint8)                # one contiguous array: forking a process that holds thousands of small
    for i in range(len(rows)): R[i] = np.frombuffer(rows[i], dtype=np.uint8); rows[i] = None   # allocations costs 20x more than one mapping
    return names, R, loci, fmt, declared, nexus


# ------------------------------------------------------------------ output
def write_masked(path, out_dir, fmt, names, rows, nexus):
    stem, ext = os.path.splitext(os.path.basename(path)); out = os.path.join(out_dir, f"{stem}.masked{ext or '.' + fmt}")
    if fmt == "nexus" and isinstance(nexus, NexusInfo):
        with open(out, "wb") as fh:
            fh.write(b"#NEXUS\nbegin data;\n" + f"dimensions ntax={len(names)} nchar={len(rows[0])};\n".encode())
            fh.write(("format datatype=" + nexus.datatype + (" missing=" + nexus.missing if nexus.missing else "")
                      + (" gap=" + nexus.gap if nexus.gap else "") + ";\nmatrix\n").encode())
            pad = max(len(_safename(t)) for t in names)
            for t, r in zip(names, rows): fh.write(_safename(t).ljust(pad + 1).encode() + bytes(r) + b"\n")
            fh.write(b";\nend;\n")
            if nexus.charsets or nexus.taxsets:
                lines = ["\nbegin sets"] + [f"charset {_safename(k)} = {_compact4nexus(v)}" for k, v in nexus.charsets.items() if v] + \
                        [f"taxset {_safename(k)} = {' '.join(_safename(t) for t in v)}" for k, v in nexus.taxsets.items() if v] + ["end;\n"]
                fh.write(";\n".join(lines).encode())
    elif fmt == "nexus":
        try: from Bio.Seq import Seq
        except ImportError: sys.exit("writing the masked alignment back into this NEXUS file needs Biopython, which is not installed")
        for t, r in zip(names, rows): nexus.matrix[t] = Seq(bytes(r).decode())
        nexus.write_nexus_data(filename=out, interleave=False, append_sets=True)
    elif fmt == "fasta":
        with open(out, "wb") as fh:
            for t, r in zip(names, rows): fh.write(f">{t}\n".encode() + bytes(r) + b"\n")
    else:
        with open(out, "wb") as fh:
            fh.write(f"{len(names)} {len(rows[0])}\n".encode())
            for t, r in zip(names, rows): fh.write(f"{t} ".encode() + bytes(r) + b"\n")
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


# ---------------------------------------------------------------- the broad-sample test's special functions
# (numpy is the only hard dependency, so the few pieces of scipy.special the scale null needs are implemented here;
# the unit tests check each against scipy where scipy is installed)

_LANCZOS = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
            -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)


def gammaln_np(x):
    """log Gamma(x) for x > 0 (Lanczos, g = 7): about 1e-13 relative, enough for log-binomials and tail prefactors"""
    x = np.asarray(x, np.float64); xm = x - 1.0
    s = np.full(x.shape, _LANCZOS[0])
    for k in range(1, 9): s = s + _LANCZOS[k] / (xm + k)
    t = xm + 7.5
    return 0.5 * math.log(2.0 * math.pi) + (xm + 0.5) * np.log(t) - t + np.log(s)


def _poly_shift(x, m=8):
    """x shifted up by m with the recurrence terms 1/(x+j) and 1/(x+j)^2, for digamma and trigamma"""
    x = np.asarray(x, np.float64); j = np.arange(m, dtype=np.float64)
    inv = 1.0 / (x[..., None] + j)
    return x + m, inv.sum(-1), (inv * inv).sum(-1)


def digamma_np(x):
    """psi(x) for x > 0: recurrence to x + 8, then the asymptotic series (about 1e-14 absolute there)"""
    y, s1, _ = _poly_shift(x); y2 = 1.0 / (y * y)
    ser = np.log(y) - 0.5 / y - y2 * (1.0 / 12 - y2 * (1.0 / 120 - y2 * (1.0 / 252 - y2 / 240)))
    return ser - s1


def trigamma_np(x):
    """psi'(x) for x > 0, by the same shift"""
    y, _, s2 = _poly_shift(x); inv = 1.0 / y; y2 = inv * inv
    ser = inv + 0.5 * y2 + inv * y2 * (1.0 / 6 - y2 * (1.0 / 30 - y2 * (1.0 / 42 - y2 / 30)))
    return ser + s2


def chi2_logsf_np(x, nu):
    """log P(chi-square_nu >= x), elementwise, x >= 0: the regularised upper incomplete gamma Q(nu/2, x/2), by the
    Kummer series for the lower tail where x is small (log1p keeps the far-from-tail side exact) and the Lentz
    continued fraction for Q where x is large; the same split as scipy's gammainc"""
    x = np.asarray(x, np.float64); a = np.asarray(nu, np.float64) / 2.0; s = np.broadcast_arrays(x / 2.0, a)[0].copy(); a = np.broadcast_arrays(x / 2.0, a)[1]
    out = np.zeros(s.shape, np.float64)
    low = s < a + 1.0
    if low.any():                                                   # P = s^a e^-s / Gamma(a+1) * sum_n s^n / prod(a+1..a+n)
        sl, al = s[low], a[low]; term = np.ones(sl.shape); total = np.ones(sl.shape)
        for k in range(1, 512):
            term = term * sl / (al + k); total = total + term
            if k % 16 == 0 and float(np.max(term / np.maximum(total, 1e-300))) < 1e-17: break
        logP = al * np.log(np.maximum(sl, 1e-300)) - sl - gammaln_np(al + 1.0) + np.log(total)
        out[low] = np.log1p(-np.minimum(np.exp(logP), 1.0 - 1e-17))
    hi = ~low
    if hi.any():                                                    # Q = s^a e^-s / Gamma(a) * CF (Lentz)
        sh, ah = s[hi], a[hi]
        tiny = 1e-300
        b = sh + 1.0 - ah; c = np.full(sh.shape, 1.0 / tiny); d = 1.0 / np.maximum(b, tiny); h = d.copy()
        for i in range(1, 256):
            an = -i * (i - ah); b = b + 2.0
            d = an * d + b; d = np.where(np.abs(d) < tiny, tiny, d)
            c = b + an / c; c = np.where(np.abs(c) < tiny, tiny, c)
            d = 1.0 / d; delta = d * c; h = h * delta
            if i % 16 == 0 and float(np.max(np.abs(delta - 1.0))) < 1e-16: break
        out[hi] = ah * np.log(sh) - sh - gammaln_np(ah) + np.log(np.maximum(h, tiny))
    z = np.broadcast_arrays(x, out)[0]
    return np.where(z <= 0.0, 0.0, np.minimum(out, 0.0))


def bs_null(x, k, chunk=200_000, threads=1):
    """log p of each row's variance ratio x (its residual variance over the locus median, k comparators) under the
    scale null: rows differ in scale, x = sigma2 chi2_(k-1)/(k-1) with log sigma2 ~ Normal(mu, tau2), mu and tau2 from
    the bulk of all rows dataset-wide (robust location and spread of log x less the chi-square's own part), the
    integral by Gauss-Hermite. Returns (log p, the fitted SD tau)"""
    x = np.asarray(x, np.float64); k = np.asarray(k, np.float64)
    nu = np.maximum(k - 1.0, 1.0)
    u = np.log(np.maximum(x, 1e-12)) - (digamma_np(nu / 2.0) + math.log(2.0) - np.log(nu))
    mu = float(np.median(u)); tau2 = max((1.4826 * float(np.median(np.abs(u - mu)))) ** 2 - float(np.median(trigamma_np(nu / 2.0))), 1e-6)
    t, w = np.polynomial.hermite.hermgauss(40); s_ = mu + math.sqrt(2.0 * tau2) * t
    logw = np.log(w) - 0.5 * math.log(math.pi); es = np.exp(-s_)
    out = np.empty(len(x))
    def do(i0):
        sl = slice(i0, min(i0 + chunk, len(x)))
        terms = logw[None, :] + chi2_logsf_np(nu[sl, None] * x[sl, None] * es[None, :], nu[sl, None] * np.ones((1, len(s_))))
        m_ = terms.max(1)
        return i0, m_ + np.log(np.exp(terms - m_[:, None]).sum(1))
    starts = list(range(0, len(x), chunk))
    if threads > 1 and len(starts) > 1 and _can_fork():                    # rows are independent, so any split gives the same bits
        chunk = max(20_000, min(chunk, -(-len(x) // (4 * threads)))); starts = list(range(0, len(x), chunk))
        for i0, v in fork_map(do, starts, min(threads, len(starts))): out[i0:i0 + len(v)] = v
    else:
        for i0 in starts: _, v = do(i0); out[i0:i0 + len(v)] = v
    return np.minimum(out, 0.0), math.sqrt(tau2)


def order_stat_bound(lp_sorted_k, nf, k):
    """log of the order-statistic bound P(p_(k) <= p) <= C(nf, k) p^k, from the k-th smallest log tail; capped at 0.
    k = 1 is the Bonferroni min-tail; the excess side uses k = ceil(nf/2), so half the far comparators must lean"""
    nf = np.asarray(nf, np.float64); k = np.asarray(k, np.float64)
    logc = gammaln_np(nf + 1.0) - gammaln_np(k + 1.0) - gammaln_np(np.maximum(nf - k, 0.0) + 1.0)
    return np.where(nf >= np.maximum(k, 1), np.minimum(k * lp_sorted_k + logc, 0.0), 0.0)


# ---------------------------------------------------------------- the broad-sample draws
# (our own arithmetic end to end: the pools and draws depend only on the taxon and locus indices and these constants,
# never on numpy's generator algorithms, threads or block sizes, so the same input gives the same comparators anywhere)

def _sm64(z):
    """splitmix64, vectorised over uint64"""
    z = (np.asarray(z, np.uint64) + np.uint64(0x9E3779B97F4A7C15))
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return z ^ (z >> np.uint64(31))


def _first_distinct(vals, need):
    """per row the first `need` distinct values in slot order, -1 padded, and whether each row has enough: the draws
    are hashes, so a duplicate is simply skipped, exactly as a rejection draw would"""
    R, S = vals.shape
    o = np.argsort(vals, 1, kind="stable")
    sv = np.take_along_axis(vals, o, 1)
    dup = np.zeros((R, S), bool); dup[:, 1:] = sv[:, 1:] == sv[:, :-1]
    keep = np.ones((R, S), bool); np.put_along_axis(keep, o, ~dup, 1)       # duplicates keep their first slot
    cnt = keep.sum(1)
    first = np.argsort(~keep, 1, kind="stable")[:, :need]                   # kept slots first, in slot order
    out = np.take_along_axis(vals, first, 1)
    out[np.arange(need)[None, :] >= cnt[:, None]] = -1
    return out, cnt >= need


def _hash_draw(keys, m, need, start=0):
    """per row of keys: `need` distinct values in [0, m) by hashing (key, slot) from slot `start`, extending the
    slots deterministically until every row has enough"""
    R = len(keys)
    if m < 2 * need:                                                         # the pool wants most candidates (few taxa): rank them all by hash
        h = _sm64(keys[:, None] ^ np.arange(m, dtype=np.uint64)[None, :])
        return np.argsort(h, 1, kind="stable")[:, :need].astype(np.int64)
    S = int(need + 4 * math.sqrt(need) + (need * need) / max(m - need, 1) + 16)
    out = np.full((R, need), -1, np.int64); todo = np.arange(R)
    while len(todo):
        slots = (np.arange(start, start + S, dtype=np.uint64))[None, :]
        vals = (_sm64(keys[todo][:, None] ^ slots) % np.uint64(m)).astype(np.int64)
        got, ok = _first_distinct(vals, need)
        out[todo[ok]] = got[ok]
        todo = todo[~ok]; start += S; S *= 2                                 # the unlucky rows draw further slots; same answer whatever the batch
        if start > (1 << 20): raise RuntimeError(f"broad-sample draw could not find {need} distinct of {m}")
    return out


def pair_median(vals, W, min_loci, fallback):
    """site-weighted median of vals across loci (axis 0), each locus weighted by the sites the
    pair shares in it. Pairs with fewer than min_loci contributing loci get `fallback`; also
    returns the count."""
    ok = np.isfinite(vals) & (W > 0)
    nl = ok.sum(0)
    v = wmedian_stack(vals, np.where(ok, W, 0.0))
    return np.where(nl >= min_loci, v, fallback), nl


# ------------------------------------------------------------------ the pair stage's compiled kernel (numba, optional)
# The site-weighted medians of every pair over its loci, written as plain Python that numba compiles on first use when
# it is installed and the stage is large enough to earn its import (NUMBA_PAIR_LOCI_MIN). Operation for operation the
# arithmetic is that of the numpy path (pair_divergences, pair_median, wmedian_stack): float32 divergences, float64 once
# the rate divides them, the integer site counts as the weights; the two paths give the same bits and the suite holds
# them to it (tests/unit/test_kernel.py). The median is found by selection, not a sort, and nothing is materialised, so
# a core does a pair-locus in 35 ns against numpy's 107 and 96 threads stay busy where numpy's chunks fight for memory
# (Zuntini's pair stages 183 s to about 22 s). prange is range until _pair_kernel rebinds it to numba's.
prange = range
_KERNEL = None          # set by _pair_kernel: (the compiled kernel, a thread setter), or False when numba is not installed


def _wmed_select_py(v, w, L):
    """the weighted median of v[:L] under the integer weights w[:L], by selection: the value at the first position of the
    sorted order whose prefix weight reaches half the total, which is what wmedian_stack returns; exact because prefix
    sums of integers are the same in any order. v and w are scratch and come back permuted. NaN when nothing weighs."""
    tot = 0.0
    for t in range(L): tot += w[t]
    if tot <= 0: return np.nan
    half = tot / 2; lo = 0; hi = L - 1; wleft = 0.0            # [lo, hi] still to search; wleft: the weight placed before lo
    while hi - lo >= 16:
        mid = (lo + hi) >> 1; a = v[lo]; b = v[mid]; c = v[hi]  # the median of three as the pivot
        if a > b: a, b = b, a
        if b > c: b = c
        if a > b: b = a
        piv = b; i = lo; j = lo; k = hi; wl = 0.0; we = 0.0    # three-way partition: [lo, i) below, [i, j) equal, (k, hi] above
        while j <= k:
            x = v[j]
            if x < piv:
                v[i], v[j] = v[j], v[i]; w[i], w[j] = w[j], w[i]; wl += w[i]; i += 1; j += 1
            elif x > piv:
                v[j], v[k] = v[k], v[j]; w[j], w[k] = w[k], w[j]; k -= 1
            else:
                we += w[j]; j += 1
        if wleft + wl >= half: hi = i - 1                        # the half-way position is below the pivot
        elif wleft + wl + we >= half: return piv                 # ... among the values equal to it
        else: wleft += wl + we; lo = j                           # ... above it
    for s in range(lo + 1, hi + 1):                              # the last few: an insertion sort, then walk the prefix weights
        x = v[s]; y = w[s]; t = s - 1
        while t >= lo and v[t] > x:
            v[t + 1] = v[t]; w[t + 1] = w[t]; t -= 1
        v[t + 1] = x; w[t + 1] = y
    cw = wleft
    for t in range(lo, hi + 1):
        cw += w[t]
        if cw >= half: return v[t]
    return v[hi]


_wmed_select = _wmed_select_py


def _pair_medians_py(S, W, c0, c1, rr, ok_r, diag, min_loci, c50, B):
    """pairs c0..c1 of the per-locus counts S (same residues) and W (shared sites), both (loci, pairs): the site-weighted
    median over loci of the rate-normalised divergence 100 (1 - S/W) / rr, the same median of the floor c50 / W / rr,
    and the number of loci voting, each locus weighted by W, zero on the diagonal and at loci without a rate (ok_r);
    NaN below min_loci. The pairs are gathered B at a time so that a pair's loci are read from one slab of rows."""
    nloc = S.shape[0]; npair = c1 - c0
    nblk = (npair + B - 1) // B
    m = np.full(npair, np.nan); m_fl = np.full(npair, np.nan); nl = np.zeros(npair, np.int64)
    one = np.float32(1.0); hundred = np.float32(100.0); c50f = np.float32(c50)
    for b in prange(nblk):
        p0 = b * B; p1 = min(p0 + B, npair); nb_ = p1 - p0
        V = np.empty((nb_, nloc)); WV = np.empty((nb_, nloc)); cnt = np.zeros(nb_, np.int64)
        F = np.empty((nb_, nloc)); WF = np.empty((nb_, nloc)); cntf = np.zeros(nb_, np.int64)
        for l in range(nloc):
            if not ok_r[l]: continue
            r = rr[l]
            for q in range(nb_):
                p = c0 + p0 + q
                if diag[p]: continue
                wc = np.float32(W[l, p])
                if wc > 0:
                    sc = np.float32(S[l, p])
                    v = np.float64(hundred * (one - sc / wc)) / r
                    if np.isfinite(v):
                        V[q, cnt[q]] = v; WV[q, cnt[q]] = np.float64(wc); cnt[q] += 1
                    vf = np.float64(c50f / wc) / r
                    if np.isfinite(vf):
                        F[q, cntf[q]] = vf; WF[q, cntf[q]] = np.float64(wc); cntf[q] += 1
        for q in range(nb_):
            nl[p0 + q] = cnt[q]
            if cnt[q] >= min_loci and cnt[q] > 0: m[p0 + q] = _wmed_select(V[q], WV[q], cnt[q])
            if cntf[q] >= min_loci and cntf[q] > 0: m_fl[p0 + q] = _wmed_select(F[q], WF[q], cntf[q])
    return m, m_fl, nl


_pair_medians = _pair_medians_py
_pair_medians_serial = _pair_medians_py


def _renamed(f, name):
    """f under another name. numba keeps one cache file per function name, whatever the options: the parallel and the
    serial compilations of one function would share it, and whichever was compiled first would be loaded for both (the
    sweep's worker processes each running a thread pool, or the parent's kernel running on one thread)"""
    g = types.FunctionType(f.__code__.replace(co_name=name), f.__globals__, name, f.__defaults__, f.__closure__); g.__qualname__ = name; return g


def _pair_kernel(max_threads=1):
    """the compiled pair kernel and a thread setter, built on first use: (kernel, set_threads), or None when numba is not
    installed or --no-numba was given. numba is imported here rather than at the top: it costs about half a second, which
    only a pair stage of NUMBA_PAIR_LOCI_MIN loci x pairs earns back. Its thread pool is capped at max_threads (the
    run's --threads) and uses the workqueue layer, the one that is safe beside the fork pools of the later stages."""
    global _KERNEL, _wmed_select, _pair_medians, _pair_medians_serial, _closest_rows_nb, prange
    if not USE_NUMBA: return None
    if _KERNEL is None:
        os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue"); os.environ.setdefault("NUMBA_NUM_THREADS", str(max(1, int(max_threads))))
        try: import numba
        except ImportError:
            _KERNEL = False; LOG("pair stage: numba is not installed, so the numpy pair stage runs (conda install numba for the compiled one: the same numbers, 3-8x faster at scale)")
        else:
            prange = numba.prange
            _wmed_select = numba.njit(cache=True)(_wmed_select_py)
            _pair_medians = numba.njit(parallel=True, cache=True)(_pair_medians_py)
            _pair_medians_serial = numba.njit(cache=True, error_model="numpy")(_renamed(_pair_medians_py, "_pair_medians_serial_py"))   # for worker processes: the same bits, no thread pool
            _closest_rows_nb = numba.njit(cache=True)(_closest_rows_py)         # the closest-match reduce of the sweep's block pairs
            def set_threads(t):
                t = max(1, min(int(t), int(numba.config.NUMBA_NUM_THREADS))); numba.set_num_threads(t); return t
            _KERNEL = (_pair_medians, set_threads)
    return _KERNEL or None


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
MEM_RESIDENT_EST = 0              # what the run will hold resident once its arrays are filled (main's estimate): the room is measured
                                  # against the larger of this and the resident set now, since empty arrays are not resident yet


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
    return max(MEM_BUDGET - max(_rss(), MEM_RESIDENT_EST), MEM_BUDGET // 16)


def fmem(b):
    """bytes as GB or MB"""
    return f"{b / float(1 << 30):.1f} GB" if b >= (1 << 30) else f"{b / float(1 << 20):.0f} MB"


def memory_estimate(n, L, loci, C=None, full=None, K_eff=5, rc=True):
    """bytes the run holds resident, from the alignment's shape: the alignment twice (the raw rows and the residue
    codes), each taxon's candidate list (C entries of 32 bytes), the closest-match band entries of every cell (24
    bytes per slot), up to FULL_M_MAX_TAXA taxa the full m matrix and its floor, up to SWEEP_BLOCK_MIN taxa the
    counts of every pair at every locus (the one block pair's slab, kept for the run), and the broad-sample
    pools with their pair divergences (BS_POOL x K_eff per taxon, 12 bytes each). The sweep's per-worker slabs
    and the stages' other temporaries are sized to the room left and are not counted. Returns the total and the six parts."""
    nloc = len(loci); C = min(C or CANDIDATES, n); full = n <= FULL_M_MAX_TAXA if full is None else full
    longest = max((len(st) for _, st in loci), default=0); slab = 2 * nloc * n * n * (2 if longest < 65536 else 4) if n <= SWEEP_BLOCK_MIN else 0
    x = 2 * n * L; lists = n * C * 32; closest = n * nloc * BAND_SLOTS * 24; mfull = 2 * n * n * 8 if full else 0
    mrc = n * min(BS_POOL * K_eff, max(n - 1, 0)) * 12 if rc else 0
    return x + lists + closest + mfull + slab + mrc, x, lists, closest, mfull, slab, mrc



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

def report_floor_bytes(n, nloc, K_eff):
    """the least the report bundle can weigh before compression: the per-cell members every run writes (five bytes per
    cell, and up to DENSE_CELLS_MAX cells the dense numeric columns and the comparators' tests); deflate leaves about a third"""
    cells = n * nloc; b = cells * 5
    if cells <= DENSE_CELLS_MAX: b += cells * (14 * 4 + 3 * K_eff * 4)
    return b


def write_tables(a, names, loci, rep, taxa_json, loci_json, manifest, K_eff):
    """the viewer-free tables, written when report.html is not: taxa.tsv (status, loci, masks, comparator shortfalls,
    the closest genomes), loci.tsv (length, status, rate, plateau, dispersion, masks), summary.json (= manifest.json)"""
    n, nloc = len(names), len(loci); side = rep["side"]; n_comp = rep["n_comp"]; tested = n_comp > 0
    div = (side & 1) > 0; sim = (side & 2) > 0; short = tested & (n_comp < K_eff); cl = np.asarray(rep["taxa_closest"])
    with open(f"{a.out}/taxa.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["taxon", "status", "loci_present", "loci_tested", "masked_divergent", "masked_similar", "short_of_comparators"] + [f"closest_{k + 1}" for k in range(cl.shape[1])]) + "\n")
        for i in range(n):
            fh.write("\t".join([names[i], taxa_json["status"][i], str(taxa_json["loci_present"][i]), str(int(tested[i].sum())), str(int(div[i].sum())), str(int(sim[i].sum())), str(int(short[i].sum()))]
                               + [names[int(c)] if c >= 0 else "" for c in cl[i]]) + "\n")
    with open(f"{a.out}/loci.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["locus", "length", "status", "rate", "plateau", "phi", "test_rate", "similar_tests", "cells_tested", "masked_divergent", "masked_similar"]) + "\n")
        for li in range(nloc):
            fh.write("\t".join([loci_json["names"][li], str(loci_json["length"][li]), loci_json["status"][li]] + ["" if v is None else str(v) for v in (loci_json["rate"][li], loci_json["plateau"][li], loci_json["phi"][li], loci_json["test_rate"][li])]
                               + [str(loci_json["similar_tests"][li]), str(int(tested[:, li].sum())), str(int(div[:, li].sum())), str(int(sim[:, li].sum()))]) + "\n")
    with open(f"{a.out}/summary.json", "w", encoding="utf-8") as fh: json.dump(manifest, fh, indent=1)


def write_report(a, names, loci, rep, dtype, fmt, K, K_eff, M_eff, t0):
    """report.npz, a ZIP of .npy members and JSON that np.load opens like an .npz (streamed to the file), and
    report.html, the viewer with the bundle inline as base64, one file, up to REPORT_ONE_FILE_MAX of bundle; above
    that the tables stand in for the page. Both files are written as .part and renamed, so a file is whole or absent.
    Returns True when report.html was written"""
    n, nloc = len(names), len(loci)
    dense = rep["dense"]
    npz = f"{a.out}/report.npz"
    fh_npz = open(npz + ".part", "wb")
    with zipfile.ZipFile(fh_npz, "w", zipfile.ZIP_DEFLATED, compresslevel=4) as zf:          # level 4: half the time of 6 for 5% more bytes
        def put(name, arr):
            b = io.BytesIO(); np.lib.format.write_array(b, np.ascontiguousarray(arr), version=(1, 0)); zf.writestr(name, b.getvalue())
        def put_json(name, obj, compact=False): zf.writestr(name, json.dumps(obj, separators=(",", ":")) if compact else json.dumps(obj, indent=1))
        for k in ("n_comp", "n_fail", "n_similar", "side", "focal_sites"): put(f"cells/{k}.npy", rep[k])
        if "side2" in rep:
            put("cells/side2.npy", rep["side2"])
            for k in ("bs_locus", "bs_taxon", "bs_var_log10p", "bs_lo_log10p", "bs_hi_log10p", "bs_partner", "bs_call", "bs_log10p", "bs_rho", "bs_stat", "bs_rep"):
                if k in rep: put(f"cells/{k}.npy", rep[k])
        numcols = ("med_identity", "med_expected", "excess_m", "closest", "closest_overlap", "closest_identity", "closest_expected",
                   "best_match", "best_match_identity", "partner", "div_log10p", "sim_log10p", "locus_closest", "locus_closest_log10p")
        if not dense: put("cells/sparse_locus.npy", rep["sparse_locus"]); put("cells/sparse_taxon.npy", rep["sparse_taxon"])
        for k in numcols: put(f"cells/{k}.npy", rep[k])
        for k in ("comparators", "comp_div_log10p", "comp_sim_log10p"): put(f"cells/{k}.npy", rep[k])
        for li, (cells_, blk) in rep["codes"].items():
            put(f"codes/{li}.npy", blk); put(f"codes/{li}_cells.npy", cells_)
        for li, tr in rep["trees"].items(): put_json(f"trees/{li}.json", tr, compact=True)         # machine-read only: no indentation
        for li, ev in rep["similar"].items(): put_json(f"similar/{li}.json", ev, compact=True)
        for li, ev in rep.get("bs_records", {}).items(): put_json(f"broadsample/{li}.json", ev, compact=True)
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
        taxa_json = dict(
            names=list(names),
            status=[("excluded: %d loci with at least %d residues, %d needed" % (int(rep["present"][i]), a.fit_overlap, MIN_LOCI_PER_TAXON)) if i in thin else "tested" for i in range(n)],
            loci_present=[int(v) for v in (rep["focal_sites"] > 0).sum(1)])
        put_json("taxa.json", taxa_json)
        loci_json = dict(
            names=[nm for nm, _ in loci], length=[int(len(st)) for _, st in loci], ranges=[_ranges(st) for _, st in loci],
            status=list(rep["loc_status"]), rate=[None if not np.isfinite(v) else round(float(v), 5) for v in rep["rate"]],
            plateau=[round(float(v), 3) for v in rep["plateaus"]], phi=[None if not np.isfinite(v) else float(v) for v in rep["phi"]],
            test_rate=[None if not np.isfinite(v) else round(float(v), 5) for v in rep["test_rate"]], similar_tests=[int(v) for v in rep["pair_tests"]])
        put_json("loci.json", loci_json)
        manifest = dict(
            format_version=REPORT_FORMAT, tool="LocusMasker", version=__version__, commit=_git_commit(),
            python=platform.python_version(), numpy=np.__version__,
            system=platform.platform(), date=time.strftime("%Y-%m-%d %H:%M:%S"), command=" ".join(sys.argv),
            options=dict(datatype=dtype, n_comp=K, K_eff=K_eff, M_eff=M_eff, fit_overlap=a.fit_overlap, alpha=a.alpha, donor_threshold=a.donor_threshold,
                         mask_similar_only=a.mask_similar_only, mask_different_only=a.mask_different_only, threads=a.threads, all_tests=a.all_tests,
                         report_codes=a.report_codes, report_sequences=rep["seq_policy"]),
            input=dict(file=os.path.basename(a.aln), bytes=os.path.getsize(a.aln), sha256=getattr(a, "input_sha256", None) or _sha256(a.aln), format=fmt,
                       partitions=os.path.basename(a.partitions) if a.partitions else None),
            dataset=a.name or (os.path.basename(os.path.dirname(os.path.abspath(a.aln))) if os.path.splitext(os.path.basename(a.aln))[0].lower() == "alignment" and os.path.basename(os.path.dirname(os.path.abspath(a.aln))) not in ("", ".", "work") else os.path.splitext(os.path.basename(a.aln))[0]),
            counts=counts, test=dict(tests=int(rep["tests"]["total"]), close_relative_tests=int(rep["tests"]["close_relative"]),
                                     broad_sample_tests=int(rep["tests"].get("broad_sample", 0)),
                                     closest_match_relative_tests=int(rep["tests"]["closest_match_relative"]), closest_match_outside_tests=int(rep["tests"]["closest_match_outside"]),
                                     similar_candidates=int(rep["tests"]["locus_sets"]), comparisons=int(rep["tests"]["close_relative"]),
                                     bonferroni_k=round(float(rep["kq"]), 4), log_bar=round(float(rep["logbar"]), 4), log10_bar=round(float(rep["logbar"]) / LN10, 4),
                                     phi_dataset=round(float(rep["phi_bar"]), 6), half=HALF),
            carried=dict(numeric_columns="dense" if dense else "sparse", codes=a.report_codes, sequences=rep["seq_policy"], trees="masked",
                         tree_partners="the pair's partner", similar="every too-similar test and every source, per locus (similar/<locus>.json)",
                         tests="every cell's comparators and both of each comparator's tests (cells/comp_div_log10p, cells/comp_sim_log10p, log10 p by comparator rank), the closest match at the locus and its too-similar test (cells/locus_closest, cells/locus_closest_log10p)",
                         best_match_cells="all" if n <= BEST_MATCH_ALL_MAX_TAXA else "fails"),
            code_map=dict(mismatch_count_max=K_eff, similar_base=CODE_SIMILAR, tick=CODE_TICK, absent=CODE_ABSENT),
            side_bits=dict(masked_divergent=1, masked_similar=2, ambiguous=4, similar_uncertain=8, flagged_not_masked=16, similar_not_masked=32, similar_tested=64, similar_source=128),
            side2_bits=dict(broad_tested=1, broad_past_bar=2, broad_close=4, broad_divergent=8, broad_excluded_sample=16, broad_fired_not_masked=32, divergent_also_close_label=64),
            broad_sample=rep.get("bs_summary"),
            restart=rep.get("restart"), lists_extended=rep.get("lists_extended"),
            timings=dict(PROG.durations, total=round(time.time() - t0, 2)))
        put_json("manifest.json", manifest)
    fh_npz.close(); os.replace(npz + ".part", npz); size = os.path.getsize(npz)
    if size > REPORT_ONE_FILE_MAX:
        with zipfile.ZipFile(npz) as Z: extra = sum(i.compress_size for i in Z.infolist() if i.filename.startswith(("codes/", "seq/")))
        write_tables(a, names, loci, rep, taxa_json, loci_json, manifest, K_eff)
        LOG(f"report.html not written: the report bundle is {fmem(size)}, above the one-file ceiling of {fmem(REPORT_ONE_FILE_MAX)} (a browser holds about four times the bundle). "
            f"report.npz holds every result (numpy.load opens it), and taxa.tsv, loci.tsv and summary.json carry the per-taxon and per-locus tables. "
            f"For a report.html: --report-codes none --report-sequences none drop the site codes and sequences ({fmem(extra)} of the bundle), or run a subset of the taxa or loci")
        return False
    tpl_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viewer"); tpl = os.path.join(tpl_dir, "viewer.html")
    if not os.path.exists(tpl):
        LOG(f"report.html not written: the viewer template is not beside the script ({tpl}); report.npz holds the data"); return False
    def text(name):
        with open(os.path.join(tpl_dir, name), encoding="utf-8") as fh: return fh.read()
    def inline_script(mt): return "<script>\n" + text(mt.group(1)) + "\n</script>"
    def inline_style(mt): return "<style>\n" + text(mt.group(1)) + "\n</style>"
    def inline_icon(mt):                                           # the logo and favicons, as data: URIs
        src = os.path.join(tpl_dir, mt.group(2)); kind = "image/svg+xml" if src.endswith(".svg") else "image/png"
        with open(src, "rb") as fh: return f'{mt.group(1)}="data:{kind};base64,' + base64.b64encode(fh.read()).decode("ascii") + '"'
    def inlined(part):
        part = re.sub(r'<script src="([^"]+)"></script>', inline_script, part)
        part = re.sub(r'<link rel="stylesheet" href="([^"]+)">', inline_style, part)
        return re.sub(r'(href|src)="(icons/[^"]+)"', inline_icon, part)
    before, after = text("viewer.html").split("@@BUNDLE@@", 1)      # the template is split at the marker first, so nothing inlined can carry one
    meta_js = json.dumps(dict(dataset=manifest["dataset"], bundle_bytes=size, version=__version__))
    before = before.replace("@@META@@", meta_js, 1)                  # the loading screen's title and total, known before the data arrives
    page = f"{a.out}/report.html"
    with open(page + ".part", "w", encoding="utf-8") as out, open(npz, "rb") as src:   # the bundle as base64 in 3-byte multiples: no copy of it in memory
        out.write(inlined(before))
        while True:
            chunk = src.read(3 << 20)
            if not chunk: break
            out.write(base64.b64encode(chunk).decode("ascii"))
        out.write(inlined(after))
    os.replace(page + ".part", page)
    LOG(f"report: report.html ({os.path.getsize(page) >> 20} MB), report.npz ({size >> 20} MB)")
    return True


def _slab_chunk_medians(Sc, Wc, rr, ok_r, diag_local):
    """m, m_fl and the loci voting for a chunk of pair positions of a (loci, positions) slab of counts: the arithmetic of
    the numpy pair stage, which the compiled kernel reproduces bit for bit"""
    Sc = Sc.astype(np.float32); Wc = Wc.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        Dv = np.where(Wc > 0, 100.0 * (1.0 - Sc / np.maximum(Wc, 1)), np.nan)
    Wfl = 100.0 * HALF / np.maximum(Wc, 1.0)
    Wm = Wc.copy(); Wm[:, diag_local] = 0.0
    if rr is not None:
        Dv = Dv / rr[:, None]; Wfl = Wfl / rr[:, None]
        Wm = np.where(ok_r[:, None], Wm, 0.0)
    nan_fb = np.full(Sc.shape[1], np.nan)
    m, nl = pair_median(Dv, Wm, MIN_LOCI_FOR_M, nan_fb); m_fl, _ = pair_median(Wfl, Wm, MIN_LOCI_FOR_M, nan_fb)
    return m, m_fl, nl


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
    kern = _pair_kernel(threads) if nloc * npair >= NUMBA_PAIR_LOCI_MIN else None
    if kern is not None:
        kernel, set_threads = kern; used = set_threads(threads)
        rr_ = np.ones(nloc) if rr is None else np.asarray(rr, np.float64); ok_ = np.ones(nloc, bool) if ok_r is None else np.asarray(ok_r, bool)
        is_diag = iu[0] == iu[1]; starts = range(0, npair, NUMBA_PAIR_SLAB)
        LOG(f"pair stage: compiled kernel (numba), {len(starts)} slab{'s' if len(starts) > 1 else ''} of {min(NUMBA_PAIR_SLAB, npair):,} pairs, {used} thread{'s' if used > 1 else ''}")
        if progress: progress.start(stage, len(starts))
        for c0 in starts:
            c1 = min(c0 + NUMBA_PAIR_SLAB, npair)
            m[c0:c1], m_fl[c0:c1], nl[c0:c1] = kernel(S, W, c0, c1, rr_, ok_, is_diag, MIN_LOCI_FOR_M, 100.0 * HALF, PAIR_BLOCK)
            if progress: progress.tick()
        if progress: progress.finish()
        return _expand(m, n, iu), _expand(m_fl, n, iu), _expand(nl, n, iu)
    per_chunk = min(PAIR_CHUNK_MAX, max(PAIR_CHUNK_MIN, _mem_room() // (4 * threads)))
    chunk = max(int(per_chunk // (64 * max(nloc, 1))), 1)
    if threads > 1: chunk = max(min(chunk, -(-npair // (4 * threads))), int((32 << 20) // (64 * max(nloc, 1))), 1)   # four chunks per thread, but never below 32 MB
    chunk = min(chunk, npair)
    threads = max(1, min(threads, -(-npair // chunk)))                       # no more workers than chunks
    LOG(f"pair stage: {-(-npair // chunk)} chunk{'s' if npair > chunk else ''} of {chunk} pairs, about {max(chunk * 64 * nloc >> 20, 1)} MB each, {threads} thread{'s' if threads > 1 else ''}")
    def do_chunk(c0):
        c1 = min(c0 + chunk, npair)
        m[c0:c1], m_fl[c0:c1], nl[c0:c1] = _slab_chunk_medians(S[:, c0:c1], W[:, c0:c1], rr, ok_r, diag[(diag >= c0) & (diag < c1)] - c0)
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



# ------------------------------------------------------------------ the pair sweep: nothing is stored per pair
# Every statement of the method is about pairs of taxa, and the old engine stored every pair's per-locus counts and ten
# taxa x taxa matrices (142 GB at 9,524 taxa; terabytes at GTDB size). The sweep walks pairs of taxon blocks: for each
# block pair it computes the per-locus counts with pass one's matrix products, the medians with the kernel, the
# identity, the ranking key and, in sweep 2, the closest match at every locus, keeps what the method asks for (each
# taxon's closest CANDIDATES by the key, the closest match per cell, the identical pairs, a few counters, and the full m
# up to FULL_M_MAX_TAXA) and drops the counts. Everything downstream reads a cell's few pairs from the residues
# (pair_counts_at) and a pair's m from the lists. The reductions are order-independent (the keys are total orders), so
# any block size, block order and worker count give the same bits; the goldens hold the sweep to the old engine.
_SWEEP = None                     # the sweep's context for worker processes (set before the pool is forked)


def encode_alignment(names, rows, dtype):
    """the alignment as residue codes (taxa x sites int8, -1 for everything outside the alphabet); warns when most
    letters are ambiguity codes"""
    n, L = len(names), len(rows[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = np.empty((n, L), dtype=np.int8)                               # one row at a time: no joined copy of the alignment
    letters = np.zeros(256, np.int64)
    for i, r in enumerate(rows): b = np.frombuffer(r, dtype=np.uint8); X[i] = lut[b]; letters += np.bincount(b, minlength=256)
    is_letter = np.zeros(256, bool); is_letter[65:91] = True; is_letter[97:123] = True
    n_letters = int(letters[is_letter].sum()); n_in = int(letters[is_letter & (lut >= 0)].sum())
    if n_letters and n_in < 0.95 * n_letters:
        LOG(f"WARNING: {n_letters - n_in:,} of the {n_letters:,} letters ({100.0 * (n_letters - n_in) / n_letters:.1f}%) are ambiguity codes and count as missing data")
    return X


PLATEAU_BIG = 1 << 20             # a locus of this many residues x taxa or more counts its letters one at a time (no cast to a wider type) and
                                  # the loci run in threads; smaller loci count in one bincount call each, serially (the per-call overhead rules there)
def locus_plateaus(X, loci, A_, threads=1):
    """each locus's random-identity plateau from its residue frequencies over every taxon, thin ones included (as pass
    one computed it). The same integer counts either way, so the same plateaus"""
    plateaus = np.empty(len(loci)); n = X.shape[0]
    def one(li):
        Xm = X[:, loci[li][1]]
        if Xm.size < PLATEAU_BIG: v = Xm.ravel(); cnt = np.bincount(v[v >= 0], minlength=A_).astype(float)
        else: cnt = np.array([np.count_nonzero(Xm == a) for a in range(A_)], float)
        f = cnt / max(cnt.sum(), 1.0); plateaus[li] = max(100.0 * (1.0 - float((f ** 2).sum())), 1.0)
    big = len(loci) > 1 and float(np.mean([len(st) for _, st in loci])) * n >= PLATEAU_BIG
    if threads > 1 and big:
        with ThreadPoolExecutor(min(threads, len(loci))) as ex: list(ex.map(one, range(len(loci)), chunksize=max(1, len(loci) // (8 * threads))))
    else:
        for li in range(len(loci)): one(li)
    return plateaus


def _rows_sites(X, idx, sites):
    return X[idx][:, sites] if isinstance(idx, slice) else X[np.ix_(idx, sites)]


def _letter_columns(XA, XB, A_):
    """for each letter the columns where it occurs in either block: the other columns add nothing to its product"""
    if A_ > 4:
        occ = np.zeros((A_ + 1, XA.shape[1]), bool); ar = np.arange(XA.shape[1]); occ[XA + 1, ar] = True
        if XB is not XA: occ[XB + 1, ar] = True
        return [np.flatnonzero(occ[i + 1]) for i in range(A_)]
    return [np.flatnonzero((XA == i).any(0) | (XB == i).any(0)) for i in range(A_)]


def _block_counts(X, sites, A, B, A_):
    """the shared sites and same residues of every (a in A, b in B) pair at one locus by pass one's matrix products,
    exact integers below 2^24 in two (bA, bB) float32 arrays. A and B are slices or index arrays"""
    XA = _rows_sites(X, A, sites); same = B is A; XB = XA if same else _rows_sites(X, B, sites)
    PA = (XA >= 0).astype(np.float32); PB = PA if same else (XB >= 0).astype(np.float32)
    sh = PA @ PB.T; sm = np.zeros(sh.shape, np.float32)
    for i, cols in enumerate(_letter_columns(XA, XB, A_)):
        if len(cols):
            Ai = (XA[:, cols] == i).astype(np.float32); Bi = Ai if same else (XB[:, cols] == i).astype(np.float32); sm += Ai @ Bi.T
    return sh, sm


def pair_counts_at(X, sites, gi, ci, keep, Xs=None):
    """shared sites and same residues of the pairs (gi[k], ci[k]) at a locus, int64; a pair with an excluded member
    counts 0 of each, as a thin taxon's stored counts were zeroed. Xs: the locus's residues, when the caller has them"""
    gi = np.asarray(gi, dtype=np.int64); ci = np.asarray(ci, dtype=np.int64); o = np.zeros(len(gi), np.int64); same = np.zeros(len(gi), np.int64)
    if len(gi) == 0: return o, same
    if Xs is None: Xs = X[:, sites]
    step = max(1, (10 << 20) // max(Xs.shape[1], 1))                           # about 10 M pair-sites of temporaries at a time
    for k0 in range(0, len(gi), step):
        XG = Xs[gi[k0:k0 + step]]; XC = Xs[ci[k0:k0 + step]]; both = (XG >= 0) & (XC >= 0)
        o[k0:k0 + step] = both.sum(1); same[k0:k0 + step] = (both & (XG == XC)).sum(1)
    ok = keep[gi] & keep[ci]; o[~ok] = 0; same[~ok] = 0
    return o, same


class LocusRows:
    """the one-hot products of one locus, built once, for rows of counts of a few taxa against every taxon (the best
    match, the outgroup search); pairs with an excluded member count 0"""
    def __init__(self, X, sites, keep, A_, Xs=None):
        self.Xs = X[:, sites] if Xs is None else Xs; self.keep = keep; self.P = (self.Xs >= 0).astype(np.float32)
        self.ones = [(i, c, (self.Xs[:, c] == i).astype(np.float32)) for i, c in enumerate(_letter_columns(self.Xs, self.Xs, A_)) if len(c)]

    def rows(self, rows, slice_rows=256):
        """shared sites and same residues of each row taxon against every taxon, two (rows, n) float64 arrays"""
        rows = np.asarray(rows, dtype=np.int64); n = self.Xs.shape[0]; W = np.empty((len(rows), n)); S = np.empty((len(rows), n))
        for r0 in range(0, len(rows), slice_rows):
            rr = rows[r0:r0 + slice_rows]; XR = self.Xs[rr]; sh = self.P[rr] @ self.P.T; sm = np.zeros(sh.shape, np.float32)
            for i, c, Ai in self.ones: sm += (XR[:, c] == i).astype(np.float32) @ Ai.T
            W[r0:r0 + slice_rows] = sh; S[r0:r0 + slice_rows] = sm
        W[~self.keep[rows]] = 0.0; S[~self.keep[rows]] = 0.0; W[:, ~self.keep] = 0.0; S[:, ~self.keep] = 0.0
        return W, S


def pair_rows_at(X, sites, rows, keep, A_, slice_rows=256):
    """shared sites and same residues of each row taxon against every taxon at a locus (LocusRows for one call)"""
    return LocusRows(X, sites, keep, A_).rows(rows, slice_rows)


class SlabCounts:
    """the pair counts of every locus, kept from the one block pair's sweep (up to SWEEP_BLOCK_MIN taxa, as the old
    engine stored them for the run): the counts of pairs and the rows of a few taxa are gathers. A pair's position in
    the slab is g * n + c; a pair with an excluded member holds 0"""
    def __init__(self, Ssl, Wsl, n): self.S = Ssl; self.W = Wsl; self.n = n

    def at(self, li, gi, ci, Xs=None):
        """shared sites and same residues of the pairs (gi[k], ci[k]) at locus li, int64"""
        idx = np.asarray(gi, np.int64) * self.n + np.asarray(ci, np.int64); return self.W[li][idx].astype(np.int64), self.S[li][idx].astype(np.int64)

    def locus(self, li, Xs=None): return _SlabRows(self.W[li].reshape(self.n, self.n), self.S[li].reshape(self.n, self.n))


class _SlabRows:
    """rows of a kept slab at one locus, with LocusRows' interface"""
    def __init__(self, W2, S2): self.W2 = W2; self.S2 = S2

    def rows(self, rows, slice_rows=256):
        r = np.asarray(rows, np.int64); return self.W2[r].astype(np.float64), self.S2[r].astype(np.float64)


class ResidueCounts:
    """the pair counts computed from the residues as they are needed (above SWEEP_BLOCK_MIN taxa, where no slab of
    every pair is kept): pairs by comparing their rows, the rows of a few taxa by the one-hot products"""
    def __init__(self, X, loci, keep, A_): self.X = X; self.loci = loci; self.keep = keep; self.A_ = A_

    def at(self, li, gi, ci, Xs=None): return pair_counts_at(self.X, self.loci[li][1], gi, ci, self.keep, Xs)

    def locus(self, li, Xs=None): return LocusRows(self.X, self.loci[li][1], self.keep, self.A_, Xs)


def slab_medians(Ssl, Wsl, diag, rr, ok_r, threads, kernel=None):
    """m, m_fl and the loci voting for every position of a (loci, positions) slab of counts: the compiled kernel when
    it was built, else the numpy chunks of the old pair stage (the same bits either way)"""
    nloc, P = Ssl.shape
    if kernel is not None:
        rr_ = np.ones(nloc) if rr is None else np.asarray(rr, np.float64); ok_ = np.ones(nloc, bool) if ok_r is None else np.asarray(ok_r, bool)
        return kernel(Ssl, Wsl, 0, P, rr_, ok_, diag, MIN_LOCI_FOR_M, 100.0 * HALF, PAIR_BLOCK)
    m = np.full(P, np.nan); m_fl = np.full(P, np.nan); nl = np.zeros(P, np.int64)
    chunk = max(int(min(PAIR_CHUNK_MAX, max(PAIR_CHUNK_MIN, _mem_room() // (4 * max(threads, 1)))) // (64 * max(nloc, 1))), 1)
    dpos = np.flatnonzero(diag)
    def do(c0):
        c1 = min(c0 + chunk, P)
        m[c0:c1], m_fl[c0:c1], nl[c0:c1] = _slab_chunk_medians(Ssl[:, c0:c1], Wsl[:, c0:c1], rr, ok_r, dpos[(dpos >= c0) & (dpos < c1)] - c0)
    starts = range(0, P, chunk)
    if threads > 1 and len(starts) > 1:
        with ThreadPoolExecutor(min(threads, len(starts))) as ex: list(ex.map(do, starts))
    else:
        for c0 in starts: do(c0)
    return m, m_fl, nl


def _block_keys(m, tot_same, tot_shared, gidx, cidx, diag):
    """for a block of pair positions: the identity flag and the comparator ranking key mj in both directions (the
    tie-break carries the OTHER taxon's pooled identity and index), by the operations the old engine used"""
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(tot_shared > 0, 100.0 * tot_same / tot_shared, np.nan)
    ani = np.where(diag, np.nan, ani)
    ident = (tot_shared > 0) & (tot_same == tot_shared) & ~diag
    base = np.nan_to_num(100.0 - ani, nan=100.0) * 1e-9
    fin = np.where(np.isfinite(m) & ~ident, m, np.inf)
    return ident, fin + (base + cidx * 1e-12), fin + (base + gidx * 1e-12)


def _topc_rows(keys, m, m_fl, cidx, C):
    """per row the C smallest finite keys in (key, column) order: columns (-1 beyond), m, m_fl, keys, and the number of
    finite keys in the row"""
    bA, bB = keys.shape; k = min(C, bB)
    order = np.lexsort((np.broadcast_to(cidx, (bA, bB)), keys), axis=1)[:, :k]
    sk = np.take_along_axis(keys, order, 1); fin = np.isfinite(sk)
    return (np.where(fin, cidx[order], -1), np.where(fin, np.take_along_axis(m, order, 1), np.nan),
            np.where(fin, np.take_along_axis(m_fl, order, 1), np.nan), np.where(fin, sk, np.inf), np.isfinite(keys).sum(1))


def _closest_rows(ov, sm, mj, fit, cidx):
    """the closest-match search of one locus for a block of rows over a block of columns (the old engine's lines, with
    bookkeeping): per row the minimum mismatch rate and, for each distinct rate within 1e-9 of it (BAND_SLOTS at
    most), the entry (d, vals, column) the rule would pick; overflow counts rows with more distinct rates than slots"""
    okc = (ov >= fit) & np.isfinite(mj)
    with np.errstate(invalid="ignore", divide="ignore"): d = np.where(okc, (ov - sm + HALF) / np.maximum(ov, 1.0), np.inf)
    dmin = d.min(1); R = BAND_SLOTS; nr = len(ov); big = np.iinfo(np.int64).max
    remaining = okc & (d <= dmin[:, None] + 1e-9); vals = d + 1e-9 * mj / 100.0
    ed = np.full((nr, R), np.inf); ev = np.full((nr, R), np.inf); ec = np.full((nr, R), -1, np.int64); cur = np.where(np.isfinite(dmin), dmin, np.inf)
    for k in range(R):
        grp = remaining & (d == cur[:, None])                           # the members at this distinct rate
        if not grp.any(): break
        v = np.where(grp, vals, np.inf); vmin = v.min(1); cmin = np.where(grp & (v == vmin[:, None]), cidx[None, :], big).min(1)
        ok = np.isfinite(cur); ed[ok, k] = cur[ok]; ev[ok, k] = vmin[ok]; ec[ok, k] = cmin[ok]
        remaining &= ~grp; cur = np.where(remaining, d, np.inf).min(1)  # the next distinct rate within the band
    return dmin, ed, ev, ec, int(np.isfinite(cur).sum())


def _closest_rows_py(ov, sm, mj, fit, cidx, R, half):
    """_closest_rows as one pass per row (compiled by _pair_kernel): the same operations in the same order, so the
    same bits. Per row: the rate of every column (inf where the pair is unusable), its minimum, then for each distinct
    rate within 1e-9 of it, R at most, the (vals, column)-smallest member; a row with a further distinct rate counts
    as an overflow"""
    nr, nc = ov.shape; dmin_out = np.empty(nr); ed = np.full((nr, R), np.inf); ev = np.full((nr, R), np.inf); ec = np.full((nr, R), -1, np.int64)
    dbuf = np.empty(nc); big = np.iinfo(np.int64).max; overflow = 0
    for i in range(nr):
        dmin = np.inf
        for j in range(nc):
            o = ov[i, j]
            if o >= fit and np.isfinite(mj[i, j]):
                dj = (o - sm[i, j] + half) / max(o, 1.0)
            else: dj = np.inf
            dbuf[j] = dj
            if dj < dmin: dmin = dj
        dmin_out[i] = dmin; cur = dmin; band = dmin + 1e-9; k = 0
        while k < R and cur < np.inf:
            vmin = np.inf; cmin = big
            for j in range(nc):
                if dbuf[j] == cur:
                    v = dbuf[j] + 1e-9 * mj[i, j] / 100.0
                    if v < vmin: vmin = v; cmin = cidx[j]
                    elif v == vmin and cidx[j] < cmin: cmin = cidx[j]
            ed[i, k] = cur; ev[i, k] = vmin; ec[i, k] = cmin
            nxt = np.inf
            for j in range(nc):
                dj = dbuf[j]
                if dj > cur and dj <= band and dj < nxt: nxt = dj
            cur = nxt; k += 1
        if cur < np.inf: overflow += 1
    return dmin_out, ed, ev, ec, overflow


_closest_rows_nb = None           # the compiled _closest_rows_py when numba is on (set by _pair_kernel), else None: the numpy rule runs


def closest_rows(ov, sm, mj, fit, cidx):
    """the closest-match reduce of a block pair at one locus: the compiled kernel when it was built, else the numpy
    rule; the same bits either way"""
    if _closest_rows_nb is not None: return _closest_rows_nb(ov, sm, mj, fit, cidx, BAND_SLOTS, HALF)
    return _closest_rows(ov, sm, mj, fit, cidx)


def _merge_closest(dd, dv, dc, sd, sv, sc):
    """merge two band-entry sets of the same cells (..., R): the entries within 1e-9 of the lower minimum, one per
    distinct rate (the (vals, column)-smallest), R at most; order-independent"""
    big = np.iinfo(np.int64).max; R = dd.shape[-1]
    if R > 1 and not (np.isfinite(dd[..., 1]).any() or np.isfinite(sd[..., 1]).any()):
        # every cell holds one entry on each side (the rule): the lower d wins, an equal d goes to the smaller (vals, column),
        # and only the rare cells with unequal d within 1e-9 of each other take the general path
        a, b = dd[..., 0], sd[..., 0]
        with np.errstate(invalid="ignore"): close = np.isfinite(a) & np.isfinite(b) & (a != b) & (np.abs(a - b) <= 1e-9)   # inf - inf where a cell has no candidate
        take_b = (b < a) | ((b == a) & ((sv[..., 0] < dv[..., 0]) | ((sv[..., 0] == dv[..., 0]) & (sc[..., 0] < dc[..., 0]))))
        od = dd.copy(); ov = dv.copy(); oc = dc.copy()
        od[..., 0] = np.where(take_b, b, a); ov[..., 0] = np.where(take_b, sv[..., 0], dv[..., 0]); oc[..., 0] = np.where(take_b, sc[..., 0], dc[..., 0])
        if not close.any(): return od, ov, oc, 0
        gd, gv, gc, o = _merge_closest(dd[close], dv[close], dc[close], sd[close], sv[close], sc[close]); od[close] = gd; ov[close] = gv; oc[close] = gc
        return od, ov, oc, o
    d = np.concatenate([dd, sd], -1); v = np.concatenate([dv, sv], -1); c = np.concatenate([dc, sc], -1)
    keep = d <= d.min(-1, keepdims=True) + 1e-9
    d = np.where(keep, d, np.inf); v = np.where(keep, v, np.inf); c = np.where(keep, c, big)
    o = np.lexsort((c, v, d), axis=-1); d = np.take_along_axis(d, o, -1); v = np.take_along_axis(v, o, -1); c = np.take_along_axis(c, o, -1)
    dup = np.zeros(d.shape, bool); dup[..., 1:] = (d[..., 1:] == d[..., :-1]) & np.isfinite(d[..., 1:])
    d = np.where(dup, np.inf, d); v = np.where(dup, np.inf, v); c = np.where(dup, big, c)
    o = np.lexsort((c, v, d), axis=-1); d = np.take_along_axis(d, o, -1)[..., :R]; v = np.take_along_axis(v, o, -1)[..., :R]; c = np.take_along_axis(c, o, -1)[..., :R]
    return d, v, np.where(np.isfinite(d), c, -1), int((np.isfinite(np.take_along_axis(d, o, -1) if False else d).sum(-1) > R).sum())


def closest_winners(cd, cv, cc):
    """the closest match per cell from its band entries: the (vals, column)-smallest entry, -1 where the cell has none"""
    big = np.iinfo(np.int64).max; c = np.where(cc >= 0, cc, big); o = np.lexsort((c, cv), axis=-1)[..., :1]
    w = np.take_along_axis(cc, o, -1)[..., 0]; return w, cd[..., 0]


def _merge_lists(L, rows, cols, mm, mf, keys, nfin):
    """the C smallest (key, column) of the running lists of these rows and the new candidates (a candidate comes from one
    block pair only, so the union has no duplicates); order-independent"""
    C = L["cand"].shape[1]; big = np.iinfo(np.int64).max
    cc = np.concatenate([L["cand"][rows], cols], 1); ck = np.concatenate([L["key"][rows], keys], 1)
    cm = np.concatenate([L["m"][rows], mm], 1); cf = np.concatenate([L["mfl"][rows], mf], 1)
    o = np.lexsort((np.where(cc >= 0, cc, big), ck), axis=1)[:, :C]
    L["cand"][rows] = np.take_along_axis(cc, o, 1); L["key"][rows] = np.take_along_axis(ck, o, 1)
    L["m"][rows] = np.take_along_axis(cm, o, 1); L["mfl"][rows] = np.take_along_axis(cf, o, 1); L["nfin"][rows] += nfin


def _sweep_block(A, B, gidx, cidx, ctx, threads=1, kernel=None, procs=1, slab=None, cl_threads=None):
    """one block pair of the sweep: the counts of its pairs at every locus (or the slab of counts kept from an earlier
    sweep), their medians, keys and identity, then the reductions the sweep keeps. With threads or procs above 1 (the
    single-block path) the per-locus phases run over the loci as pass one did. Returns a dict of per-row results for
    A's rows and, when B differs, B's rows; the slab too when the context asks for it. cl_threads: the threads of the
    closest reduce when they differ from the counts' (the one block pair: the run's whole count)"""
    X, loci, keep, A_, fit, C, rr, ok_r, cdt, want_closest, want_full, key_kind = (ctx[k] for k in ("X", "loci", "keep", "A_", "fit", "C", "rr", "ok_r", "cdt", "closest", "full", "key"))
    nloc = len(loci); bA, bB = len(gidx), len(cidx); P = bA * bB; same = B is A
    kpos = (keep[gidx][:, None] & keep[cidx][None, :]).ravel()               # pairs with both members in the run
    if slab is not None: Ssl, Wsl = slab                                     # kept from an earlier sweep: pairs with an excluded member already 0
    else:
        Ssl = np.empty((nloc, P), cdt); Wsl = np.empty((nloc, P), cdt)
        def counts(li):
            sh, sm = _block_counts(X, loci[li][1], A, B, A_); Wsl[li] = sh.ravel(); Ssl[li] = sm.ravel()
        if procs > 1:
            def chunk_counts(lis):
                out = []
                for li in lis: sh, sm = _block_counts(X, loci[li][1], A, B, A_); out.append((li, sh.ravel().astype(cdt), sm.ravel().astype(cdt)))
                return out
            step = -(-nloc // (4 * procs)); chunks = [list(range(b0, min(b0 + step, nloc))) for b0 in range(0, nloc, step)]
            for res in fork_map(chunk_counts, chunks, procs):
                for li, wrow, srow in res: Wsl[li] = wrow; Ssl[li] = srow
        elif threads > 1:
            with ThreadPoolExecutor(threads) as ex: list(ex.map(counts, range(nloc), chunksize=max(1, nloc // (8 * threads))))
        else:
            for li in range(nloc): counts(li)
        if not kpos.all(): Ssl[:, ~kpos] = 0; Wsl[:, ~kpos] = 0             # a pair with an excluded member has no counts, as before
    G = np.repeat(gidx, bB); Cc = np.tile(cidx, bA); diag = G == Cc                   # a self pair: also when the rows of a second pass meet their own column
    m, m_fl, nl = slab_medians(Ssl, Wsl, diag, rr, ok_r, threads, kernel)
    tot_same = Ssl.sum(0, dtype=np.float64); tot_shared = Wsl.sum(0, dtype=np.float64)
    ident, mj_row, mj_col = _block_keys(m, tot_same, tot_shared, G, Cc, diag)
    if key_kind == "m0":                                                 # sweep 1: the raw divergence of resolved, non-identical pairs
        elig = np.isfinite(m) & (m > m_fl) & ~ident & ~diag; key_row = np.where(elig, m, np.inf); key_col = key_row
    elif key_kind == "m_all":                                            # rows of m on demand: every finite m, identical pairs included
        key_row = np.where(np.isfinite(m) & ~diag, m, np.inf); key_col = key_row
    else: key_row, key_col = mj_row, mj_col
    m2 = m.reshape(bA, bB); f2 = m_fl.reshape(bA, bB); out = dict(overflow=0)
    out["rows"] = (gidx, _topc_rows(key_row.reshape(bA, bB), m2, f2, cidx, C))
    if not same: out["cols"] = (cidx, _topc_rows(key_col.reshape(bA, bB).T, m2.T, f2.T, gidx, C))
    upper = kpos & ~diag & ((G < Cc) if same else np.ones(P, bool))
    out["counters"] = (int((nl[upper] < MIN_LOCI_FOR_M).sum()), int(upper.sum()), int(np.isfinite(m[upper]).sum()))
    ip = ident & (G < Cc) if same else ident; out["ident"] = np.stack([G[ip], Cc[ip]], 1)
    if want_full: out["full"] = (m2, f2)
    if want_closest:
        R = BAND_SLOTS; mr2 = mj_row.reshape(bA, bB); mc2 = mj_col.reshape(bA, bB)
        rd = np.full((bA, nloc, R), np.inf); rv = np.full((bA, nloc, R), np.inf); rc = np.full((bA, nloc, R), -1, np.int64)
        if not same: qd = np.full((bB, nloc, R), np.inf); qv = np.full((bB, nloc, R), np.inf); qc = np.full((bB, nloc, R), -1, np.int64)
        ovf = [0]
        def closest(li):
            ov = Wsl[li].reshape(bA, bB).astype(np.float64); sm = Ssl[li].reshape(bA, bB).astype(np.float64)
            _, ed, ev, ec, o1 = closest_rows(ov, sm, mr2, fit, cidx); rd[:, li] = ed; rv[:, li] = ev; rc[:, li] = ec
            if not same:
                _, ed, ev, ec, o2 = closest_rows(np.ascontiguousarray(ov.T), np.ascontiguousarray(sm.T), np.ascontiguousarray(mc2.T), fit, gidx); qd[:, li] = ed; qv[:, li] = ev; qc[:, li] = ec; o1 += o2
            return o1
        small = bA * bB < CLOSEST_THREADS_MIN                                # measured 31 Oct 2026: Ivan (701 taxa) 2.7 s in threads against 5.2 s forked, Shen (343) 2.6 s forked against 2.8 s in threads
        cl_procs = max(1, min(max(cl_threads or 0, threads, procs), nloc // 16)) if same and small and max(cl_threads or 0, threads, procs) > 1 and _can_fork() else 1
        if cl_procs > 1:                                                     # one block pair of few taxa: the loci over worker processes (its per-locus numpy work is too fine for threads, which queue on the interpreter lock; a fork is a few ms)
            def closest_chunk(lis):
                out = []
                for li in lis:
                    ov = Wsl[li].reshape(bA, bB).astype(np.float64); sm = Ssl[li].reshape(bA, bB).astype(np.float64)
                    _, ed, ev, ec, o1 = closest_rows(ov, sm, mr2, fit, cidx); out.append((li, ed, ev, ec, o1))
                return out
            step = -(-nloc // (4 * cl_procs)); chunks = [list(range(b0, min(b0 + step, nloc))) for b0 in range(0, nloc, step)]
            for res in fork_map(closest_chunk, chunks, cl_procs):
                for li, ed, ev, ec, o1 in res: rd[:, li] = ed; rv[:, li] = ev; rc[:, li] = ec; ovf[0] += o1
        elif (cl_threads or threads) > 1:                                    # no fork here (a worker process, or Windows): the loci in threads
            t_ = min(cl_threads or threads, nloc)
            with ThreadPoolExecutor(t_) as ex: ovf[0] = sum(ex.map(closest, range(nloc), chunksize=max(1, nloc // (8 * t_))))
        else: ovf[0] = sum(closest(li) for li in range(nloc))
        out["closest_rows"] = (gidx, rd, rv, rc)
        if not same: out["closest_cols"] = (cidx, qd, qv, qc)
        out["overflow"] = ovf[0]
    if ctx.get("keep_slab"): out["slab"] = (Ssl, Wsl)
    return out


def sweep_block_bytes(b, nloc, cdt, C, Lmax):
    """bytes a block pair of b taxa holds in flight: the slab of counts, the products, medians and keys (about 200 bytes
    a pair), the one-hot rows of the longest locus, and the lists and band entries it returns"""
    return b * b * (2 * nloc * np.dtype(cdt).itemsize + 200) + b * (2 * C * 32 + 2 * nloc * BAND_SLOTS * 24) + 4 * b * Lmax * 4


def sweep_block_size(n, nloc, cdt, threads, C, Lmax, L=None, A_=4):
    """the taxon block of the sweep: one block up to SWEEP_BLOCK_MIN taxa (its loci run in threads); above it about
    SWEEP_BLOCKS * sqrt(workers) blocks (the workers capped by the sweep's work, SWEEP_FLOPS_PER_WORKER), so every worker
    gets a few block pairs and the last ones are short, never under SWEEP_BLOCK_FLOOR taxa, and small enough that the
    workers' block pairs in flight take at most half the room. L: the alignment's sites (the work estimate)"""
    if n <= SWEEP_BLOCK_MIN: return n
    flops = float(n) * n * float(L if L is not None else Lmax * nloc) * (A_ + 1)       # the products of one sweep: n^2 x sites per letter and for the shared sites
    workers = max(1, min(threads, int(flops / SWEEP_FLOPS_PER_WORKER))); b = max(SWEEP_BLOCK_FLOOR, -(-n // int(math.ceil(SWEEP_BLOCKS * math.sqrt(workers)))))
    room = _mem_room() // 2
    while b > SWEEP_BLOCK_FLOOR and workers * sweep_block_bytes(b, nloc, cdt, C, Lmax) > room: b = max(SWEEP_BLOCK_FLOOR, int(b * 0.9))
    return min(b, n)


def _flatten(v, prefix, out):
    """a block pair's results (a dict of arrays, scalars and tuples of them) as flat name -> array, for one npz file"""
    if isinstance(v, dict):
        for k, x in v.items(): _flatten(x, f"{prefix}{k}.", out)
    elif isinstance(v, tuple):
        for i, x in enumerate(v): _flatten(x, f"{prefix}{i}.", out)
    else: out[prefix[:-1]] = np.asarray(v)


def _unflatten(flat):
    """the inverse: dotted names back to the dict of tuples, 0-d arrays back to scalars"""
    tree = {}
    for key, v in flat.items():
        node = tree; parts = key.split(".")
        for q in parts[:-1]: node = node.setdefault(q, {})
        node[parts[-1]] = v
    def build(node):
        if not isinstance(node, dict): return node.item() if node.ndim == 0 else node
        if node and all(k.isdigit() for k in node): return tuple(build(node[str(i)]) for i in range(len(node)))
        return {k: build(x) for k, x in node.items()}
    return build(tree)


_SWEEP_SEQ = 0
def _results_dir(need):
    """a fresh directory on the memory-backed filesystem for a sweep's results, when it exists and has room for about
    twice `need` bytes; None when the results should come back through the pool's pipe instead"""
    global _SWEEP_SEQ
    try:
        if not os.path.isdir(SWEEP_RESULTS_DIR) or shutil.disk_usage(SWEEP_RESULTS_DIR).free < 2 * need: return None
        _SWEEP_SEQ += 1; d = os.path.join(SWEEP_RESULTS_DIR, f"locusmasker-{os.getpid()}-{_SWEEP_SEQ}"); os.makedirs(d); return d
    except OSError: return None


def _sweep_worker(pair):
    """a block pair in a worker process: the context was set in the parent before the fork"""
    ctx = _SWEEP; blocks = ctx["blocks"]; i, j = pair; A, B = blocks[i], (blocks[i] if i == j else blocks[j])
    return _sweep_block(A, B if i != j else A, ctx["gidx"][i], ctx["gidx"][j], ctx, 1, ctx["kernel_serial"])


def sweep_pairs(X, loci, keep, A_, fit, C, b, rr, ok_r, threads, want_closest, want_full, key, cdt, progress=None, stage="pairs", kernel=None, kernel_serial=None, rows_only=None, slab=None, keep_slab=False):
    """the sweep over block pairs (rows_only: the given taxa against everyone, for the second passes): returns the
    lists (cand, key, m, mfl, nfin), the closest-match band entries, the identical pairs, the counters and the full m.
    slab: the counts of every pair kept from the one block pair's first sweep (nothing is recomputed from the residues);
    keep_slab: return that slab from a one-block sweep"""
    global _SWEEP
    n, nloc = X.shape[0], len(loci); kept = None
    def part(ga, gb):
        """the kept counts of these rows against these columns, or None when counts are computed"""
        if slab is None: return None
        pos = (np.asarray(ga, np.int64)[:, None] * n + np.asarray(gb, np.int64)[None, :]).ravel(); return (slab[0][:, pos], slab[1][:, pos])
    if rows_only is not None:
        rows = np.asarray(rows_only, np.int64); blocks_r = [rows[i:i + b] for i in range(0, len(rows), b)]; blocks_c = [slice(i, min(i + b, n)) for i in range(0, n, b)]
        items = [(i, j) for i in range(len(blocks_r)) for j in range(len(blocks_c))]
        gidx_r = [np.asarray(r) for r in blocks_r]; gidx_c = [np.arange(c.start, c.stop) for c in blocks_c]
        nrows = len(rows)
    else:
        blocks_c = [slice(i, min(i + b, n)) for i in range(0, n, b)]; gidx_c = [np.arange(c.start, c.stop) for c in blocks_c]
        items = [(i, j) for i in range(len(blocks_c)) for j in range(i, len(blocks_c))]; nrows = n
    Cmax = min(C, n)
    L = dict(cand=np.full((nrows, Cmax), -1, np.int64), key=np.full((nrows, Cmax), np.inf), m=np.full((nrows, Cmax), np.nan), mfl=np.full((nrows, Cmax), np.nan), nfin=np.zeros(nrows, np.int64))
    R = BAND_SLOTS
    cd = np.full((nrows, nloc, R), np.inf) if want_closest else None; cv = np.full((nrows, nloc, R), np.inf) if want_closest else None; cc = np.full((nrows, nloc, R), -1, np.int64) if want_closest else None
    full = (np.full((n, n), np.nan), np.full((n, n), np.nan)) if want_full else None
    ident = []; counters = [0, 0, 0]; overflow = 0
    pos = {int(g): k for k, g in enumerate(rows)} if rows_only is not None else None     # row taxon -> row of the results
    def rowpos(g): return np.array([pos[int(x)] for x in g]) if pos is not None else g
    def absorb(out):
        nonlocal overflow
        for side in ("rows", "cols"):
            if side in out:
                g, (cols, mm, mf, keys, nfin) = out[side]
                if rows_only is not None and side == "cols": continue      # the second pass keeps the rows' results only
                _merge_lists(L, rowpos(g), cols, mm, mf, keys, nfin)
        for side in ("closest_rows", "closest_cols"):
            if side in out:
                g, rd, rv, rc = out[side]
                if rows_only is not None and side == "closest_cols": continue
                rp = rowpos(g); cd[rp], cv[rp], cc[rp], o = _merge_closest(cd[rp], cv[rp], cc[rp], rd, rv, rc); overflow += o
        if "full" in out and full is not None:
            m2, f2 = out["full"]; g, c = out["gidx"], out["cidx"]
            full[0][np.ix_(g, c)] = m2; full[1][np.ix_(g, c)] = f2; full[0][np.ix_(c, g)] = m2.T; full[1][np.ix_(c, g)] = f2.T
        ident.append(out["ident"]); overflow += out["overflow"]
        for k in range(3): counters[k] += out["counters"][k]
    ctx = dict(X=X, loci=loci, keep=keep, A_=A_, fit=fit, C=Cmax, rr=rr, ok_r=ok_r, cdt=cdt, closest=want_closest, full=want_full, key=key,
               blocks=blocks_c if rows_only is None else None, gidx=gidx_c if rows_only is None else None, kernel_serial=kernel_serial, keep_slab=keep_slab)
    def block_pair(i, j):
        if rows_only is not None: A, B = blocks_r[i], blocks_c[j]; ga, gb = gidx_r[i], gidx_c[j]
        else: A, B = blocks_c[i], blocks_c[j]; ga, gb = gidx_c[i], gidx_c[j]
        return (A, (A if (rows_only is None and i == j) else B), ga, gb)
    bmax = max(len(g) for g in gidx_c); rdir = None
    if len(items) > 1 and threads > 1 and _can_fork():               # the worker path: where its results go
        ret = 2 * bmax * (Cmax * 32 + nloc * BAND_SLOTS * 24) + (2 * bmax * bmax * 8 if want_full else 0)     # bytes one block pair returns
        rdir = _results_dir(min(threads, len(items)) * ret)
    if rows_only is None:
        inflight = sweep_block_bytes(bmax, nloc, cdt, Cmax, max(len(st) for _, st in loci))
        LOG(f"pair sweep {stage}: {len(items)} block pair{'s' if len(items) > 1 else ''} of {bmax:,} taxa, medians by {'the compiled kernel (numba)' if kernel else 'numpy chunks'}, "
            f"about {fmem(inflight)} per block pair in flight" + (f", {min(threads, len(items))} worker processes" if len(items) > 1 and threads > 1 else "")
            + (f", results through {SWEEP_RESULTS_DIR}" if rdir else ""))
    if progress: progress.start(stage, len(items))
    if len(items) == 1 and rows_only is None:                            # one block pair: the per-locus phases over the loci, as pass one did
        A, B, ga, gb = block_pair(*items[0])
        per_locus = 3 * n * n * 4 + 2 * n * max(len(st) for _, st in loci) * 4
        sample = [loci[i][1] for i in np.linspace(0, nloc - 1, min(nloc, 64)).astype(int)]
        kbar = float(np.mean([np.mean([(X[:, st] == i).any(0).mean() for i in range(A_)]) * A_ for st in sample]))
        flops = 2.0 * n * n * float(np.mean([len(st) for _, st in loci])) * (1.0 + kbar)
        workers = max(1, min(threads, nloc, int((_mem_room() // 2) // per_locus), int(round(flops / PASS1_FLOPS_PER_THREAD))))
        procs = min(threads, max(1, nloc // PASS1_LOCI_PER_WORKER)) if workers <= 2 and threads > 1 and 2 * nloc * n * n * np.dtype(cdt).itemsize <= PASS1_RETURN_MAX else 1
        out = _sweep_block(A, B, ga, gb, ctx, threads=max(workers, 1) if procs == 1 else 1, kernel=kernel, procs=procs, slab=slab, cl_threads=threads)
        out["gidx"], out["cidx"] = ga, gb; absorb(out); kept = out.get("slab")
        if progress: progress.tick()
    else:
        workers = min(threads, len(items)) if threads > 1 and _can_fork() else 1
        per_worker = max(1, threads // max(workers, 1))                        # block pairs fewer than the threads: each runs its loci in threads
        order = sorted(items, key=lambda ij: (ij[0] == ij[1], ij))            # diagonal blocks last: about half the work of the others (symmetric products, one closest direction), so the tail is short
        if workers > 1:
            _SWEEP = ctx
            try:
                def run(ij):
                    A, B, ga, gb = block_pair(*ij); out = _sweep_block(A, B, ga, gb, ctx, per_worker, kernel_serial, slab=part(ga, gb)); out["gidx"], out["cidx"] = ga, gb
                    if rdir is None: return out
                    flat = {}; _flatten(out, "", flat); path = os.path.join(rdir, f"{ij[0]}-{ij[1]}.npz"); np.savez(path, **flat); return path
                for r in fork_map(run, order, workers):
                    if rdir is not None:
                        with np.load(r) as z: out = _unflatten({k: z[k] for k in z.files})
                        os.remove(r)
                    else: out = r
                    absorb(out); progress and progress.tick()
            finally:
                _SWEEP = None
                if rdir is not None: shutil.rmtree(rdir, ignore_errors=True)
        else:
            for ij in order:
                A, B, ga, gb = block_pair(*ij); out = _sweep_block(A, B, ga, gb, ctx, threads, kernel, slab=part(ga, gb)); out["gidx"], out["cidx"] = ga, gb; absorb(out)
                if progress: progress.tick()
    if progress: progress.finish()
    ident = np.concatenate(ident) if ident else np.zeros((0, 2), np.int64)
    return dict(lists=L, closest=(cd, cv, cc), ident=ident, counters=tuple(counters), full=full, overflow=overflow, blocks=len(blocks_c), b=b, items=len(items), slab=kept)


def classes_from_pairs(n, pairs):
    """each taxon's class of identical genomes (its lowest index), the fixed point the old engine iterated on the
    full relation, run on the list of identical pairs"""
    cls = np.arange(n)
    if len(pairs) == 0: return cls
    g, c = pairs[:, 0], pairs[:, 1]
    while True:
        new = cls.copy(); np.minimum.at(new, g, cls[c]); np.minimum.at(new, c, cls[g])
        if np.array_equal(new, cls): return cls
        cls = new


def list_lookup(L, gi, ci):
    """the position of each pair's partner in the focal taxon's list, -1 where it is not there"""
    gi = np.asarray(gi, np.int64); ci = np.asarray(ci, np.int64); cand = L["cand"][gi]
    hit = cand == ci[:, None]; pos = hit.argmax(1); return np.where(hit.any(1), pos, -1)


def _restart_header(a, n, nloc, L, dtype, Kn, C, thin, full):
    """what a restart point must have been made from: the input (name, bytes, sha256), the partitions, and every option
    and constant the sweeps' results depend on"""
    if getattr(a, "input_sha256", None) is None: a.input_sha256 = _sha256(a.aln)
    return dict(tool="LocusMasker", version=__version__, input=dict(file=os.path.basename(a.aln), bytes=os.path.getsize(a.aln), sha256=a.input_sha256),
                partitions_sha256=_sha256(a.partitions) if getattr(a, "partitions", None) else None, datatype=dtype, fit_overlap=a.fit_overlap, K_eff=int(Kn),
                candidates=int(C), full_m_max_taxa=FULL_M_MAX_TAXA, full=bool(full), taxa=int(n), loci=int(nloc), sites=int(L), thin=[int(t) for t in thin])


def write_restart(a, point, header, arrays):
    """a restart point: the arrays and then the header, into a .tmp directory that is renamed into place, so the point
    is whole or absent. A write error warns and the run goes on"""
    d = os.path.join(a.out, "restart", point); tmp = d + ".tmp"
    try:
        shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
        for k, v in arrays.items(): np.save(os.path.join(tmp, k + ".npy"), np.ascontiguousarray(v))
        with open(os.path.join(tmp, "header.json"), "w", encoding="utf-8") as fh: json.dump(dict(header, point=point, written=time.strftime("%Y-%m-%d %H:%M:%S")), fh, indent=1)
        shutil.rmtree(d, ignore_errors=True); os.rename(tmp, d)
        LOG(f"restart point written: {d} ({'after the locus rates' if point == 'rates' else 'after the pair sweeps'}): a rerun of this input with these options in this --out starts from it")
    except OSError as e:
        LOG(f"WARNING: the restart point {d} could not be written ({e.strerror or e}); the run goes on without it")
    if STOP_AFTER == point: raise SystemExit(f"stopped after the {point} restart point (LOCUSMASKER_STOP_AFTER)")


def read_restart(a, point, header):
    """the arrays of a matching restart point and when it was written, or None: absent, or incomplete (said so and
    removed). One from another input or other options stops the run naming what differs; --restart ignores them all"""
    if getattr(a, "restart", False): return None
    d = os.path.join(a.out, "restart", point)
    if not os.path.isdir(d): return None
    try:
        with open(os.path.join(d, "header.json"), encoding="utf-8") as fh: got = json.load(fh)
        names_ = RESTART_ARRAYS[point] + (("m_full", "mfl_full") if point == "pairs" and got.get("full") else ())
        arrays = {k: np.load(os.path.join(d, k + ".npy")) for k in names_}
    except (OSError, ValueError, KeyError) as e:
        LOG(f"note: the restart point {d} is incomplete ({type(e).__name__}: {e}): removed; the sweeps run"); shutil.rmtree(d, ignore_errors=True); return None
    diff = [k for k in header if got.get(k) != header[k]]
    if diff: sys.exit(f"the restart point {d} is from a different run ({', '.join(diff)} differ): --restart ignores it, or use another --out")
    LOG(f"restarted after {'the locus rates' if point == 'rates' else 'the pair sweeps'} from the restart point written {got['written']}: "
        f"{'sweep 1 and the rates are' if point == 'rates' else 'the sweeps are'} not repeated")
    return arrays, got["written"]


def prepare(a, names, rows, loci, dtype, K_eff, t0, full=None):
    """everything the test needs before it runs: the alignment as residue codes, the locus plateaus and rates, and the
    pair structures of the sweep (each taxon's candidate list with m, its floor and the ranking key; the closest match
    at every locus; the classes of identical genomes; the full m up to FULL_M_MAX_TAXA), and the counts provider: the
    kept slab of every pair's counts when the sweep was one block pair, else the residues. a needs fit_overlap and
    threads; a.candidates and a.block_taxa are read when present"""
    n, L = len(names), len(rows[0]); A_ = 4 if dtype == "dna" else 20; nloc = len(loci)
    X = encode_alignment(names, rows, dtype)
    focal_sites = np.stack([(X[:, st] >= 0).sum(1) for _, st in loci], 1).astype(np.uint32)     # each taxon's residues at each locus
    present = (focal_sites >= a.fit_overlap).sum(1); thin = np.flatnonzero(present < MIN_LOCI_PER_TAXON)
    keep = np.ones(n, bool); keep[thin] = False
    if len(thin):
        LOG(f"{len(thin)} taxa with fewer than {MIN_LOCI_PER_TAXON} loci present are excluded from the test and left unmasked: "
            + ", ".join(names[i] for i in thin[:20]) + (f" and {len(thin) - 20} more (the report lists them)" if len(thin) > 20 else ""))
    plateaus = locus_plateaus(X, loci, A_, min(a.threads, 8))
    longest = max(len(st) for _, st in loci); assert longest < (1 << 24), "a locus of 2^24 sites or more: the float32 products would not be exact"
    cdt = np.uint16 if longest < 65536 else np.uint32
    C = max(1, int(getattr(a, "candidates", None) or CANDIDATES)); b = int(getattr(a, "block_taxa", None) or sweep_block_size(n, nloc, cdt, a.threads, C, longest, L, A_)); b = max(1, min(b, n))
    if full is None: full = n <= FULL_M_MAX_TAXA
    kern = _pair_kernel(a.threads) if nloc * n * n >= NUMBA_PAIR_LOCI_MIN else None
    kernel = kern[0] if kern else None; kernel_serial = _pair_medians_serial if kern else None
    if kern:                                                             # compile both kernels here, before any worker is forked
        tiny_s = np.zeros((nloc, 1), cdt); tiny_w = np.ones((nloc, 1), cdt); dg = np.zeros(1, bool)
        kernel(tiny_s, tiny_w, 0, 1, np.ones(nloc), np.ones(nloc, bool), dg, MIN_LOCI_FOR_M, 100.0 * HALF, PAIR_BLOCK)
        kernel_serial(tiny_s, tiny_w, 0, 1, np.ones(nloc), np.ones(nloc, bool), dg, MIN_LOCI_FOR_M, 100.0 * HALF, PAIR_BLOCK)
        if _closest_rows_nb is not None: _closest_rows_nb(np.ones((1, 1)), np.ones((1, 1)), np.ones((1, 1)), 1, np.zeros(1, np.int64), BAND_SLOTS, HALF)
    LOG(f"pair sweep: blocks of {b:,} taxa, {C} candidates per taxon, the full m matrix {'kept' if full else 'not kept'}")
    common = dict(X=X, loci=loci, keep=keep, A_=A_, fit=a.fit_overlap, C=C, b=b, threads=a.threads, cdt=cdt, kernel=kernel, kernel_serial=kernel_serial)
    Kn = min(K_eff, n - 1); restart = None; n_starved = 0
    points = getattr(a, "out", None) is not None and getattr(a, "aln", None) is not None       # restart points need an output folder and the input file
    header = _restart_header(a, n, nloc, L, dtype, Kn, C, thin, full) if points else None
    pt2 = read_restart(a, "pairs", header) if points else None
    pt1 = read_restart(a, "rates", header) if points and pt2 is None else None
    if pt2 is not None:                                                  # both sweeps done by an earlier run: nothing is recomputed
        R_, when = pt2; restart = dict(after="the pair sweeps", written=when); ident = R_["ident"]; rate = R_["rate"]; ok_r = R_["ok_r"]; rr = np.where(ok_r, rate, np.nan)
        S2 = dict(lists={k: R_[k] for k in ("cand", "key", "m", "mfl", "nfin")}, closest=(R_["cd"], R_["cv"], R_["cc"]), counters=tuple(int(x) for x in R_["counters"]), overflow=0,
                  full=(R_["m_full"], R_["mfl_full"]) if full else None)
        counts = ResidueCounts(X, loci, keep, A_)
    elif pt1 is not None:                                                # sweep 1 and the rates done by an earlier run
        R_, when = pt1; restart = dict(after="the locus rates", written=when); ident = R_["ident"]; rate = R_["rate"]; ok_r = R_["ok_r"]; rr = np.where(ok_r, rate, np.nan)
        counts = ResidueCounts(X, loci, keep, A_)
    else:
        # ---- sweep 1: the raw divergences, for the rate-vote neighbours
        S1 = sweep_pairs(rr=None, ok_r=None, want_closest=False, want_full=False, key="m0", progress=PROG, stage="pairs0", keep_slab=b >= n, **common)
        if S1["slab"] is not None:                                       # one block pair: its counts of every pair serve the whole run, as the old engine's did
            common["slab"] = S1["slab"]; LOG(f"the pair counts of the one block pair ({fmem(2 * S1['slab'][0].nbytes)}) are kept for the run")
        counts = SlabCounts(*common["slab"], n) if "slab" in common else ResidueCounts(X, loci, keep, A_)
        ident = S1["ident"]
    cls = classes_from_pairs(n, ident); has_dup = len(ident) > 0
    if has_dup: LOG(f"{int((cls != np.arange(n)).sum())} genomes are identical to another at every shared site: an identical genome is never a comparator, and no two comparators are identical to each other")
    if restart is None: LOG(f"per-locus identities done, {time.time() - t0:.0f}s")
    rep = cls == np.arange(n)
    def neighbours(L1, rows_of):
        gi, ci, starved = [], [], []
        for k, g in enumerate(rows_of):
            cands = L1["cand"][k]; cands = cands[cands >= 0]; picked = cands[rep[cands]][:Kn]
            gi += [g] * len(picked); ci += picked.tolist()
            if len(picked) < Kn and L1["nfin"][k] > L1["cand"].shape[1]: starved.append(g)
        return np.array(gi, np.int64), np.array(ci, np.int64), starved
    if restart is None:                                                  # the locus rates, voted by each taxon's nearest resolved non-identical neighbours
        gi0, ci0, starved = neighbours(S1["lists"], range(n)); n_starved = len(starved)
        if starved:                                                          # the list ran out before Kn representatives: their full rows
            LOG(f"rate neighbours: {len(starved)} taxa ran out of candidates (copies of one genome fill the list): full rows computed")
            E = sweep_pairs(rr=None, ok_r=None, want_closest=False, want_full=False, key="m0", rows_only=np.array(starved), **{**common, "C": n})
            gi_e, ci_e, _ = neighbours(E["lists"], starved)
            keepmask = ~np.isin(gi0, starved); gi0 = np.concatenate([gi0[keepmask], gi_e]); ci0 = np.concatenate([ci0[keepmask], ci_e])
            lk = list_lookup(E["lists"], np.array([starved.index(g) for g in gi_e]), ci_e); m0n_e = np.maximum(E["lists"]["m"][np.array([starved.index(g) for g in gi_e]), lk], E["lists"]["mfl"][np.array([starved.index(g) for g in gi_e]), lk])
        lk = list_lookup(S1["lists"], gi0[:len(gi0) - (len(gi_e) if starved else 0)], ci0[:len(ci0) - (len(ci_e) if starved else 0)])
        base_rows = gi0[:len(gi0) - (len(gi_e) if starved else 0)]
        m0n = np.maximum(S1["lists"]["m"][base_rows, lk], S1["lists"]["mfl"][base_rows, lk])
        if starved: m0n = np.concatenate([m0n, m0n_e])
        del S1
        rate = np.full(nloc, np.nan)
        def locus_rate(li):
            o_i, s_i = counts.at(li, gi0, ci0); o_ = o_i.astype(np.float32); s_ = s_i.astype(np.float32)
            with np.errstate(invalid="ignore", divide="ignore"):
                d_ = np.where(o_ > 0, 100.0 * (1.0 - s_ / np.maximum(o_, 1)), np.nan)
            val = np.isfinite(d_) & (o_ >= a.fit_overlap)
            if val.sum() < MIN_FIT_PAIRS: return
            P_ = plateaus[li]
            corr = lambda p: -P_ * np.log(1.0 - np.minimum(np.asarray(p, np.float64), 0.98 * P_) / P_)
            rate[li] = float(np.median(corr(d_[val] + 100.0 * HALF / o_[val]) / corr(m0n[val])))
        PROG.start("rates", nloc)
        if a.threads > 1 and nloc > 1 and "slab" not in common:             # the counts from the residues release the interpreter lock; gathers from the kept slab are quicker serial
            with ThreadPoolExecutor(min(a.threads, nloc)) as ex:
                for _ in ex.map(locus_rate, range(nloc), chunksize=max(1, nloc // (8 * a.threads))): PROG.tick()
        else:
            for li in range(nloc): locus_rate(li); PROG.tick()
        PROG.finish()
        lsites = np.array([len(s) for _, s in loci], dtype=float)
        ok_r = np.isfinite(rate) & (rate > 0)
        if ok_r.sum() >= 3:
            scale = wmedian_stack(rate[ok_r, None], lsites[ok_r, None])[0]
            rate = rate / max(scale, 1e-6)
        rr = np.where(ok_r, rate, np.nan)
        if points: write_restart(a, "rates", header, dict(ident=ident, rate=rate, ok_r=ok_r))
    if pt2 is None:
        # ---- sweep 2: the rate-normalised divergences, the lists, the closest match at every locus
        S2 = sweep_pairs(rr=rr, ok_r=ok_r, want_closest=True, want_full=full, key="mj", progress=PROG, stage="pairs", **common)
        if points: write_restart(a, "pairs", header, dict(ident=ident, rate=rate, ok_r=ok_r, **S2["lists"], cd=S2["closest"][0], cv=S2["closest"][1], cc=S2["closest"][2],
                                                           counters=np.array(S2["counters"], np.int64), **({"m_full": S2["full"][0], "mfl_full": S2["full"][1]} if full else {})))
    nl_low, npk, fin_pairs = S2["counters"]; fb = nl_low / npk if npk else 0.0
    LOG(f"pair divergences done, {time.time() - t0:.0f}s"
          + (f"; note: no locus has a rate (no nearest-neighbour pairs resolved beyond their resolution floor), so no pair has a divergence" if ok_r.sum() == 0 else
             f"; note: {fb:.0%} of taxon pairs overlap in <{MIN_LOCI_FOR_M} loci with a rate and are never compared" if fb > 0.05 else ""))
    if S2["overflow"]: LOG(f"note: {S2['overflow']:,} cells had more distinct closest-match rates within 1e-9 of the minimum than {BAND_SLOTS} slots")
    m_full, mfl_full = S2["full"] if full else (None, None)
    return argparse.Namespace(X=X, keep=keep, thin=thin, present=present, focal_sites=focal_sites, plateaus=plateaus, rate=rate, rr=rr, ok_r=ok_r, cls=cls, has_dup=has_dup,
                              lists=S2["lists"], closest=S2["closest"], m_full=m_full, mfl_full=mfl_full, fin_pairs=fin_pairs, C=C, b=b, cdt=cdt, A_=A_, fit=a.fit_overlap,
                              sweep_args=common, ext={}, counts=counts, restart=restart, n_starved=n_starved)


def pair_m_lists(P, gi, ci):
    """m and its floor for the pairs (gi[k], ci[k]) from the focal taxa's lists (a comparator is always in its focal's
    list), the extended rows of the second pass, or the full matrix when it is kept"""
    gi = np.asarray(gi, np.int64); ci = np.asarray(ci, np.int64); m_ = np.full(len(gi), np.nan); f_ = np.full(len(gi), np.nan)
    if len(gi) == 0: return m_, f_
    if P.m_full is not None: return P.m_full[gi, ci], P.mfl_full[gi, ci]       # the lists were filled from the same block results: the same numbers
    pos = list_lookup(P.lists, gi, ci); hit = pos >= 0
    m_[hit] = P.lists["m"][gi[hit], pos[hit]]; f_[hit] = P.lists["mfl"][gi[hit], pos[hit]]
    if P.ext:
        for k in np.flatnonzero(~hit & np.isin(gi, np.array(list(P.ext)))).tolist():
            cand, mm, mf, _ = P.ext[int(gi[k])]; j = np.flatnonzero(cand == ci[k])
            if len(j): m_[k] = mm[j[0]]; f_[k] = mf[j[0]]; hit[k] = True
    if P.m_full is not None and not hit.all():
        miss = ~hit; m_[miss] = P.m_full[gi[miss], ci[miss]]; f_[miss] = P.mfl_full[gi[miss], ci[miss]]
    return m_, f_


def pair_m(P, gi, ci, step=2048):
    """m and its floor for arbitrary pairs (gi[k], ci[k]): from the full matrix when it is kept, else by the sweep's own
    counts and medians on the pairs' rows (the same bits; NaN for a taxon with itself, as in m_rows). Serial, with the
    serial kernel: it runs inside the test stage's worker processes"""
    gi = np.asarray(gi, np.int64); ci = np.asarray(ci, np.int64)
    if P.m_full is not None: return P.m_full[gi, ci], P.mfl_full[gi, ci]
    m_ = np.full(len(gi), np.nan); f_ = np.full(len(gi), np.nan); loci = P.sweep_args["loci"]; kern = P.sweep_args["kernel_serial"]
    for k0 in range(0, len(gi), step):
        g, c = gi[k0:k0 + step], ci[k0:k0 + step]; U, inv = np.unique(np.r_[g, c], return_inverse=True); XU = P.X[U]; kU = P.keep[U]
        Wsl = np.empty((len(loci), len(g)), P.cdt); Ssl = np.empty((len(loci), len(g)), P.cdt)
        for li, (_, st) in enumerate(loci):
            o, s = pair_counts_at(XU, st, inv[:len(g)], inv[len(g):], kU); Wsl[li] = o; Ssl[li] = s
        diag = g == c; mm, mf, _ = slab_medians(Ssl, Wsl, diag, P.rr, P.ok_r, 1, kern)
        mm[diag] = np.nan; mf[diag] = np.nan; m_[k0:k0 + step] = mm; f_[k0:k0 + step] = mf
    return m_, f_


def m_rows(P, rows):
    """full rows of m for a few taxa: from the matrix when it is kept, else computed on demand by the sweep (their pairs
    against everyone, every finite m, identical pairs included)"""
    rows = np.asarray(rows, np.int64)
    if P.m_full is not None: return P.m_full[rows]
    n = P.X.shape[0]; E = sweep_pairs(rr=P.rr, ok_r=P.ok_r, want_closest=False, want_full=False, key="m_all", rows_only=rows, **{**P.sweep_args, "C": n})
    out = np.full((len(rows), n), np.nan)
    for k in range(len(rows)):
        c = E["lists"]["cand"][k]; fin = c >= 0; out[k, c[fin]] = E["lists"]["m"][k][fin]
    return out


def bs_keys_taxa(n):
    """each taxon's hash key for its pool draw"""
    return _sm64(np.arange(1, n + 1, dtype=np.uint64) ^ np.uint64(BS_SEED))


def bs_keys_cells(gidx, li):
    """each cell's hash key for its comparator draw at locus li"""
    with np.errstate(over="ignore"):                                   # uint64 wraparound is the point of the mix
        mix = np.uint64(li + 1) * np.uint64(0x9E3779B97F4A7C15)
    return _sm64(_sm64(np.asarray(gidx, np.int64).astype(np.uint64) + np.uint64(1)) ^ mix ^ np.uint64(BS_SEED))


def bs_prepare(P, K_eff, threads, progress=None):
    """the broad-sample pools: per taxon min(BS_POOL x K_eff, n-1) candidates drawn by hash (self excluded), and
    the genome-wide m* of every (taxon, pool member) pair, from the kept matrix or computed per pair above
    FULL_M_MAX_TAXA. Everything downstream reads mv by pool slot, so no other pair's m is ever needed"""
    n = P.X.shape[0]
    k = min(BS_PER_COMP * K_eff, n - 1); Pn = min(BS_POOL * K_eff, n - 1)
    if k < BS_MIN: return None
    if Pn >= n - 1:
        pool = np.tile(np.arange(n - 1, dtype=np.int64), (n, 1))              # every other taxon (few taxa: n - 1 <= the pool)
    else:
        pool = _hash_draw(bs_keys_taxa(n), n - 1, Pn)
    pool += (pool >= np.arange(n)[:, None]).astype(np.int64)                  # values in [0, n) skipping self
    if P.m_full is not None:
        g = np.repeat(np.arange(n), Pn); c = pool.ravel()
        m_ = P.m_full[g, c]; f_ = P.mfl_full[g, c]
        ms = np.where(np.isfinite(f_), np.fmax(m_, f_), m_).reshape(n, Pn)
    else:                                                                     # per-pair m for the pool pairs, forked by taxon slices
        if progress: progress.start("broad-sample pool m", n)
        def rows_ms(sl):
            out = np.empty((len(sl), Pn))
            for i, g_ in enumerate(sl):
                m_, f_ = pair_m(P, np.full(Pn, g_, np.int64), pool[g_])
                out[i] = np.where(np.isfinite(f_), np.fmax(m_, f_), m_)
            return sl[0], out
        workers = max(1, min(threads, 96)) if threads > 1 else 1
        slices = [np.arange(i, min(i + 256, n)) for i in range(0, n, 256)]
        ms = np.empty((n, Pn))
        for s0, out in fork_map(rows_ms, slices, min(workers, len(slices))):
            ms[s0:s0 + len(out)] = out
            if progress: progress.tick(len(out))
        if progress: progress.finish()
    return argparse.Namespace(k=k, Pn=Pn, pool=pool.astype(np.int32), mstar=ms)


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
    if n > FULL_M_MAX_TAXA: raise ValueError(f"{n} taxa: the full matrix is kept up to {FULL_M_MAX_TAXA}")
    P = prepare(argparse.Namespace(fit_overlap=fit_overlap, threads=max(1, min(threads, _cpus()))), names, rows, loci, datatype, K_eff, time.time(), full=True)
    m_hi = np.maximum(P.m_full, P.mfl_full); np.fill_diagonal(m_hi, np.nan)
    return names, P.m_full, m_hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True); ap.add_argument("--partitions"); ap.add_argument("--out", required=True)
    ap.add_argument("--datatype", required=True, choices=("dna", "aa"))
    ap.add_argument("--n-comp", type=int, default=5); ap.add_argument("--fit-overlap", type=int, default=20)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--donor-threshold", type=float, default=2.0, help="the direction of a too-similar pair is called (the copy masked, its source not) when the larger of the "
                    "two members' test-1 excesses, in absolute value, is more than this many times the smaller; otherwise the direction is uncertain and both members are masked (default 2)")
    ap.add_argument("--mask-similar-only", action="store_true"); ap.add_argument("--mask-different-only", action="store_true")
    ap.add_argument("--no-broad-sample", action="store_true", help="skip the broad-sample test (the third test: a broad sample of taxa against each sequence, catching transfers "
                    "and contamination that moved a whole clade, which the closest-relative tests cannot see); this is for comparison runs")
    ap.add_argument("--threads", type=int, default=1, help="at most this many threads or worker processes (0: every CPU this process may use); each stage takes as many as pay off")
    ap.add_argument("--candidates", type=int, help=f"the closest genomes each taxon keeps as its candidate comparators (default {CANDIDATES}); a taxon whose list runs out at a locus gets its full row computed afterwards, so any value gives the same results")
    ap.add_argument("--block-taxa", type=int, help="the taxon block of the pair sweep (default from the taxa count, the threads and the memory room); any block size gives the same results, so this is for testing")
    ap.add_argument("--no-numba", action="store_true", help="the pair stage in numpy even when numba is installed (the compiled kernel gives the same numbers, 3-8x faster at scale; this is for checking that)")
    ap.add_argument("--memory", type=float, help="GB of memory the run may use (default: half of what is available when it starts; 0: no budget); the stages size their threads and chunks to stay within it, and a run that cannot fit stops before it starts")
    ap.add_argument("--progress", choices=("auto", "bar", "none"), default="auto", help="auto: a progress bar on a terminal, progress lines at 2%% steps otherwise")
    ap.add_argument("--report-codes", choices=("fails", "masked", "none"), default="fails", help="site codes in the report: for every cell with a failing comparator (default), masked cells only, or none")
    ap.add_argument("--report-sequences", choices=("auto", "all", "fails", "masked", "none"), default="auto",
                    help="sequences in the report: the whole alignment, the rows of every cell with a failing comparator, of masked cells only, or none; auto: all when the alignment is small, else fails")
    ap.add_argument("--no-report", action="store_true", help="write no report.html / report.npz")
    ap.add_argument("--restart", action="store_true", help="ignore the restart points an earlier run left in --out (restart/: the pair sweeps' results, reused by a rerun of the same input and options) and recompute everything")
    ap.add_argument("--all-tests", action="store_true", help="also write tests.tsv.gz: one row per test the run performed (every comparison, every pair of taxa at every locus, every direction test) with its statistic, exact p-value and result; large")
    ap.add_argument("--name", help="dataset name shown in the report (default: the alignment file's stem, or its folder when the file is called alignment.*)")
    a = ap.parse_args(); K = a.n_comp
    global PROG, STOP_AFTER
    if os.environ.get("LOCUSMASKER_STOP_AFTER") in ("rates", "pairs"): STOP_AFTER = os.environ["LOCUSMASKER_STOP_AFTER"]
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    stale = [f for f in ("tests.tsv.gz",) if not a.all_tests and os.path.exists(os.path.join(a.out, f))] + \
            [f for f in ("report.html", "report.npz") if a.no_report and os.path.exists(os.path.join(a.out, f))] + \
            [f for f in ("cells.tsv", "pairs.npy", "similar.tsv", "per_locus.tsv", "report_data.js") if os.path.exists(os.path.join(a.out, f))]
    for f in ("taxa.tsv", "loci.tsv", "summary.json", "report.html.part", "report.npz.part"):        # the tables of a run above the one-file ceiling, and torn files: never stale
        if os.path.exists(os.path.join(a.out, f)): os.remove(os.path.join(a.out, f))
    try: LOG.open(os.path.join(a.out, "run.log"))
    except OSError as e: sys.exit(f"cannot write in {a.out}: {e.strerror or e}")
    PROG = Progress(a.progress, t0); LOG.bar = PROG.mode == "bar"
    LOG(f"LocusMasker {__version__} ({_git_commit()}) on {time.strftime('%Y-%m-%d %H:%M:%S')}; python {platform.python_version()}, numpy {np.__version__}, {platform.platform()}")
    LOG("command: " + " ".join(sys.argv))
    if stale: LOG(f"note: {a.out} holds files this run does not write ({', '.join(stale)}): they are from an earlier run")
    ncpu = _cpus(); asked = a.threads; a.threads = ncpu if asked <= 0 else max(1, min(asked, ncpu))
    if asked > ncpu: LOG(f"threads: {asked} asked, {ncpu} CPUs available to this process: using {a.threads}")
    global MEM_BUDGET, USE_NUMBA
    USE_NUMBA = not a.no_numba
    avail = _mem_available(); MEM_BUDGET = None if a.memory is not None and a.memory <= 0 else int(a.memory * (1 << 30)) if a.memory else avail // 2
    if MEM_BUDGET is not None and MEM_BUDGET > avail: LOG(f"WARNING: --memory {a.memory:g} GB is more than the {avail / (1 << 30):.1f} GB available now")
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    if a.candidates is not None and a.candidates < 1: sys.exit("--candidates must be at least 1")
    if a.block_taxa is not None and a.block_taxa < 1: sys.exit("--block-taxa must be at least 1")
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
    floor = report_floor_bytes(n, len(loci), K_eff) // 3
    if not a.no_report and floor > REPORT_ONE_FILE_MAX:
        LOG(f"report.html is likely to exceed the one-file ceiling ({fmem(REPORT_ONE_FILE_MAX)} of bundle: this run's per-cell members alone come to about {fmem(floor)} "
            f"compressed): report.npz and the tables (taxa.tsv, loci.tsv, summary.json) will be written, and the page only if the bundle fits")
    global MEM_RESIDENT_EST
    need, mx, mlists, mclosest, mfull, mslab, mrc = memory_estimate(n, L, loci, a.candidates, K_eff=K_eff, rc=not a.no_broad_sample); MEM_RESIDENT_EST = need
    parts = (f"{fmem(mx)} of alignment (the rows and the residue codes), {fmem(mlists)} of candidate lists for {n:,} taxa, {fmem(mclosest)} of closest matches for "
             f"{n:,} taxa x {len(loci):,} loci" + (f", {fmem(mfull)} of the full m matrix (kept up to {FULL_M_MAX_TAXA:,} taxa)" if mfull else "")
             + (f", {fmem(mslab)} of pair counts at every locus (kept up to {SWEEP_BLOCK_MIN:,} taxa)" if mslab else "")
             + (f", {fmem(mrc)} of broad-sample pools" if mrc else ""))
    if MEM_BUDGET is None: LOG(f"memory: about {fmem(need)} resident ({parts}); no budget (--memory 0)")
    else:
        LOG(f"memory: about {fmem(need)} resident ({parts}); budget {fmem(MEM_BUDGET)}" + ("" if a.memory else " (half of the memory available now; --memory sets it)"))
        if need > MEM_BUDGET:
            sys.exit(f"this run needs about {fmem(need)} of memory ({parts}), more than the budget of {fmem(MEM_BUDGET)}"
                     + ("" if a.memory else ", half of what is available now") + ": --memory raises the budget (--memory 0 removes it), or give fewer taxa or loci")

    P = prepare(a, names, rows, loci, dtype, K_eff, t0)
    X = P.X; run_z.cls, run_z.has_dup = P.cls, P.has_dup
    focal_sites, thin, present = P.focal_sites, P.thin, P.present
    masked = rows                                                      # blanked in place; X is its own copy

    seq_policy = a.report_sequences if a.report_sequences != "auto" else ("all" if n * L <= REPORT_SEQ_ALL_MAX else "fails")
    rep = dict(focal_sites=focal_sites, thin=thin, present=present, seq_policy=seq_policy if not a.no_report else "none",
               codes_policy=a.report_codes if not a.no_report else "none", no_report=a.no_report, uncovered_sites=nuncov, restart=P.restart)
    nmasked = run_z(a, names, loci, P, masked, K_eff, M_eff, rep)

    PROG.start("write", 2)
    out = write_masked(a.aln, a.out, fmt, names, masked, nexus)
    PROG.tick()
    page = write_report(a, names, loci, rep, dtype, fmt, K, K_eff, M_eff, t0) if not a.no_report else False
    PROG.finish()
    shutil.rmtree(os.path.join(a.out, "restart"), ignore_errors=True)   # every output is written: the restart points have served
    LOG(f"outputs: {out}, {a.out}/masked.csv" + (f", {a.out}/tests.tsv.gz" if a.all_tests else "")
        + ("" if a.no_report else f", {a.out}/report.html" if page else f", {a.out}/report.npz, {a.out}/taxa.tsv, {a.out}/loci.tsv, {a.out}/summary.json") + f", {a.out}/run.log"
        + f"\ndone {time.time() - t0:.0f}s")
    LOG.close()


# ------------------------------------------------------------------ the counting test
LSET_CHUNK = 256        # the fewest taxa per slice when a locus's closest matches are found (LSET_ELEMENTS sets the rest)
LN10 = math.log(10.0)
_FORK_FN = None
def _fork_call(x): return _FORK_FN(x)
def _can_fork():
    """worker processes are possible: fork exists, and this is not itself a pool worker (a daemonic process may not
    have children, and a rows-only sweep can be asked for from inside the test stage's workers)"""
    return "fork" in multiprocessing.get_all_start_methods() and not multiprocessing.current_process().daemon
def fork_map(fn, items, workers, chunksize=1):
    """fn over items in forked worker processes, the results in order, one at a time as they arrive (a generator: the
    workers run on while a result is consumed). The workers inherit fn and everything it closes over at the fork, so
    nothing is pickled on the way in; only the results come back. A plain map when fork is not available (Windows),
    with one worker, or for a single item. Threads cannot do this job: the per-locus work is Python around numpy
    calls of modest size, and threads only queue on the interpreter lock (Chen_2015: 15 s alone, 32 s with 8 threads)"""
    global _FORK_FN
    items = list(items); workers = min(workers, len(items))
    if workers <= 1 or not _can_fork():
        for x in items: yield fn(x)
        return
    _FORK_FN = fn
    try:
        with multiprocessing.get_context("fork").Pool(workers) as mp:
            for r in mp.imap(_fork_call, items, chunksize=chunksize): yield r
    finally: _FORK_FN = None


MASKED_COLUMNS = ["locus", "taxon", "reason", "mask_test", "mask_p", "residues_masked",
                  "close_relative_fails", "close_relative_comparators", "close_relative_observed", "close_relative_expected", "close_relative_excess_m", "close_relative_p",
                  "closest_match_partner", "closest_match_arm", "closest_match_role", "closest_match_observed", "closest_match_expected", "closest_match_row_factor", "closest_match_p", "closest_match_relatives", "closest_match_partner_excess_m",
                  "broad_sample_partner", "broad_sample_side", "broad_sample_observed", "broad_sample_expected", "broad_sample_p",
                  "comparators", "p_bar"]


def pfmt(lnp):
    """a p-value for the text outputs, from a NATURAL-log input: plain p at six decimals (1e-7 and below print
    0.000000 - the extreme end is deliberately flattened; the report's arrays keep the exact log-space values)"""
    return "" if lnp is None or not math.isfinite(lnp) else f"{min(math.exp(lnp), 1.0):.6f}"


def pfmt10(l10):
    """the same, from a log10 input"""
    return "" if l10 is None or not math.isfinite(l10) else f"{min(10.0 ** l10, 1.0):.6f}"


def fnum(v, nd):
    """a number for the text outputs: fixed decimals, no negative zero, empty when not finite"""
    v = float(v)
    return "" if not math.isfinite(v) else f"{round(v, nd) + 0.0:.{nd}f}"
TESTS_COLUMNS = ["test", "locus", "taxon", "comparator", "overlap", "mismatches", "expected", "row_factor", "p", "fail"]


def run_z(a, names, loci, P, masked, K_eff, M_eff, rep):
    n = len(names); nloc = len(loci); X, plateaus, rate, keep, lists = P.X, P.plateaus, P.rate, P.keep, P.lists
    mask_hi = not a.mask_similar_only
    mask_lo = not a.mask_different_only
    views = [np.frombuffer(mb, dtype=np.uint8) for mb in masked]        # writable views for blanking
    bs = None if getattr(a, "no_broad_sample", False) else bs_prepare(P, K_eff, a.threads, PROG)
    if bs is not None: LOG(f"broad sample: {bs.k} per sequence per locus, drawn from each taxon's pool of {bs.Pn}"
                           + ("" if P.m_full is not None else f" (the pools' pair divergences computed: {n * bs.Pn:,} pairs)"))

    # each taxon's candidates in order of pair divergence are its list from the sweep (the tie-break makes every row a
    # strict order): a locus's comparators are the first K candidates that share enough sites with the taxon at that
    # locus. A taxon whose list runs out at some locus gets its full row afterwards and is walked again
    Lbar = float(np.mean([len(st) for _, st in loci])) if nloc else 0.0
    cmp_cpu = nloc * ((5e-4 + 4e-6 * n) if isinstance(P.counts, SlabCounts) else 20.0 * n * Lbar / 1.5e8)     # seconds of walk, about
    cmp_workers = max(1, min(a.threads, nloc, int(cmp_cpu / COMPARATOR_SECONDS_PER_WORKER) + 1)) if a.threads > 1 else 1
    if a.threads > 1: LOG(f"comparators: {cmp_workers} worker process{'es' if cmp_workers > 1 else ''} ({n:,} taxa x {nloc:,} loci)")
    C_len = lists["cand"].shape[1]; nfin = lists["nfin"]
    cl_win, cl_dmin = closest_winners(*P.closest)                           # the closest match at every cell, from the sweep
    def row_cands(g):
        if g in P.ext: return P.ext[g][0]
        c = lists["cand"][g]; return c[c >= 0]
    def choose_comparators(li, only=None):
        """for every taxon at this locus: its first K_eff candidates (in order of m) that have the overlap the taxon
        needs here, walking its list in windows that grow until every taxon has K or the list ends; then the closest
        match at the locus and whether it is in the set. only: these taxa alone (the second pass)"""
        sites = loci[li][1]; Xs = X[:, sites]; nres = np.where(keep, P.focal_sites[:, li], 0).astype(np.int64); need = np.minimum(a.fit_overlap, nres)
        top = np.full((n, K_eff), -1, np.int64); cnt = np.zeros(n, np.int64); exhausted = np.zeros(n, bool)
        rows = np.flatnonzero(nres > 0) if only is None else np.asarray(only, np.int64)[nres[np.asarray(only, np.int64)] > 0]
        dup = getattr(run_z, "has_dup", False); cls = run_z.cls if dup else None
        seen_cls = np.full((n, K_eff), -1, np.int64) if dup else None          # with identical genomes in the data: the classes each row's comparators already cover
        ext_ids = np.array(list(P.ext), np.int64); longest = max([C_len] + [len(P.ext[g][0]) for g in P.ext]) if P.ext else C_len
        def window(rows, start, stop):
            """the candidates [start, stop) of these rows, from the lists or the extended rows, -1 beyond a row's end"""
            out = np.full((len(rows), stop - start), -1, np.int64)
            if start < C_len: seg = lists["cand"][rows, start:stop]; out[:, :seg.shape[1]] = seg
            if len(ext_ids):
                for i in np.flatnonzero(np.isin(rows, ext_ids)).tolist():
                    seg = P.ext[int(rows[i])][0][start:stop]; out[i] = -1; out[i, :len(seg)] = seg
            return out
        start = 0; width = 4 * K_eff
        while len(rows) and start < longest:
            stop = min(start + width, longest); cand = window(rows, start, stop); valid = cand >= 0
            ov, _ = P.counts.at(li, np.repeat(rows, stop - start), np.where(valid, cand, 0).ravel(), Xs); ov = ov.reshape(len(rows), -1)
            ok = valid & (ov >= need[rows][:, None]) & (cand != rows[:, None])
            if dup:                                                        # one comparator per class of copies: column by column, a class already covered is skipped
                k_ = cnt[rows].copy(); acc = seen_cls[rows]
                for j in range(cand.shape[1]):
                    cj = cls[np.maximum(cand[:, j], 0)]; good = ok[:, j] & (k_ < K_eff) & ~(acc == cj[:, None]).any(1); r_ = np.flatnonzero(good)
                    top[rows[r_], k_[r_]] = cand[r_, j]; acc[r_, k_[r_]] = cj[r_]; k_[r_] += 1
                cnt[rows] = k_; seen_cls[rows] = acc
            else:
                pos = cnt[rows][:, None] + np.cumsum(ok, axis=1) - ok             # the comparator slot each accepted candidate takes
                r_, j_ = np.nonzero(ok & (pos < K_eff))
                top[rows[r_], pos[r_, j_]] = cand[r_, j_]
                cnt[rows] = np.minimum(cnt[rows] + ok.sum(1), K_eff)
            rows = rows[(cnt[rows] < K_eff) & valid[:, -1]]                      # the list ended inside the window: nothing beyond
            start = stop; width *= 4
        short = np.flatnonzero((cnt < K_eff) & (nres > 0) & (nfin > C_len)); exhausted[short[~np.isin(short, ext_ids)]] = True
        if only is not None: exhausted[~np.isin(np.arange(n), np.asarray(only, np.int64))] = False
        # the closest match AT THE LOCUS from the sweep; it is tested (test 2) when it is not one of the K_eff closest
        # genomes (a tie for closest that includes one counts as in-set)
        closest = np.full(n, -1, np.int64); fork = np.zeros(n, bool); rows = np.flatnonzero(cnt > 0) if only is None else np.asarray(only, np.int64)[cnt[np.asarray(only, np.int64)] > 0]
        if len(rows):
            closest[rows] = cl_win[rows, li]; has = closest[rows] >= 0; rr = rows[has]
            if len(rr):
                Ct = top[rr]; valid = Ct >= 0
                ov, sm = P.counts.at(li, np.repeat(rr, K_eff), np.where(valid, Ct, 0).ravel(), Xs)
                ov = ov.reshape(len(rr), K_eff).astype(np.float64); sm = sm.reshape(len(rr), K_eff).astype(np.float64)
                with np.errstate(invalid="ignore", divide="ignore"):
                    d = np.where(valid & (ov >= a.fit_overlap), (ov - sm + HALF) / np.maximum(ov, 1.0), np.inf)
                inset = (d <= cl_dmin[rr, li][:, None]).any(1); fork[rr] = ~inset
        return top, top >= 0, cnt, nres, closest, fork, exhausted

    PROG.start("comparators", nloc)
    sel = []
    for r in fork_map(choose_comparators, range(nloc), cmp_workers): sel.append(r); PROG.tick()
    PROG.finish()
    exhausted = sorted({int(g) for r in sel for g in np.flatnonzero(r[6])})
    if exhausted:                                                        # the bounded second pass: full rows for the taxa whose list ran out, then their comparators again
        t1 = time.time(); E = sweep_pairs(rr=P.rr, ok_r=P.ok_r, want_closest=False, want_full=False, key="mj", rows_only=np.array(exhausted, np.int64), **{**P.sweep_args, "C": n})
        for k, g in enumerate(exhausted):
            c = E["lists"]["cand"][k]; fin = c >= 0; P.ext[g] = (c[fin], E["lists"]["m"][k][fin], E["lists"]["mfl"][k][fin], E["lists"]["key"][k][fin])
        for li in range(nloc):
            r = choose_comparators(li, only=exhausted)
            for old, new in zip(sel[li][:6], r[:6]): old[exhausted] = new[exhausted]
        LOG(f"comparators: {len(exhausted)} taxa ran out of their {C_len} candidates at some locus: full rows computed and walked ({time.time() - t1:.0f}s)")
    rep["lists_extended"] = dict(rate_neighbours=int(P.n_starved), comparators=len(exhausted))      # taxa whose candidate lists ran out and were extended to full rows
    Npairs = 0                                                        # comparisons over all loci: sizes the thread decision only
    for sel_ in sel: Npairs += int(sel_[1][sel_[2] > 0].sum())

    def locus_pairs(li):
        """the comparator pairs of a locus and the locus's rate; None when the locus is skipped"""
        top, valid, ncomp, nres, closest, fork = sel[li][:6]
        idx = np.flatnonzero(ncomp > 0)
        if len(idx) == 0: return "skipped-no-cells"
        gi, ci, rk = flatten_pairs(idx, top, valid, K_eff)
        Xs = X[:, loci[li][1]]; oi, si = P.counts.at(li, gi, ci, Xs); o = oi.astype(np.float32); same = si.astype(np.float32); mis = o - same
        ident = 100.0 * same / o; dob = 100.0 - ident
        mp, mfl_ = pair_m_lists(P, gi, ci); mph = np.maximum(mp, mfl_); mpx = np.maximum(mp - mfl_, 0.0)
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
                    ncomp=ncomp, nres=nres, expect=expect, corr=corr, rate_c=rate_c, Xs=Xs)

    # worker threads pay off in the test only when a locus's arrays are large (thousands of comparisons, or the pairs of
    # a thousand taxa and more): with small loci the per-call interpreter-lock traffic of many small numpy operations
    # makes threads slower than one
    ntested = sum(1 for t in sel if (t[2] > 0).any())
    # dispersion pre-pass: the robust spread of observed/expected among each locus's resolved comparator pairs with
    # enough expected mismatches, with the Poisson part removed: Var(O/E) = 1/E + phi. The same pass counts the tests
    # each locus will run (its distinct comparisons and its testable pairs), so the Bonferroni bar is known before any
    # test runs.
    def locus_phi(li):
        L = locus_pairs(li)
        if isinstance(L, str): return np.nan, 0.0, 0, 0
        gmin = np.minimum(L["gi"], L["ci"]).astype(np.int64); gmax = np.maximum(L["gi"], L["ci"]).astype(np.int64)
        ndist = len(np.unique(gmin * n + gmax))
        fk = np.flatnonzero(sel[li][5]); Y_ = sel[li][4][fk]; top_ = sel[li][0][fk]
        wo, _ = P.counts.at(li, np.where(top_ >= 0, top_, 0).ravel(), np.repeat(Y_, K_eff))
        okr = (top_ >= 0) & (top_ != Y_[:, None]) & (wo.reshape(top_.shape) >= a.fit_overlap)
        npt = int((okr.sum(1) >= 3).sum())
        E = L["o"] * L["expect"](L["mph"]) / 100.0
        use = L["fit"] & (E >= PHI_MIN_E)
        if use.sum() < 10: return np.nan, 0.0, ndist, npt
        r = np.log((L["mis"][use] + HALF) / E[use])
        s_ = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        return max(s_ ** 2 - float(np.median(1.0 / E[use])), 0.0), float(use.sum()), ndist, npt
    phi_hat = np.full(nloc, np.nan); phi_n = np.zeros(nloc); ndist_l = np.zeros(nloc, np.int64); npt_l = np.zeros(nloc, np.int64)
    PROG.start("dispersion", nloc)
    disp_workers = max(1, min(a.threads, nloc, Npairs // DISPERSION_PAIRS_PER_WORKER)) if a.threads > 1 else 1
    for li, (h, c, nd_, np_) in enumerate(fork_map(locus_phi, range(nloc), disp_workers, chunksize=max(1, nloc // (8 * disp_workers)))):
        phi_hat[li] = h; phi_n[li] = c; ndist_l[li] = nd_; npt_l[li] = np_; PROG.tick()
    PROG.finish()
    if int(ndist_l.sum()) == 0 and int(npt_l.sum()) == 0:
        # no test can run: say why, in the words of the thresholds, and stop before writing outputs that would mean nothing
        reasons = collections.Counter(locus_pairs(li) for li in range(nloc))
        reasons = collections.Counter(v if isinstance(v, str) else "tested" for v in reasons.elements())
        thin = rep["thin"]
        nk = int(keep.sum()); fin = int(P.fin_pairs); npairs = nk * (nk - 1) // 2; nrated = int(np.isfinite(rate).sum())
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
    D_tests = int(ndist_l.sum()); S_tests = int(npt_l.sum())
    BS_tests = int(((P.focal_sites > 0) & keep[:, None]).sum()) if bs is not None else 0    # one per present sequence of a kept taxon, frozen before any test
    N_tests = 2 * D_tests + S_tests + BS_tests
    logbar = math.log(a.alpha / max(N_tests, 1))
    kq = norm_quantile(1 - a.alpha / max(N_tests, 1))
    lbar10 = logbar / LN10
    nfork = int(sum(int(sel_[5].sum()) for sel_ in sel))
    LOG(f"tests: {D_tests:,} close-relative comparisons + {D_tests:,} closest-match comparisons (relative arm) + {S_tests:,} closest-match tests (outside arm: "
        f"{nfork:,} sequences whose closest match at the locus is not among their comparators, {S_tests:,} of them with three relatives "
        f"to set the expectation)" + (f" + {BS_tests:,} broad-sample tests (one per present sequence)" if bs is not None else "")
        + f" = N {N_tests:,}; every test is held to alpha / N = {a.alpha / max(N_tests, 1):.3g} (log10 {lbar10:.2f})")

    # ---- the report's per-cell columns: dense taxa x loci while the grid is small enough, else rows for cells with a fail
    dense = n * nloc <= DENSE_CELLS_MAX
    NUMCOLS = (("med_identity", np.float32, np.nan), ("med_expected", np.float32, np.nan), ("excess_m", np.float32, np.nan),
               ("closest", np.int32, -1), ("closest_overlap", np.uint32, 0), ("closest_identity", np.float32, np.nan), ("closest_expected", np.float32, np.nan),
               ("best_match", np.int32, -1), ("best_match_identity", np.float32, np.nan), ("partner", np.int32, -1),
               ("div_log10p", np.float32, np.nan), ("sim_log10p", np.float32, np.nan), ("locus_closest", np.int32, -1), ("locus_closest_log10p", np.float32, np.nan))
    for k in ("n_comp", "n_fail", "n_similar", "side"): rep[k] = np.zeros((n, nloc), np.uint8)
    rep["side2"] = np.zeros((n, nloc), np.uint8)       # the broad-sample test per cell: 1 tested, 2 past the bar, 4 shortfall at
    bsL = [None] * nloc                                # alpha, 8 excess at alpha, 16 excluded as a comparator here, 32 the shortfall is
                                                       # not reported (a too-similar mask covers the cell), 64 the divergent-cell label
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
        L = locus_pairs(li)
        if isinstance(L, str): return dict(skipped=True, log=f"{name}: {L}", status=L)
        idx, gi, ci, rk, o, same, mis, ident, dob, mp, mph, fit, ncomp, nres, expect, corr, rate_c = \
            (L[k] for k in ("idx", "gi", "ci", "rk", "o", "same", "mis", "ident", "dob", "mp", "mph", "fit", "ncomp", "nres", "expect", "corr", "rate_c"))
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
        cpos = {int(g): q for q, g in enumerate(idx.tolist())}; Xs = L["Xs"]; LR = None
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
            wo, ws = P.counts.at(li, C0.ravel(), np.repeat(Yv, K_eff), Xs)
            w = wo.reshape(C0.shape).astype(np.float64); ws = ws.reshape(C0.shape).astype(np.float64); valid &= w >= a.fit_overlap
            with np.errstate(invalid="ignore", divide="ignore"):
                rates = np.where(valid, (w - ws) / np.maximum(w, 1.0), np.nan)
            ok3 = valid.sum(1) >= 3
            if ok3.any():
                q3 = fq[ok3]; go, gs = P.counts.at(li, G[ok3], Yv[ok3], Xs)
                o2[q3] = go.astype(np.float64); mis2[q3] = o2[q3] - gs.astype(np.float64)
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

        def pex_of(other):
            qo = cpos.get(other); return fnum(exm[qo], 3) if qo is not None else ""
        def cm_cols(q):
            """the closest-match test's columns: the fired pair when it fired (either arm), else the strongest
            comparison it made - the outside test when it ran, else the most similar relative, with an empty role.
            The arm names the expectation: "relative" (the partner is one of the focal's relatives: the pair
            machinery, with the row factor) or "outside" (the empirical expectation from the relatives' own
            distances to the partner, which closest_match_relatives lists)"""
            if q in pair_of:
                rec, other, role = pair_of[q]; in_set = rec["arm"] == "comparator"
                return [names[other], "relative" if in_set else "outside", role,
                        fnum(100 * rec["mis"] / rec["o"], 2), fnum(100 * rec["E"] / rec["o"], 2),
                        fnum(rec["f"], 3) if in_set else "", pfmt(rec["lnp"]),
                        "" if in_set else "; ".join(names[c] for c in rec["rel"]), pex_of(other)]
            if tested2[q]:
                return [names[int(Y[q])], "outside", "", fnum(100 * mis2[q] / o2[q], 2), fnum(100 * float(E2[q]) / o2[q], 2),
                        "", pfmt(float(lp_lo[q])), "; ".join(names[c] for c in rel2[q]), pex_of(int(Y[q]))]
            s_, k_ = int(starts[q]), int(kc[q])
            if k_ == 0: return [""] * 9
            j = s_ + int(np.argmin(lp_in[s_:s_ + k_]))
            return [names[int(ci[j])], "relative", "", fnum(100 * mis[j] / o[j], 2), fnum(100 * E_lo[j] / o[j], 2),
                    fnum(ff[j], 3), pfmt(lp_in[j]), "", pex_of(int(ci[j]))]

        mrows_out = []
        for q in np.flatnonzero(masked_q).tolist():
            g = int(idx[q]); dv = bool(dflag[q] and mask_hi); sm = bool(sflag[q] and mask_lo)
            reason = "both" if dv and sm else ("divergent" if dv else "similar")
            mc = []                                                             # the cell's one significance: the most significant
            if dv and np.isfinite(div_lp[q]): mc.append((float(div_lp[q]), "close_relative"))   # bar-held p over the tests that masked it
            if sm and q in pair_of: mc.append((float(pair_of[q][0]["lnp"]), "closest_match"))
            mt, mp = (min(mc)[1], pfmt(min(mc)[0])) if mc else ("", "")
            s = int(starts[q]); comps = [names[int(c)] for c in ci[s:s + kc[q]]]
            row_ = [name, names[g], reason, mt, mp, int(nres[g]), int(nf[q]), int(kc[q]), fnum(100 - med_i[q], 2), fnum(100 - med_e[q], 2),
                    fnum(exm[q], 3), pfmt(div_lp[q])] + cm_cols(q) + ["", "", "", "", ""] + ["; ".join(comps), pfmt10(lbar10)]
            mrows_out.append((g, reason, row_))

        # ---- the report's view of this locus: cell columns, best matches, site codes, sequence rows, trees, pairs
        side_ = np.zeros(ncell, np.uint8)
        side_ |= np.where(dflag, 1 if mask_hi else 16, 0).astype(np.uint8); side_ |= (amb.astype(np.uint8) << 2)
        side_ |= np.where(sflag, 2 if mask_lo else 32, 0).astype(np.uint8); side_ |= np.where(tested2, 64, 0).astype(np.uint8)
        side_ |= np.where(src_q, 128, 0).astype(np.uint8); side_ |= np.where(unc_q, 8, 0).astype(np.uint8)
        # ---- the broad-sample test's stage 1 (the third test): every present cell of a kept taxon except the
        # divergent-masked, against bs.k comparators drawn by hash from the taxon's pool; expected counts through the
        # sequence's own rate rho; the statistic is the row's standardised-residual variance over the locus median.
        # The dataset-wide scale null, the bar and both stage-2 sides run in the parent once every locus
        # has reported. Divergent cells and masked copies are never comparators; a comparator whose own row variance
        # is extreme is dropped from every cell that drew it
        bs_blk = None; bsrep_full = np.full(n, -1, np.int64)      # every tested cell's broad representative, for the carried rows
        if bs is not None:
            dmask = np.zeros(n, bool); cpmask = np.zeros(n, bool)
            if mask_hi: dmask[idx[dflag]] = True
            if mask_lo:
                for q_, role in role_.items():
                    if role == "copy" and sflag[q_]: cpmask[int(idx[q_])] = True
            rows_bs = np.flatnonzero(keep & (P.focal_sites[:, li] > 0) & ~dmask)
            if len(rows_bs):
                P_ = plateaus[li]
                pos_bs = _hash_draw(bs_keys_cells(rows_bs, li), bs.Pn, bs.k)
                comp = bs.pool[rows_bs[:, None], pos_bs].astype(np.int64); mv0 = bs.mstar[rows_bs[:, None], pos_bs]
                mv = np.maximum(mv0, 1e-4)
                O_ = np.empty((len(rows_bs), bs.k)); W_ = np.empty_like(O_)
                for r0 in range(0, len(rows_bs), 512):
                    sl_ = slice(r0, min(r0 + 512, len(rows_bs)))
                    XG = Xs[rows_bs[sl_]]; XA = Xs[comp[sl_]]
                    both_ = (XG[:, None, :] >= 0) & (XA >= 0)
                    W_[sl_] = both_.sum(2); O_[sl_] = ((XG[:, None, :] != XA) & both_).sum(2)
                ok0 = (W_ >= a.fit_overlap) & np.isfinite(mv0) & keep[comp] & ~dmask[comp] & ~cpmask[comp]
                cexp = rate_c * corr(mv)
                def nmed(v, ok_):
                    out_ = np.full(v.shape[0], np.nan); some = ok_.any(1)
                    if some.any(): out_[some] = np.nanmedian(np.where(ok_, v, np.nan)[some], axis=1)
                    return out_
                def bs_var(ok_):
                    with np.errstate(invalid="ignore", divide="ignore"):
                        ratio = corr(100.0 * (O_ + HALF) / np.maximum(W_, 1.0)) / cexp
                        rho_ = nmed(ratio, ok_)
                    rho_ = np.where(np.isfinite(rho_) & (rho_ > 0), rho_, 1.0)
                    dhat_ = np.maximum(P_ * (1.0 - np.exp(-rho_[:, None] * cexp / P_)), 1e-4)
                    E_ = np.clip(W_ * dhat_ / 100.0, 1e-9, W_ * (1.0 - 1e-9))
                    with np.errstate(invalid="ignore", divide="ignore"):
                        z_ = np.where(ok_, (O_ - E_) / np.sqrt(E_ * (1.0 - dhat_ / 100.0) + phi * E_ * E_), np.nan)
                    kc_ = ok_.sum(1); s2 = np.full(len(rows_bs), np.nan)
                    two = kc_ >= 2
                    if two.any(): s2[two] = np.nanvar(z_[two], axis=1, ddof=1)
                    t_ = (kc_ >= BS_MIN) & np.isfinite(s2)
                    ref = float(np.median(s2[t_])) if t_.sum() >= 5 else np.nan
                    x_ = np.where(t_, s2 / ref, np.nan) if np.isfinite(ref) and ref > 0 else np.full(len(rows_bs), np.nan)
                    az = np.where(ok_ & np.isfinite(z_), np.abs(z_), -np.inf)           # the representative: the most extreme
                    jr = az.argmax(1)                                                   # standardised residual, either side
                    rep_ = np.where(np.isfinite(az[np.arange(len(rows_bs)), jr]), comp[np.arange(len(rows_bs)), jr], -1)
                    return x_, kc_, rho_, rep_
                x_bs, k_bs, rho_bs, rep_bs = bs_var(ok0)
                fin_bs = np.isfinite(x_bs); badt = np.zeros(n, bool)
                if fin_bs.sum() >= 5:
                    medx = float(np.median(x_bs[fin_bs]))
                    badt[rows_bs[fin_bs & (x_bs > BS_VAR_SCREEN * medx)]] = True
                    if badt.any(): x_bs, k_bs, rho_bs, rep_bs = bs_var(ok0 & ~badt[comp])
                fin_bs = np.isfinite(x_bs)
                bsrep_full[rows_bs[fin_bs]] = rep_bs[fin_bs]
                if fin_bs.any():
                    bs_blk = dict(rows=rows_bs[fin_bs].astype(np.int32), x=x_bs[fin_bs], k=np.minimum(k_bs[fin_bs], 32767).astype(np.int16),
                                  rho=rho_bs[fin_bs].astype(np.float32), rep=rep_bs[fin_bs].astype(np.int32),
                                  excl=np.flatnonzero(dmask | cpmask).astype(np.int32), varbad=np.flatnonzero(badt).astype(np.int32))
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
            LR = LR or P.counts.locus(li, Xs); wrow, srow = LR.rows(rows_bm)
            okc = wrow >= BEST_MATCH_MIN_SITES; okc[np.arange(len(rows_bm)), rows_bm] = False
            with np.errstate(invalid="ignore", divide="ignore"):
                dvg = np.where(okc, (wrow - srow) / np.maximum(wrow, 1.0), np.inf)
            mcache = {}
            def mb(r_, cols):
                """m between the r_-th row's taxon and the taxa cols, inf where there is none: asked for only where the
                tie-breaks below reach it (above FULL_M_MAX_TAXA each is computed, and remembered)"""
                g = int(rows_bm[r_]); cols = np.asarray(cols, np.int64)
                if P.m_full is not None: v = P.m_full[g, cols]
                else:
                    miss = sorted({int(t) for t in cols.tolist() if (g, int(t)) not in mcache})
                    if miss: mcache.update(zip(((g, t) for t in miss), pair_m(P, np.full(len(miss), g), np.array(miss, np.int64))[0].tolist()))
                    v = np.array([mcache[(g, int(t))] for t in cols.tolist()], np.float64)
                return np.where(np.isfinite(v), v, np.inf)
            wmax = np.where(okc, wrow, -1.0).max(1)
            dtier = np.where(okc & (wrow >= wmax[:, None]), dvg, np.inf); lo_t = dtier.min(1)
            tied = dtier == lo_t[:, None]; champ = tied.argmax(1)                              # the first tied; more than one: the closest by m
            for r_ in np.flatnonzero(np.isfinite(lo_t) & (tied.sum(1) > 1)).tolist():
                tc = np.flatnonzero(tied[r_]); tk = mb(r_, tc)
                if np.isfinite(tk).any(): champ[r_] = tc[tk.argmin()]                          # none of the tied has an m: the first tied
            okb = np.isfinite(lo_t)
            Xs_ = Xs if Xs is not None else X[:, sites]
            for r_ in np.flatnonzero(okb & (dvg < lo_t[:, None]).any(1)).tolist():
                g = int(rows_bm[r_]); fx = Xs_[g]; fp = fx >= 0; c = int(champ[r_])
                T0 = np.flatnonzero(dvg[r_] < lo_t[r_]); kw = -wrow[r_, T0]; kd = dvg[r_, T0]; km = np.zeros(len(T0))
                _, gid, gsz = np.unique(np.stack([kw, kd], 1), axis=0, return_inverse=True, return_counts=True); gid = gid.ravel()
                tie = gsz[gid] > 1
                if tie.any(): km[tie] = mb(r_, T0[tie])                                          # m only breaks ties of coverage and divergence
                T = T0[np.lexsort((km, kd, kw))]                                                 # broadest coverage first (a stable sort, as sorted())
                while len(T):                                                                  # every remaining challenger against the champion at once;
                    XT = Xs_[T]; common = fp[None, :] & (Xs_[c] >= 0)[None, :] & (XT >= 0)      # the first that wins takes over, and the rest face it
                    dt_ = ((XT != fx[None, :]) & common).sum(1); dc_ = ((Xs_[c] != fx)[None, :] & common).sum(1)
                    eqd = dt_ == dc_; mt = np.full(len(T), np.inf); mc_ = np.inf
                    if eqd.any(): mt[eqd] = mb(r_, T[eqd]); mc_ = mb(r_, [c])[0]
                    wins = (common.sum(1) >= BEST_MATCH_MIN_SITES) & ((dt_ < dc_) | (eqd & (mt < mc_)))
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
            cells_t = []
            for q in np.flatnonzero(masked_q).tolist():
                g = int(idx[q]); tips = [g] + [int(c) for c in top[g] if c >= 0]
                if partner_[q] >= 0 and int(partner_[q]) not in tips: tips.append(int(partner_[q]))
                if len(tips) >= 2: cells_t.append((g, tips))
            Wu = mu = None
            if cells_t and P.m_full is not None:                                   # the outgroup search's rows of counts and of m for every tip at the locus, in one call
                u = np.array(sorted({t for _, tips in cells_t for t in tips})); upos = {int(t): k for k, t in enumerate(u)}
                LR = LR or P.counts.locus(li, Xs); Wu = LR.rows(u)[0]; mu = m_rows(P, u)
            for g, tips in cells_t:
                og = -1; cand = np.ones(n, bool); cand[tips] = False; cidx = np.flatnonzero(cand)
                if len(cidx) and Wu is not None:                                   # above FULL_M_MAX_TAXA there is no outgroup search: the trees are rooted at their midpoint
                    ti = [upos[t] for t in tips]; Wsub = Wu[ti].T[cidx]; msub = mu[ti].T[cidx]
                    okog = (Wsub >= a.fit_overlap).all(1) & np.isfinite(msub).all(1)
                    if okog.any(): og = int(cidx[okog][np.argmax(msub[okog].mean(1))]); tips.append(og)
                tarr = np.array(tips); nt = len(tips); rows_t = Xs_[tarr]; pres_t = rows_t >= 0
                both = pres_t[:, None, :] & pres_t[None, :, :]; nb = both.sum(2)
                with np.errstate(invalid="ignore", divide="ignore"):
                    Do = np.where(nb > 0, 100.0 * ((rows_t[:, None, :] != rows_t[None, :, :]) & both).sum(2) / np.maximum(nb, 1), np.nan)
                    mv = pair_m(P, np.repeat(tarr, nt), np.tile(tarr, nt))[0].reshape(nt, nt); De =np.where(np.isfinite(mv), expect(np.where(np.isfinite(mv), mv, 1.0)), np.nan)
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
                    if bsrep_full[g] >= 0: rs.add(int(bsrep_full[g]))
                rows_s = np.array(sorted(rs), dtype=np.int64)
            if len(rows_s):
                seq_out = (rows_s.astype(np.int32), np.stack([views[g][sites] for g in rows_s.tolist()]))
        sim_json = ([dict(test="relative" if rec["arm"] == "comparator" else "outside", taxon=g, partner=y, overlap=int(rec["o"]), mismatches=int(rec["mis"]), expected=round(rec["E"], 3),
                          row_factor=round(rec["f"], 4), relatives=rec["rel"], log10p=round(rec["lnp"] / LN10, 3), fail=True,
                          copy=(-1 if rec.get("uncertain") else int(g if role_.get(cpos.get(g)) == "copy" else y)), uncertain=bool(rec.get("uncertain", False)),
                          excess_m=[(round(float(exm[cpos[x]]), 4) if x in cpos and np.isfinite(exm[cpos[x]]) else None) for x in (g, y)])
                     for (g, y), rec in pairs.items()]
                    + [dict(test="outside", taxon=int(idx[q]), partner=int(Y[q]), overlap=int(o2[q]), mismatches=int(mis2[q]), expected=round(float(E2[q]), 3),
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
            p6 = lambda v: (np.round(np.exp(np.minimum(np.asarray(v, float), 0.0)), 6) + 0.0).tolist()
            tests_block = ("close_relative\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t\t%.6f\t%d\n") * len(gi) % tuple(chain.from_iterable(zip(
                [tnames[x] for x in gi.tolist()], [tnames[x] for x in ci.tolist()], o.astype(np.int64).tolist(), mis.astype(np.int64).tolist(),
                r3(E), p6(lp_hi), fail_hi.astype(np.int64).tolist())))
            tests_block += ("closest_match_relative\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t%.4f\t%.6f\t%d\n") * len(gi) % tuple(chain.from_iterable(zip(
                [tnames[x] for x in gi.tolist()], [tnames[x] for x in ci.tolist()], o.astype(np.int64).tolist(), mis.astype(np.int64).tolist(),
                r3(E_lo), ff.tolist(), p6(lp_in), fail_in.astype(np.int64).tolist())))
            tq_ = np.flatnonzero(tested2)
            if len(tq_):
                tests_block += ("closest_match_outside\t" + lname + "\t%s\t%s\t%d\t%d\t%.3f\t\t%.6f\t%d\n") * len(tq_) % tuple(chain.from_iterable(zip(
                    [tnames[int(idx[q])] for q in tq_.tolist()], [tnames[int(Y[q])] for q in tq_.tolist()], o2[tq_].astype(np.int64).tolist(),
                    mis2[tq_].astype(np.int64).tolist(), r3(E2[tq_]), p6(lp_lo[tq_]), fail2[tq_].astype(np.int64).tolist())))

        log = (f"{name}: {ncell} cells, {len(gi)} comparator slots, {int(tested2.sum())} closest-match tests; too divergent {int(flag.sum())}"
               f" (+{int(amb.sum())} with fewer comparators), too-similar pairs {len(pairs)} ({sum(1 for r in pairs.values() if r['arm'] == 'closest')} by the closest match)")
        return dict(skipped=False, mrows=mrows_out, tests=tests_block, log=log, rep=repblk, phi=phi, rate_c=rate_c, rc=bs_blk)

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
        if R.get("rc") is not None:
            b = R["rc"]; bsL[li] = b
            rep["side2"][b["rows"], li] |= 1; rep["side2"][b["excl"], li] |= 16
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
    for lis, Rs in zip(batches, fork_map(run_batch, batches, workers)):   # results in locus order, consumed as they arrive
        for li, R in zip(lis, Rs): consume(li, R)
    PROG.finish()
    if tests_fh is not None: tests_fh.close()

    # ---- the broad-sample verdicts: the scale null over every tested row dataset-wide, the run's single bar,
    # then stage 2 for the rows past it (and the divergent-cell labels) in a second forked pass over only the loci
    # involved. A mask is a far shortfall: a coalescence younger than the species split allows (a transfer or
    # contamination: too close, the viewer's blue); a divergent cell with the shortfall gets a label (yellow)
    bs_summary = None
    if bs is not None and any(b is not None for b in bsL):
        loc_a = np.concatenate([np.full(len(b["rows"]), li, np.int64) for li, b in enumerate(bsL) if b is not None])
        tax_a = np.concatenate([b["rows"] for b in bsL if b is not None]).astype(np.int64)
        x_a = np.concatenate([b["x"] for b in bsL if b is not None])
        k_a = np.concatenate([b["k"] for b in bsL if b is not None])
        rho_a = np.concatenate([b["rho"] for b in bsL if b is not None])
        PROG.start("broad-sample null", 1)
        lp1, bs_tau = bs_null(x_a, k_a, threads=a.threads)
        PROG.tick(); PROG.finish()
        past = lp1 < logbar
        rep["side2"][tax_a[past], loc_a[past]] |= 2
        if dense:
            rep["bs_log10p"] = np.full((n, nloc), np.nan, np.float32); rep["bs_log10p"][tax_a, loc_a] = (lp1 / LN10).astype(np.float32)
            rep["bs_rho"] = np.full((n, nloc), np.nan, np.float32); rep["bs_rho"][tax_a, loc_a] = rho_a
            rep["bs_stat"] = np.full((n, nloc), np.nan, np.float32); rep["bs_stat"][tax_a, loc_a] = x_a.astype(np.float32)   # observed/expected variance: the stage-1 statistic itself
            rep["bs_rep"] = np.full((n, nloc), -1, np.int32)                                                              # the stage-1 representative of every tested cell:
            rep["bs_rep"][tax_a, loc_a] = np.concatenate([b["rep"] for b in bsL if b is not None])                        # the comparator with the most extreme standardised residual
        key_a = np.sort(loc_a * np.int64(n) + tax_a)                           # stage-1 p lookup for the report's sparse rows
        lp1_by = lp1[np.argsort(loc_a * np.int64(n) + tax_a, kind="stable")]
        flag_by_loc = {int(li): tax_a[past & (loc_a == li)] for li in np.unique(loc_a[past])}
        div_by_loc = {}
        if mask_hi:
            for li in range(nloc):
                if rep["test_rate"] is not None and np.isfinite(rep["test_rate"][li]):
                    dv = np.flatnonzero((rep["side"][:, li] & 1) > 0)
                    if len(dv): div_by_loc[li] = dv
        todo = sorted(set(flag_by_loc) | set(div_by_loc))
        la_bs = math.log(a.alpha)

        def bs_locus_b(li):
            """stage 2 at one locus, for the rows past the bar and the divergent-masked cells: the far shortfall with
            the comparator screen, exactly the per-locus pass's draws and fit. The fork also decides what is reported
            and masked here, and builds the masked cells' viewer artifacts (the codes row, the face-to-face tree with
            the partner as a tip, and the too-similar record)"""
            rows_f = np.asarray(flag_by_loc.get(li, []), np.int64); rows_d = np.asarray(div_by_loc.get(li, []), np.int64)
            rate_c = float(rep["test_rate"][li])
            if not np.isfinite(rate_c) or (len(rows_f) == 0 and len(rows_d) == 0): return li, None
            P_ = plateaus[li]; phi = float(phi_l[li]); sites = loci[li][1]; Xs = X[:, sites]
            corr = lambda p: -P_ * np.log(1.0 - np.minimum(np.asarray(p, float), 0.98 * P_) / P_)
            expect_ = lambda mv_: np.maximum(P_ * (1.0 - np.exp(-rate_c * corr(np.maximum(np.asarray(mv_, float), 1e-4)) / P_)), 1e-4)
            excl16 = (rep["side2"][:, li] & 16) > 0
            vb = np.zeros(n, bool); b_ = bsL[li]
            if b_ is not None: vb[b_["varbad"]] = True
            def draw(rows_):
                pos = _hash_draw(bs_keys_cells(rows_, li), bs.Pn, bs.k)
                comp = bs.pool[rows_[:, None], pos].astype(np.int64); mv0 = bs.mstar[rows_[:, None], pos]
                O_ = np.empty((len(rows_), bs.k)); W_ = np.empty_like(O_)
                for r0 in range(0, len(rows_), 512):
                    sl_ = slice(r0, min(r0 + 512, len(rows_)))
                    XG = Xs[rows_[sl_]]; XA = Xs[comp[sl_]]
                    both_ = (XG[:, None, :] >= 0) & (XA >= 0)
                    W_[sl_] = both_.sum(2); O_[sl_] = ((XG[:, None, :] != XA) & both_).sum(2)
                return comp, np.maximum(mv0, 1e-4), np.isfinite(mv0) & (W_ >= a.fit_overlap) & keep[comp] & ~excl16[comp] & ~vb[comp], O_, W_
            def nmed(v, ok2):
                out_ = np.full(v.shape[0], np.nan); some = ok2.any(1)
                if some.any(): out_[some] = np.nanmedian(np.where(ok2, v, np.nan)[some], axis=1)
                return out_
            def far_tails(rows_, use_rho, screen_lo, screen_hi):
                comp, mv, ok_, O_, W_ = draw(rows_)
                cexp = rate_c * corr(mv)
                def one(screen, upper):
                    ok2 = ok_ if screen is None else ok_ & ~screen[comp]
                    with np.errstate(invalid="ignore", divide="ignore"):
                        rho_ = nmed(corr(100.0 * (O_ + HALF) / np.maximum(W_, 1.0)) / cexp, ok2)
                    rho_ = np.where(np.isfinite(rho_) & (rho_ > 0), rho_, 1.0)
                    if not use_rho: rho_ = np.ones(len(rows_))
                    dhat_ = np.maximum(P_ * (1.0 - np.exp(-rho_[:, None] * cexp / P_)), 1e-4)
                    E_ = np.clip(W_ * dhat_ / 100.0, 1e-9, W_ * (1.0 - 1e-9))
                    med = nmed(mv, ok2)
                    far = ok2 & (mv > np.where(np.isfinite(med), med, np.inf)[:, None]); nf = far.sum(1)
                    lp = np.full(far.shape, np.inf); ii, jj = np.nonzero(far)
                    if len(ii): lp[ii, jj] = bb_tails(O_[ii, jj], W_[ii, jj], E_[ii, jj], np.full(len(ii), phi), upper)
                    with np.errstate(invalid="ignore", divide="ignore"):                            # the standardised residual of every far member:
                        z_ = np.where(far, (O_ - E_) / np.sqrt(E_ * (1.0 - dhat_ / 100.0) + phi * E_ * E_), np.nan)   # the representative is the most extreme aligned one
                    return dict(lp=lp, nf=nf, E=E_, rho=rho_, z=z_)
                lo_ = one(screen_lo, False)
                hi_ = one(screen_hi, True) if screen_hi is not False else None
                return comp, O_, W_, lo_, hi_
            def min_stat(d):
                return np.where(d["nf"] > 0, np.minimum(d["lp"].min(1) + np.log(np.maximum(d["nf"], 1)), 0.0), 0.0)
            def excess_stat(d):
                srt = np.sort(d["lp"], 1); kk = np.maximum(np.ceil(d["nf"] / 2.0), 1.0).astype(np.int64)
                pk = srt[np.arange(len(kk)), np.minimum(kk, d["lp"].shape[1]) - 1]
                return order_stat_bound(pk, d["nf"], kk)
            badS = np.zeros(n, bool); badE = np.zeros(n, bool)
            if len(rows_f):                                                    # the flagged rows' comparators, screened by their own rows
                pos_f = _hash_draw(bs_keys_cells(rows_f, li), bs.Pn, bs.k)
                anc = np.unique(bs.pool[rows_f[:, None], pos_f].astype(np.int64).ravel())
                anc = anc[keep[anc] & ~excl16[anc] & ~vb[anc]]
                if len(anc):
                    _, _, _, alo, ahi = far_tails(anc, True, None, None)
                    badS[anc[min_stat(alo) < BS_SCREEN_LOGP]] = True
                    badE[anc[min_stat(ahi) < BS_SCREEN_LOGP]] = True
            top = sel[li][0]
            out = dict(mask=[], note=[], label=[], lean=[], dlean=[], codes=None, trees=None, broad=[])
            tree_cells = []
            if len(rows_f):
                comp2, O2, W2, LO, HI = far_tails(rows_f, True, badS, badE)
                sp = min_stat(LO); ep = excess_stat(HI); ar_ = np.arange(len(rows_f))
                jlo = np.where(np.isfinite(LO["z"]), LO["z"], np.inf).argmin(1)         # the representatives: the most extreme
                plo = np.where(LO["nf"] > 0, comp2[ar_, jlo], -1)                       # standardised residual on each side's far set
                jhi = np.where(np.isfinite(HI["z"]), HI["z"], -np.inf).argmax(1)
                phi_c = np.where(HI["nf"] > 0, comp2[ar_, jhi], -1)
                for i, g in enumerate(rows_f.tolist()):
                    close_call = bool(sp[i] <= ep[i])                                   # one direction at most: the stronger side is the call,
                    D_ = LO if close_call else HI                                        # significant or not; the other side stays in the alignment rows
                    stat = float(sp[i] if close_call else ep[i])
                    j_ = int(jlo[i] if close_call else jhi[i]); partner = int(plo[i] if close_call else phi_c[i])
                    out["lean"].append((g, float(sp[i]), float(ep[i]), partner, 0 if partner < 0 else (1 if close_call else -1)))
                    if not stat < la_bs or partner < 0: continue
                    side_txt = "close" if close_call else "divergent"
                    sd0 = int(rep["side"][g, li])
                    obs_ = 100.0 * O2[i, j_] / max(W2[i, j_], 1.0); exp_ = 100.0 * float(D_["E"][i, j_]) / max(W2[i, j_], 1.0)
                    rec_ = dict(side=side_txt, taxon=int(g), partner=partner, overlap=int(W2[i, j_]), mismatches=int(O2[i, j_]),
                                expected=round(float(D_["E"][i, j_]), 3), rho=round(float(D_["rho"][i]), 4), log10p=round(stat / LN10, 3))
                    out["broad"].append(rec_)
                    b2 = 4 if close_call else 8
                    if sd0 & 3:                                                          # a mask already covers the cell: the firing is reported
                        out["note"].append((g, b2 | 32, side_txt, stat, partner, obs_, exp_)); continue   # in its columns (a divergent call still turns the reason "both")
                    if not (mask_lo if close_call else mask_hi):
                        out["note"].append((g, b2, side_txt, stat, partner, obs_, exp_)); continue
                    out["mask"].append((g, b2, side_txt, stat, partner, obs_, exp_))
                    tree_cells.append((g, partner, close_call))
            if len(rows_d):                                                    # the divergent cells' label: no rho, held at the bar
                comp3, Od, Wd, DLO, _ = far_tails(rows_d, False, badS, False)
                dp = min_stat(DLO)
                jd = np.where(np.isfinite(DLO["lp"]), DLO["lp"], np.inf).argmin(1)
                pld = np.where(DLO["nf"] > 0, comp3[np.arange(len(rows_d)), jd], -1)
                for i, g in enumerate(rows_d.tolist()):
                    if dp[i] < logbar:
                        jd_ = int(np.where(np.isfinite(DLO["lp"][i]), DLO["lp"][i], np.inf).argmin())
                        out["label"].append((g, float(dp[i]), int(pld[i]), 100.0 * Od[i, jd_] / max(Wd[i, jd_], 1.0), 100.0 * float(DLO["E"][i, jd_]) / max(Wd[i, jd_], 1.0)))
                    else:                                                        # no label: the cell still names its representative,
                        zi = DLO["z"][i]; fin_ = np.isfinite(zi)                 # the far member with the most extreme standardised residual
                        out["dlean"].append((g, float(dp[i]), int(comp3[i, int(np.where(fin_, zi, np.inf).argmin())]) if fin_.any() else -1))
            # the masked cells' artifacts, exactly the per-locus pass's constructions
            if tree_cells:
                qs = np.array([g for g, _, _ in tree_cells], np.int64)
                FX = Xs[qs]; pres_ = FX >= 0
                Bq = np.array([p if cc else -1 for _, p, cc in tree_cells], np.int64)   # only a close-side partner wears the agreement marks
                C_ = top[qs]; cv = (C_ >= 0) & (C_ != Bq[:, None]); XC = Xs[np.maximum(C_, 0)]
                k_ = ((XC != FX[:, None, :]) & (XC >= 0) & cv[:, :, None]).sum(1)
                rows_c = np.where(pres_, np.minimum(k_, 127), CODE_ABSENT).astype(np.uint8)
                hb = Bq >= 0
                if hb.any():
                    XB = Xs[np.maximum(Bq, 0)]; shared = hb[:, None] & pres_ & (XB >= 0) & (FX == XB) & (k_ >= 2)
                    rows_c[shared] = (CODE_SIMILAR + np.minimum(k_[shared], 125)).astype(np.uint8)
                out["codes"] = (qs.astype(np.int32), rows_c)
                tcells, texp, tobs, ttips, togs, tclos = [], [], [], [], [], []
                cells_t = []
                for g, partner, _cc in tree_cells:
                    tips = [int(g)] + [int(c) for c in top[g] if c >= 0]
                    if partner >= 0 and int(partner) not in tips: tips.append(int(partner))
                    if len(tips) >= 2: cells_t.append((int(g), tips))
                Wu = mu = None
                if cells_t and P.m_full is not None:
                    u = np.array(sorted({t for _, tips in cells_t for t in tips})); upos = {int(t): k2_ for k2_, t in enumerate(u)}
                    LR = P.counts.locus(li, Xs); Wu = LR.rows(u)[0]; mu = m_rows(P, u)
                for g, tips in cells_t:
                    og = -1; cand = np.ones(n, bool); cand[tips] = False; cidx = np.flatnonzero(cand)
                    if len(cidx) and Wu is not None:
                        ti_ = [upos[t] for t in tips]; Wsub = Wu[ti_].T[cidx]; msub = mu[ti_].T[cidx]
                        okog = (Wsub >= a.fit_overlap).all(1) & np.isfinite(msub).all(1)
                        if okog.any(): og = int(cidx[okog][np.argmax(msub[okog].mean(1))]); tips.append(og)
                    tarr = np.array(tips); nt = len(tips); rows_t = Xs[tarr]; pres_t = rows_t >= 0
                    both = pres_t[:, None, :] & pres_t[None, :, :]; nb = both.sum(2)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        Do = np.where(nb > 0, 100.0 * ((rows_t[:, None, :] != rows_t[None, :, :]) & both).sum(2) / np.maximum(nb, 1), np.nan)
                        mv_ = pair_m(P, np.repeat(tarr, nt), np.tile(tarr, nt))[0].reshape(nt, nt)
                        De = np.where(np.isfinite(mv_), expect_(np.where(np.isfinite(mv_), mv_, 1.0)), np.nan)
                    np.fill_diagonal(De, 0.0); np.fill_diagonal(Do, 0.0)
                    if np.isnan(De).any(): De = np.where(np.isnan(De), np.nanmax(De) if np.isfinite(np.nanmax(De)) else 50.0, De)
                    if np.isnan(Do).any(): Do = np.where(np.isnan(Do), De, Do)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        De = corr(De); Do = corr(Do)
                    np.fill_diagonal(De, 0.0); np.fill_diagonal(Do, 0.0)
                    t_ = cell_trees(De, Do, [str(x) for x in tips], tips.index(og) if og >= 0 else None)
                    if t_:
                        cand = np.ones(nt, bool); cand[0] = False
                        if og >= 0: cand[tips.index(og)] = False
                        ce = int(np.argmin(np.where(cand, De[0], np.inf))) if cand.any() else -1; co = int(np.argmin(np.where(cand, Do[0], np.inf))) if cand.any() else -1
                        tcells.append(g); texp.append(t_[0]); tobs.append(t_[1]); ttips.append(tips); togs.append(tips.index(og) if og >= 0 else -1); tclos.append([ce, co])
                if tcells: out["trees"] = dict(cells=tcells, expected=texp, observed=tobs, tips=ttips, outgroup=togs, closest=tclos)
                if rep["seq_policy"] not in ("none", "all"):
                    rs_ = set()
                    for g, partner, _cc in tree_cells:
                        rs_.add(int(g)); rs_.update(int(c) for c in top[g] if c >= 0)
                        if partner >= 0: rs_.add(int(partner))
                    for g, _dp, rp_, _o, _e in out["label"]:                    # the divergent cells' broad rows: the label's partner
                        if rp_ >= 0: rs_.add(int(rp_))
                    for g, _dp, rp_ in out["dlean"]:                            # or the representative without one
                        if rp_ >= 0: rs_.add(int(rp_))
                    if rs_:
                        rows_s = np.array(sorted(rs_), np.int64)
                        out["seq"] = (rows_s.astype(np.int32), np.stack([views[g_][sites] for g_ in rows_s.tolist()]))
            return li, out

        mrow_at = {(r[0], r[1]): j for j, r in enumerate(mrows)}
        sparse_at = {}
        if not dense:
            for piece_l, piece_t in zip(sparse_loc, sparse_tax):
                for l_, t_ in zip(piece_l.tolist(), piece_t.tolist()): sparse_at.setdefault(l_, set()).add(t_)
        bs_rows_rep = []; tally_bs = collections.Counter()
        if todo:
            bsw = min(workers, len(todo)) if a.threads > 1 else 1
            PROG.start("broad sample", len(todo))
            for li, out in fork_map(bs_locus_b, todo, bsw):
                PROG.tick()
                if out is None: continue
                name, sites = loci[li]
                def var_lp10(li_, g_):
                    """the cell's stage-1 variance p (log10): the broad test's one bar-held number"""
                    k_ = np.int64(li_) * np.int64(n) + np.int64(g_); p_ = int(np.searchsorted(key_a, k_))
                    return float(lp1_by[p_]) / LN10 if p_ < len(key_a) and key_a[p_] == k_ else float("nan")
                for g, sp_i, ep_i, partner, call_i in out["lean"]:              # every gate-passer: both sides' statistics, the call and
                    bs_rows_rep.append((li, g, sp_i, ep_i, partner, call_i))     # its representative, significant or not (the viewer's lean)
                for g, b2, side_txt, lp_use, partner, obs_, exp_ in out["note"]:
                    rep["side2"][g, li] |= b2; tally_bs["noted"] += 1
                    if b2 & 8 and (b2 & 32) and mask_hi and (rep["side"][g, li] & 3) == 2:   # masked too similar already, now also too
                        rep["side"][g, li] |= 1                                              # divergent from the broad sample: the directions differ
                        j = mrow_at.get((name, names[g]))
                        if j is not None and mrows[j][2] == "similar": tally["similar"] -= 1; tally["both"] += 1; mrows[j][2] = "both"
                    j = mrow_at.get((name, names[g]))
                    if j is not None:
                        mrows[j][21] = names[partner] if partner >= 0 else ""; mrows[j][22] = side_txt
                        mrows[j][23] = fnum(obs_, 2); mrows[j][24] = fnum(exp_, 2); mrows[j][25] = pfmt(lp_use)
                        v10 = var_lp10(li, g)                                   # the broad firing joins the cell's one significance when it is the most significant
                        if np.isfinite(v10) and (mrows[j][4] == "" or 10.0 ** v10 < float(mrows[j][4]) - 1e-12):
                            mrows[j][3] = "broad_sample"; mrows[j][4] = pfmt10(v10)
                for g, b2, side_txt, lp_use, partner, obs_, exp_ in out["mask"]:
                    rep["side2"][g, li] |= b2; rep["side"][g, li] |= 2 if b2 & 4 else 1
                    views[g][sites] = 45
                    reason = "similar" if b2 & 4 else "divergent"
                    tally[reason] += 1; tally_bs["close" if b2 & 4 else "divergent"] += 1
                    mrows.append([name, names[g], reason, "broad_sample", pfmt10(var_lp10(li, g)), int(P.focal_sites[g, li])] + [""] * 6 + [""] * 9
                                 + [names[partner] if partner >= 0 else "", side_txt, fnum(obs_, 2), fnum(exp_, 2), pfmt(lp_use)]
                                 + ["", pfmt10(lbar10)])
                    mrow_at[(name, names[g])] = len(mrows) - 1
                    if dense:
                        rep["partner"][g, li] = partner
                    elif g not in sparse_at.get(li, ()):                        # the masked cell joins the sparse rows with its partner
                        sparse_loc.append(np.full(1, li, np.int32)); sparse_tax.append(np.full(1, g, np.int32))
                        for k_, dt_, fill_ in NUMCOLS:
                            v_ = np.full(1, fill_, dt_)
                            if k_ == "partner": v_[0] = partner
                            rep[k_].append(v_)
                        rep["comparators"].append(sel[li][0][g:g + 1].astype(np.int32))
                        for k_ in ("comp_div_log10p", "comp_sim_log10p"): rep[k_].append(np.full((1, K_eff), np.nan, np.float32))
                        sparse_at.setdefault(li, set()).add(g)
                for g, dp_i, rp_ in out["dlean"]:                              # every divergent cell's representative, label or not
                    bs_rows_rep.append((li, g, dp_i, np.nan, rp_, 0))
                for g, dp_i, partner, obs_, exp_ in out["label"]:
                    rep["side2"][g, li] |= 64
                    bs_rows_rep.append((li, g, dp_i, np.nan, partner, 1))
                    j = mrow_at.get((name, names[g]))
                    if j is not None:                                           # a label on the divergent row: the mask, its reason and its significance unchanged
                        mrows[j][21] = names[partner] if partner >= 0 else ""; mrows[j][22] = "close"
                        mrows[j][23] = fnum(obs_, 2); mrows[j][24] = fnum(exp_, 2); mrows[j][25] = pfmt(dp_i)
                        tally_bs["labels"] += 1
                if out["broad"]:
                    rep.setdefault("bs_records", {}).setdefault(li, []).extend(out["broad"])
                if out["codes"] is not None:
                    qs_, blk_ = out["codes"]
                    if li in rep["codes"]:                                      # a cell the single pass carried (a stray failing comparator) is replaced: its row now knows the partner
                        c0, b0 = rep["codes"][li]; keep_ = ~np.isin(c0, qs_)
                        rep["codes"][li] = (np.concatenate([c0[keep_], qs_]), np.concatenate([b0[keep_], blk_]))
                    else: rep["codes"][li] = (qs_, blk_)
                if out.get("seq") is not None:
                    rows_s, blk_s = out["seq"]
                    if li in rep["seqs"]:
                        r0, b0 = rep["seqs"][li]; keep_ = ~np.isin(r0, rows_s)
                        o_ = np.argsort(np.concatenate([r0[keep_], rows_s]), kind="stable")
                        rep["seqs"][li] = (np.concatenate([r0[keep_], rows_s])[o_], np.concatenate([b0[keep_], blk_s])[o_])
                    else: rep["seqs"][li] = (rows_s, blk_s)
                if out["trees"] is not None:
                    T_ = out["trees"]
                    if li in rep["trees"]:                                      # merge, the cells in ascending taxon order as the single pass writes them
                        T0 = rep["trees"][li]
                        order = np.argsort(np.concatenate([T0["cells"], T_["cells"]]), kind="stable")
                        rep["trees"][li] = {k_: [ (T0[k_] + T_[k_])[o_] for o_ in order.tolist() ] for k_ in ("cells", "expected", "observed", "tips", "outgroup", "closest")}
                    else: rep["trees"][li] = T_
            PROG.finish()
        if bs_rows_rep:
            bs_rows_rep.sort()
            rl = np.array([r[0] for r in bs_rows_rep], np.int32); rt = np.array([r[1] for r in bs_rows_rep], np.int32)
            pos_ = np.searchsorted(key_a, rl.astype(np.int64) * np.int64(n) + rt.astype(np.int64))
            hit_ = (pos_ < len(key_a)) & (key_a[np.minimum(pos_, len(key_a) - 1)] == rl.astype(np.int64) * np.int64(n) + rt.astype(np.int64))
            rep["bs_locus"] = rl; rep["bs_taxon"] = rt
            rep["bs_var_log10p"] = np.where(hit_, lp1_by[np.minimum(pos_, len(lp1_by) - 1)] / LN10, np.nan).astype(np.float32)
            rep["bs_lo_log10p"] = np.array([r[2] / LN10 for r in bs_rows_rep], np.float32)
            rep["bs_hi_log10p"] = np.array([r[3] / LN10 for r in bs_rows_rep], np.float32)
            rep["bs_partner"] = np.array([r[4] for r in bs_rows_rep], np.int32)
            rep["bs_call"] = np.array([r[5] for r in bs_rows_rep], np.int8)
        bs_summary = dict(tested=int(len(x_a)), bar_passed=int(past.sum()), scale_sd=round(float(bs_tau), 4),
                          per_cell=int(bs.k), pool=int(bs.Pn), family=BS_tests,
                          close=int(tally_bs["close"]), divergent=int(tally_bs["divergent"]),
                          noted=int(tally_bs["noted"]), divergent_labels=int(tally_bs["labels"]))
        LOG(f"broad sample: {len(x_a):,} sequences tested (scale SD {bs_tau:.3f}), {int(past.sum()):,} past the bar; "
            f"new masks: {tally_bs['close']} too close + {tally_bs['divergent']} too divergent; "
            f"{tally_bs['noted']} fired on cells already masked, {tally_bs['labels']} divergent cells labelled also-too-close")
    rep["bs_summary"] = bs_summary

    # the report's remaining pieces: sparse rows stacked, the closest taxa by m, the test's constants
    if not dense:
        rep["sparse_locus"] = np.concatenate(sparse_loc) if sparse_loc else np.zeros(0, np.int32)
        rep["sparse_taxon"] = np.concatenate(sparse_tax) if sparse_tax else np.zeros(0, np.int32)
        for k, dt, _ in NUMCOLS: rep[k] = np.concatenate(rep[k]).astype(dt) if rep[k] else np.zeros(0, dt)
        rep["comparators"] = np.concatenate(rep["comparators"]) if rep["comparators"] else np.zeros((0, K_eff), np.int32)
        for k in ("comp_div_log10p", "comp_sim_log10p"): rep[k] = np.concatenate(rep[k]) if rep[k] else np.zeros((0, K_eff), np.float32)
    cl = np.full((n, min(5, n)), -1, np.int32); cl[:, :min(5, lists["cand"].shape[1])] = lists["cand"][:, :5]      # the closest non-identical relatives
    for g, (cand, _, _, _) in P.ext.items(): cl[g] = -1; cl[g, :min(5, len(cand))] = cand[:5]
    rep.update(taxa_closest=cl, rate=rate, plateaus=plateaus, kq=kq, logbar=logbar, phi_bar=phi_bar,
               tests=dict(close_relative=D_tests, closest_match_relative=D_tests, closest_match_outside=S_tests, broad_sample=BS_tests, total=N_tests, locus_sets=nfork))
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
