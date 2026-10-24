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
  --alpha 0.05  --min-ratio 2  --mask-similar-only  --mask-different-only  --threads 1
  (constants in the code: HALF = 0.5 mismatch for the resolution floor and the rate
  pseudocount; dispersion from pairs with >= 10 expected mismatches, shrinkage weight 30)
"""
import argparse, collections, csv, math, os, re, sys, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):   # BLAS stays single-threaded (its matrices are
    os.environ.setdefault(_v, "1")                                          # small); --threads parallelises over loci and pairs
import numpy as np
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
LR_BAR = 10.0           # |2 ln LR| above which a similar event's direction is attributed
TEST_THREAD_PAIRS = 4000  # --threads applies to the test stage only when a locus has this many comparisons on average


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


class _Unsupported(Exception):
    """the fast NEXUS reader met something it does not handle; Bio.Nexus takes over"""


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
    data = open(path, "rb").read()
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
        e = _RE_END.search(data, body_start)
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
    if len(names) != ntax or any(len(rows[t]) != nchar for t in names): raise _Unsupported("matrix dimensions")
    del data
    seen = np.zeros(256, bool)
    for t in names: seen |= np.bincount(np.frombuffer(rows[t], dtype=np.uint8), minlength=256) > 0
    valid = _NEX_VALID[dt]; valid = valid + valid.lower() + (missing or "") + (gap or "")
    if any(chr(c) not in valid for c in np.flatnonzero(seen)): raise _Unsupported("characters outside the datatype")
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


def read_input(path, partitions):
    """returns names, rows (one bytearray per taxon), loci, format, datatype declared by a NEXUS file (or None), and
    the NEXUS record the writer needs (a Bio.Nexus object, or NexusInfo from the fast reader)"""
    head = open(path).read(64).lstrip()[:7].upper()
    nexus = None; declared = None
    if head.startswith("#NEXUS"):
        try:
            names, rows, charsets, nexus = read_nexus_fast(path); dt = nexus.datatype
        except _Unsupported:
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
            names = list(nexus.taxlabels); rows = [bytearray(str(nexus.matrix[t]).encode()) for t in names]
            charsets = nexus.charsets; dt = (nexus.datatype or "").lower()
        fmt = "nexus"
        if charsets and partitions:
            sys.exit("the NEXUS file has charsets and a partition file was given: use one or the other, not both")
        if charsets:
            loci = [(k, np.array(sorted(v), dtype=int)) for k, v in charsets.items()]
        elif partitions:
            loci = parse_raxml_partitions(partitions, len(rows[0]))
        else:
            sys.exit("the NEXUS file has no charsets: add a SETS block or give --partitions")
        if dt in ("dna", "rna", "nucleotide"): declared = "dna"
        elif dt == "protein": declared = "aa"
    else:
        if not partitions: sys.exit("a FASTA or PHYLIP alignment needs --partitions")
        txt = open(path).read()
        names, seqs = read_fasta(txt) if head.startswith(">") else read_phylip(txt)
        fmt = "fasta" if head.startswith(">") else "phylip"
        rows = [bytearray(s.encode()) for s in seqs]; del seqs, txt
        loci = parse_raxml_partitions(partitions, len(rows[0]))
    L = len(rows[0])
    bad = [t for t, r in zip(names, rows) if len(r) != L]
    if bad: sys.exit(f"sequences differ in length (e.g. {bad[0]}): the input must be an alignment")
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
        from Bio.Seq import Seq
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
def _block_terms(k0, cnt, n, d, a_, b_, cpair, phi, B):
    """log pmf of the mismatch counts k0[i] .. k0[i] + cnt[i] - 1 for each pair i, as a (pairs, B) array with -inf
    beyond cnt. Binomial coefficients from the log-factorial table; the beta-binomial's two gamma terms from one
    log-gamma per pair at the block's first term and the recurrence lgamma(x + 1) = lgamma(x) + log x along the
    block, which costs a log per term instead of a log-gamma."""
    nint = n.astype(np.int64); lf = log_fact(int(n.max()) + 1)
    j = np.arange(B); kk = k0[:, None] + j[None, :]; valid = j[None, :] < cnt[:, None]
    ksafe = np.where(valid, kk, 0); nk = np.where(valid, nint[:, None] - kk, 0)
    lchoose = lf[nint][:, None] - lf[ksafe] - lf[nk]
    kf = ksafe.astype(float); nf_ = n[:, None]
    lbin = lchoose + kf * np.log(d)[:, None] + (nf_ - kf) * np.log(1 - d)[:, None]
    if phi < 1e-8: lp = lbin
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
        lbb = lchoose + lg_a + lg_b + cpair[:, None]
        lp = np.where((n >= 2)[:, None], lbb, lbin)
    return np.where(valid, lp, -np.inf)


def _logsumexp_rows(lp):
    mx = np.max(lp, axis=1); mx = np.where(np.isfinite(mx), mx, 0.0)
    return mx + np.log(np.sum(np.exp(lp - mx[:, None]), axis=1))


def bb_tails(o, n, E, phi, upper):
    """log tail probability of the mismatch count among n shared sites: binomial with mean E, or, when the locus
    is overdispersed, beta-binomial with mean E and variance E + phi E^2. upper -> log P(X >= o), else
    log P(X <= o). All pairs at once, in blocks of 64 terms per pair. The lower tail is the exact sum of its
    terms; the upper tail is summed directly from the observed count upward until a block adds less than 1e-20
    of the running sum while the pmf is falling (or n is reached), so far tails come out to full precision
    instead of as 1 - P(X < o). Both are floored at log(1e-300) and capped at log(1)."""
    o = np.asarray(o, float); n = np.asarray(n, float); E = np.asarray(E, float); m = len(o)
    if m == 0: return np.empty(0)
    E = np.minimum(np.maximum(E, 1e-9), n * (1 - 1e-9)); d = E / n
    a_ = b_ = cpair = None
    if phi >= 1e-8:
        with np.errstate(divide="ignore", invalid="ignore"):
            rho = np.where(n >= 2, np.minimum((d + phi * n * d) / ((1 - d) * np.maximum(n - 1, 1)), 1 - 1e-6), 0.5)
        a_ = d * (1 / rho - 1); b_ = (1 - d) * (1 / rho - 1)
        cpair = -lgamma_v(n + a_ + b_) - (lgamma_v(a_) + lgamma_v(b_) - lgamma_v(a_ + b_))
    oi = o.astype(np.int64); nint = n.astype(np.int64); B = 64
    out = np.full(m, -np.inf)
    if upper: active = np.flatnonzero(oi <= nint); kstart = oi.copy(); last = nint
    else: active = np.arange(m); kstart = np.zeros(m, np.int64); last = oi
    while len(active):
        cnt = np.minimum(B, last[active] - kstart[active] + 1); has = cnt > 0
        active, cnt = active[has], cnt[has]
        if not len(active): break
        lp = _block_terms(kstart[active], cnt, n[active], d[active], None if a_ is None else a_[active],
                          None if b_ is None else b_[active], None if cpair is None else cpair[active], phi, B)
        blk = _logsumexp_rows(lp)
        out[active] = np.logaddexp(out[active], blk)
        kstart[active] += cnt
        done = kstart[active] > last[active]
        if upper:
            rows = np.arange(len(active)); falling = lp[rows, cnt - 1] < lp[rows, 0]
            done |= falling & (blk < out[active] - 46.0)
        active = active[~done]
    return np.maximum(np.minimum(out, 0.0), _LOG_FLOOR)


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


def pair_divergences(S, W, n, iu, rr=None, ok_r=None, chunk_bytes=256 << 20, threads=1):
    """for every taxon pair (flattened n*n): the site-weighted median over loci of the pair's divergence
    (rate-normalised by rr when given), the same median of the resolution floor (half a mismatch per locus),
    and the number of loci contributing; pairs with fewer than MIN_LOCI_FOR_M loci get NaN. Self pairs carry
    no vote; loci without a rate (ok_r False) carry none either. Computed in chunks of pairs straight from
    the per-locus counts, so the peak memory is a few chunk-sized temporaries rather than several copies of
    the whole loci x pairs matrix; the per-pair arithmetic is exactly that of pair_median."""
    nloc, npair = S.shape                                               # the pairs are the upper triangle, diagonal included
    m = np.full(npair, np.nan); m_fl = np.full(npair, np.nan); nl = np.zeros(npair, dtype=np.int64)
    diag = np.flatnonzero(iu[0] == iu[1])
    chunk = max(int(chunk_bytes // (8 * 8 * max(nloc, 1) * max(threads, 1))), 1)   # the budget is shared by the threads' chunks
    if threads > 1: chunk = max(min(chunk, -(-npair // (4 * threads))), 1)      # enough chunks to keep the workers busy
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
    starts = range(0, npair, chunk)
    if threads > 1:
        with ThreadPoolExecutor(threads) as ex: list(ex.map(do_chunk, starts))
    else:
        for c0 in starts: do_chunk(c0)
    return _expand(m, n, iu), _expand(m_fl, n, iu), _expand(nl, n, iu)


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
    ap.add_argument("--threads", type=int, default=1, help="worker threads for the per-locus counts and the pair divergences")
    a = ap.parse_args(); K = a.n_comp; a.threads = max(1, a.threads)
    if K < 3 or K % 2 == 0: sys.exit("--n-comp must be an odd number of at least 3")
    if a.mask_similar_only and a.mask_different_only:
        sys.exit("--mask-similar-only and --mask-different-only are mutually exclusive")
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    names, rows, loci, fmt, declared, nexus = read_input(a.aln, a.partitions)
    if declared and declared != a.datatype:
        sys.exit(f"--datatype {a.datatype} was given but the NEXUS file declares {declared}")
    dtype = a.datatype
    if len(loci) < 10:
        sys.exit(f"only {len(loci)} loci defined: this is a phylogenomic test and needs at least 10")
    n, L = len(names), len(rows[0]); alphabet = DNA if dtype == "dna" else PROTEIN
    if n < 5: sys.exit(f"only {n} taxa: the test needs at least 5")
    K_eff = min(K, n - 1); K_eff -= 1 - K_eff % 2                    # largest odd <= min(K, n-1)
    M_eff = (K_eff + 1) // 2                                          # fails that make a majority
    if K_eff < K: print(f"n-comp shrunk to {K_eff} ({n} taxa)", file=sys.stderr)
    print(f"{n} taxa, {L} sites, {len(loci)} loci, {dtype}, {fmt}", file=sys.stderr)
    open(f"{a.out}/taxa.txt", "w").write("\n".join(names) + "\n")

    lut = np.full(256, -1, dtype=np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = alphabet.index("T")
    X = np.empty((n, L), dtype=np.int8)                               # one row at a time: no joined copy of the alignment
    for i, r in enumerate(rows): X[i] = lut[np.frombuffer(r, dtype=np.uint8)]

    # ---- pass 1: per-locus shared/same counts; alignment-wide identity falls out of the sums
    nloc = len(loci)
    # the same-residue and shared-site counts of every pair, per locus: symmetric, so the upper triangle (diagonal
    # included, row-major) is kept and a locus's full matrix is rebuilt when the test needs it
    iu = np.triu_indices(n); ntri = len(iu[0]); diag_pos = np.flatnonzero(iu[0] == iu[1])
    tri_idx = np.empty((n, n), np.int64); tri_idx[iu] = np.arange(ntri); tri_idx[iu[1], iu[0]] = np.arange(ntri)
    cdt = np.uint16 if max(len(sites) for _, sites in loci) < 65536 else np.uint32   # counts never exceed a locus's length
    S = np.empty((nloc, ntri), cdt); W = np.empty((nloc, ntri), cdt)
    plateaus = np.empty(nloc)
    def locus_counts(li):
        sites = loci[li][1]
        Xm = X[:, sites]; Pm = (Xm >= 0).astype(np.float32)
        sh = Pm @ Pm.T; sm = np.zeros((n, n), np.float32)
        for i in range(len(alphabet)):
            Ai = (Xm == i).astype(np.float32); sm += Ai @ Ai.T
        S[li] = sm[iu]; W[li] = sh[iu]                                  # exact: float32 counts below 2^24
        cnt = np.bincount(Xm[Xm >= 0], minlength=len(alphabet)).astype(float)
        f = cnt / max(cnt.sum(), 1.0)
        plateaus[li] = max(100.0 * (1.0 - float((f ** 2).sum())), 1.0)
    if a.threads > 1:
        with ThreadPoolExecutor(a.threads) as ex: list(ex.map(locus_counts, range(nloc), chunksize=max(1, nloc // (8 * a.threads))))
    else:
        for li in range(nloc): locus_counts(li)
    present = (W[:, diag_pos] >= a.fit_overlap).sum(0)
    thin = np.flatnonzero(present < MIN_LOCI_PER_TAXON)
    open(f"{a.out}/excluded_taxa.txt", "w").write("".join(f"{names[i]}\t{int(present[i])}\n" for i in thin))
    if len(thin):
        drop = np.isin(iu[0], thin) | np.isin(iu[1], thin)
        S[:, drop] = 0; W[:, drop] = 0                                     # out of the test entirely
        print(f"{len(thin)} taxa with fewer than {MIN_LOCI_PER_TAXON} loci present are excluded from the test "
              f"(listed in excluded_taxa.txt, left unmasked)", file=sys.stderr)
    tot_same = _expand(S.sum(0, dtype=np.float64), n, iu); tot_shared = _expand(W.sum(0, dtype=np.float64), n, iu)
    with np.errstate(invalid="ignore", divide="ignore"):
        ani = np.where(tot_shared > 0, 100.0 * tot_same / tot_shared, np.nan)
    np.fill_diagonal(ani, np.nan); np.save(f"{a.out}/ani.npy", ani)
    print(f"per-locus identities done, {time.time() - t0:.0f}s", file=sys.stderr)

    # ---- pair divergence m: capped site-weighted median across loci of rate-normalised divergence
    # ONE rate per locus, estimated once and used everywhere (to normalise loci inside m and as the
    # expectation in the test). It is estimated where the test looks: each taxon's nearest neighbours
    # by the provisional pair divergence m0 (the raw site-weighted median over loci, available before rates exist), one equal
    # vote per neighbour pair resolved beyond its floor, of the multiple-hit-corrected ratio
    # c(d + half a mismatch) / c(m0), m0 the provisional pair divergence (raw site-weighted median).
    m0, m0_fl, _ = pair_divergences(S, W, n, iu, threads=a.threads)   # the floor: the imputed half mismatch per locus, as a percentage
    m0_hi = np.maximum(m0, m0_fl); res0 = m0 > m0_fl
    dist = np.where(np.isfinite(m0), m0, np.inf); np.fill_diagonal(dist, np.inf)   # neighbours by the provisional pair divergence (a median over loci, not swayed by one long fast locus)
    Kn = min(K_eff, n - 1)
    nb = np.argsort(dist, axis=1)[:, :Kn]
    gi0 = np.repeat(np.arange(n), Kn); ci0 = nb.ravel()
    okn = np.isfinite(dist[gi0, ci0]) & res0[gi0, ci0]
    gi0, ci0 = gi0[okn], ci0[okn]
    m0n = m0_hi[gi0, ci0]
    rate = np.full(nloc, np.nan)
    t0i = tri_idx[gi0, ci0]
    for li in range(nloc):
        o_, s_ = W[li][t0i].astype(np.float32), S[li][t0i].astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            d_ = np.where(o_ > 0, 100.0 * (1.0 - s_ / np.maximum(o_, 1)), np.nan)
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
    m, m_fl, nl = pair_divergences(S, W, n, iu, rr, ok_r, threads=a.threads)
    np.fill_diagonal(m, np.nan); np.save(f"{a.out}/m.npy", m)
    # resolution floor of the pair: the median over its loci of half a mismatch per locus (rate-normalised);
    # the divergence used on the too-divergent side is the larger of the measured median and this floor
    m_hi = np.maximum(m, m_fl); np.fill_diagonal(m_hi, np.nan); np.save(f"{a.out}/m_hi.npy", m_hi)
    m_xs = np.maximum(m - m_fl, 0.0)                                  # divergence resolved beyond the floor: the rate-vote weight
    keep = np.ones(n, bool); keep[thin] = False
    nlk = nl[np.ix_(keep, keep)][np.triu_indices(int(keep.sum()), 1)]
    fb = float((nlk < MIN_LOCI_FOR_M).mean()) if len(nlk) else 0.0
    print(f"pair divergences done, {time.time() - t0:.0f}s"
          + (f"; note: {fb:.0%} of taxon pairs overlap in <{MIN_LOCI_FOR_M} loci and are never compared"
             if fb > 0.05 else ""), file=sys.stderr)

    # deterministic tie-break for comparator ranking: pooled divergence, then taxon order
    tie = np.nan_to_num(100.0 - ani, nan=100.0) * 1e-9 + np.arange(n)[None, :] * 1e-12
    masked = rows                                                      # blanked in place; X is its own copy

    nmasked = run_z(a, names, loci, X, S, W, iu, diag_pos, m, m_hi, m_xs, plateaus, rate, tie, masked, K_eff, M_eff)

    out = write_masked(a.aln, a.out, fmt, names, masked, nexus)
    print(f"masked {nmasked} cells -> {out}\ndone {time.time() - t0:.0f}s", file=sys.stderr)


# ------------------------------------------------------------------ the counting test
def run_z(a, names, loci, X, S, W, iu, diag_pos, m, m_hi, m_xs, plateaus, rate, tie, masked, K_eff, M_eff):
    n = len(names); nloc = len(loci)
    mj = np.where(np.isfinite(m), m, np.inf) + tie
    mask_hi = not a.mask_similar_only
    mask_lo = not a.mask_different_only
    views = [np.frombuffer(mb, dtype=np.uint8) for mb in masked]        # writable views for blanking

    # pre-pass: count distinct comparisons so the Bonferroni bar is known before any test runs
    Npairs = 0; sel = []
    for li in range(nloc):
        ov = _expand(W[li].astype(np.float32), n, iu)
        nres = W[li][diag_pos].astype(int)
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
        ov, sm = _expand(W[li].astype(np.float32), n, iu), _expand(S[li].astype(np.float32), n, iu)
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

    # worker threads pay off in the test only when a locus's arrays are large (thousands of pairs): with small loci
    # the per-call interpreter-lock traffic of many small numpy operations makes threads slower than one
    ntested = sum(1 for t in sel if (t[2] > 0).any())
    pool = ThreadPoolExecutor(a.threads) if a.threads > 1 and Npairs / max(ntested, 1) >= TEST_THREAD_PAIRS else None
    # dispersion pre-pass: the robust spread of observed/expected among each locus's resolved comparator
    # pairs with enough expected mismatches, with the Poisson part removed: Var(O/E) = 1/E + phi
    def locus_phi(li):
        L = locus_pairs(li)
        if isinstance(L, str): return np.nan, 0.0
        E = L["o"] * L["expect"](L["mph"]) / 100.0
        use = L["fit"] & (E >= PHI_MIN_E)
        if use.sum() < 10: return np.nan, 0.0
        r = np.log((L["mis"][use] + HALF) / E[use])
        s_ = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        return max(s_ ** 2 - float(np.median(1.0 / E[use])), 0.0), float(use.sum())
    phi_hat = np.full(nloc, np.nan); phi_n = np.zeros(nloc)
    for li, (h, c) in enumerate(pool.map(locus_phi, range(nloc)) if pool else map(locus_phi, range(nloc))):
        phi_hat[li] = h; phi_n[li] = c
    okp = phi_n > 0
    phi_bar = wmedian1(phi_hat[okp], phi_n[okp]) if okp.any() else 0.0
    print(f"dispersion: dataset-wide phi {phi_bar:.4f} (robust SD of O/E beyond Poisson {math.sqrt(phi_bar):.2f}); "
          f"per-locus phi 10/50/90%: {np.nanpercentile(phi_hat[okp], [10, 50, 90]).round(4).tolist() if okp.any() else 'n/a'}", file=sys.stderr)
    logbar = math.log(a.alpha / max(2 * Npairs, 1))
    zeros = "\t".join("0" for _ in range(M_eff, K_eff + 1))

    def process_locus(li):
        """everything the test decides and writes for one locus, as text blocks and mask actions; touches no shared
        state, so loci can be processed in parallel and written in order"""
        name, sites = loci[li]
        L = locus_pairs(li)
        if isinstance(L, str):
            nres = sel[li][3]
            return dict(skipped=True, log=f"{name}: {L}",
                        per=f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t0\t0\t{rate[li]:.3f}\t{plateaus[li]:.1f}\t{kq:.2f}\t0\t{zeros}\t0\t0\t{L}\n")
        idx, gi, ci, rk, o, same, mis, ident, dob, mp, mph, fit, ov, sm, ncomp, nres, expect, corr, rate_c = \
            (L[k] for k in ("idx", "gi", "ci", "rk", "o", "same", "mis", "ident", "dob", "mp", "mph", "fit", "ov", "sm", "ncomp", "nres", "expect", "corr", "rate_c"))
        phi = (phi_n[li] * (phi_hat[li] if okp[li] else 0.0) + PHI_N0 * phi_bar) / (phi_n[li] + PHI_N0)

        # the count test with a quasi-Poisson variance: Var(O) = E + phi E^2, phi the locus's dispersion
        # (shrunk toward the dataset's); phi = 0 is the plain Poisson test
        dhat = expect(mph)                                          # divergent side: the divergence we cannot exclude
        dlo = expect(mp)                                            # similar side: the divergence we can
        E = o * dhat / 100.0
        E_lo = o * dlo / 100.0
        # exact binomial tails over the shared sites (beta-binomial when the locus is overdispersed), evaluated
        # only for the pairs that clear the effect floors; a pair fails when its tail probability is below alpha / 2N
        excess = mis - E
        fail_hi = excess > (a.min_ratio - 1) * E
        if fail_hi.any():
            fail_hi[fail_hi] = bb_tails(mis[fail_hi], o[fail_hi], E[fail_hi], phi, True) < logbar
        defi = E_lo - mis
        fail_lo = defi > (2.0 / 3.0) * E_lo
        if fail_lo.any():
            fail_lo[fail_lo] = bb_tails(mis[fail_lo], o[fail_lo], E_lo[fail_lo], phi, False) < logbar

        # pairs.tsv: one %-format of the whole locus's rows (the locus name is part of the format)
        lname = name.replace("%", "%%"); c_g = [names[i] for i in gi.tolist()]; c_c = [names[i] for i in ci.tolist()]
        pairs_block = (lname + "\t%s\t%d\t%s\t%.3f\t%d\t%.2f\t%.2f\t%.1f\t%d\t%d\n") * len(gi) % tuple(chain.from_iterable(zip(
            c_g, (rk + 1).tolist(), c_c, mph.tolist(), o.astype(np.int64).tolist(), ident.tolist(), (100 - dhat).tolist(), excess.tolist(),
            fail_hi.astype(np.int64).tolist(), fail_lo.astype(np.int64).tolist())))

        # ---- divergent side: majority over comparators, every cell at once
        kc = ncomp[idx]; ncell = len(idx)
        starts = np.r_[0, np.cumsum(kc)[:-1]]
        nf = np.add.reduceat(fail_hi.astype(np.int64), starts); nsim = np.add.reduceat(fail_lo.astype(np.int64), starts)
        flag = (kc >= K_eff) & (nf >= M_eff)
        amb = (kc < K_eff) & (nf >= kc // 2 + 1)
        fc = np.bincount(nf[flag], minlength=K_eff + 1); fcount = {j: int(fc[j]) for j in range(M_eff, K_eff + 1)}
        namb = int(amb.sum())
        # medians over each cell's comparators: the pairs laid out one row per cell, NaN beyond the cell's count
        row = np.repeat(np.arange(ncell), kc); col = np.arange(len(gi)) - np.repeat(starts, kc)
        def cell_median(v):
            P = np.full((ncell, K_eff), np.nan); P[row, col] = v; P.sort(axis=1)            # NaN sorts last
            ar = np.arange(ncell)
            return (P[ar, (kc - 1) // 2] + P[ar, kc // 2]) / 2                             # np.median's mean of the middle two
        med_i = cell_median(ident); med_e = cell_median(100 - dhat); med_x = cell_median(excess)
        # the cell's extra divergence beyond expectation, as corrected distance in units of m (a typical
        # locus's divergence): a foreign or paralogous sequence adds a constant here whatever the comparator
        exm = (corr(100 - med_i) - corr(100 - med_e)) / rate_c
        s0 = starts
        c_g = [names[i] for i in idx.tolist()]; c_c = [names[i] for i in ci[s0].tolist()]
        cells_block = (lname + "\t%s\t%d\t%s\t%.3f\t%d\t%d\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.1f\t%d\t%d\t%d\t%d\t%.3f\n") * ncell % tuple(chain.from_iterable(zip(
            c_g, kc.tolist(), c_c, mph[s0].tolist(), nres[idx].tolist(), o[s0].astype(np.int64).tolist(), ident[s0].tolist(), (100 - dhat[s0]).tolist(),
            excess[s0].tolist(), med_i.tolist(), med_e.tolist(), med_x.tolist(), nf.tolist(), nsim.tolist(), flag.astype(np.int64).tolist(),
            amb.astype(np.int64).tolist(), exm.tolist())))
        masks = []
        if mask_hi:
            for q in np.flatnonzero(flag | amb).tolist():
                masks.append((int(idx[q]), "divergent",
                              [int(nf[q]), int(nsim[q]), c_c[q], "%.3f" % mph[s0[q]], "%.2f" % med_i[q], "%.2f" % med_e[q], "%.1f" % med_x[q], "", "%.3f" % exm[q]]))

        # ---- similar side: any failing pair; LR direction attribution on shared sites
        simpairs = {}
        for j in np.flatnonzero(fail_lo):
            key = (int(min(gi[j], ci[j])), int(max(gi[j], ci[j])))
            simpairs.setdefault(key, j)
        Xs = X[:, sites] if simpairs else None
        sim_lines = []
        for (A, B) in sorted(simpairs):
            lam, donor, recip = attribute(A, B, Xs, m, ov, expect, a.fit_overlap)
            o_ab = float(ov[A, B]); id_ab = 100.0 * sm[A, B] / max(o_ab, 1.0)
            exp_ab = float(expect(m[A, B])); defi_ab = o_ab * exp_ab / 100 - (o_ab - sm[A, B])
            targets = ([recip] if recip is not None else [A, B]) if mask_lo else []
            side = ("both" if recip is None else names[recip]) if mask_lo else "none"
            for g in targets:
                masks.append((g, "similar",
                              ["", "", names[B if g == A else A], f"{m[A, B]:.3f}", f"{id_ab:.2f}", f"{100 - exp_ab:.2f}",
                               f"{-defi_ab:.1f}", f"{lam:.1f}", f"{float((corr(100 - id_ab) - corr(exp_ab)) / rate_c):.3f}"]))
            sim_lines.append(f"{name}\t{names[A]}\t{names[B]}\t{m[A, B]:.3f}\t{int(o_ab)}\t{id_ab:.2f}\t{100 - exp_ab:.2f}"
                             f"\t{defi_ab:.1f}\t{lam:.1f}\t{names[donor] if donor is not None else 'NA'}"
                             f"\t{names[recip] if recip is not None else 'NA'}\t{side}\n")

        gd = expect(grid_m)
        trends_block = "".join(f"{name}\t{gm:.4f}\t{100 - gdv:.4f}\t{100 - min(a.min_ratio * gdv, plateaus[li]):.4f}\t{100 - gdv / 3:.4f}\n"
                               for gm, gdv in zip(grid_m, gd))
        nflag = sum(fcount.values()); fstr = "/".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
        per_line = (f"{name}\t{len(sites)}\t{int((nres > 0).sum())}\t{len(idx)}\t{len(gi)}\t"
                    f"{rate[li]:.3f}\t{plateaus[li]:.1f}\t{kq:.2f}\t{nflag}\t"
                    + "\t".join(str(fcount[j]) for j in range(M_eff, K_eff + 1))
                    + f"\t{namb}\t{len(simpairs)}\ttested\n")
        disp_line = f"{name}\t{(phi_hat[li] if okp[li] else float('nan')):.4f}\t{int(phi_n[li])}\t{phi:.4f}\n"
        log = (f"{name}: {len(idx)} cells, {len(gi)} pairs, flagged {nflag} "
               f"(fails {M_eff}..{K_eff} of {K_eff}: {fstr}), ambiguous {namb}, similar {len(simpairs)}")
        return dict(skipped=False, pairs=pairs_block, cells=cells_block, masks=masks, simf="".join(sim_lines),
                    trends=trends_block, per=per_line, disp=disp_line, log=log)

    def consume(li, R):
        name, sites = loci[li]
        if R["skipped"]:
            per.write(R["per"]); print(R["log"], file=sys.stderr); return
        pairs.write(R["pairs"]); cells.write(R["cells"])
        for g, side, ev in R["masks"]: mask_cell(name, g, sites, side, ev)
        simf.write(R["simf"]); trends.write(R["trends"]); per.write(R["per"]); disp.write(R["disp"])
        print(R["log"], file=sys.stderr)

    if pool:
        batch = 8 * a.threads                                           # loci in flight, so the text blocks do not pile up
        for b0 in range(0, nloc, batch):
            for li, R in zip(range(b0, min(b0 + batch, nloc)), pool.map(process_locus, range(b0, min(b0 + batch, nloc)))): consume(li, R)
        pool.shutdown()
    else:
        for li in range(nloc): consume(li, process_locus(li))
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
