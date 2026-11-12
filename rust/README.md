# The Rust port of LocusMasker

`rust/` holds a port of `locusmasker.py` to Rust. It takes the same options and writes the same outputs; the report
viewer (`viewer/`: HTML, CSS, JavaScript, fonts, icons) is not ported: its files are embedded in the binary unchanged
and inlined into `report.html` exactly as the Python tool inlines them.

## Build and run

    cd rust && cargo build --release
    ./target/release/locusmasker --aln alignment.fasta --partitions partitions.txt --datatype dna --out scan --threads 0

The binary stands alone (it does not need the repository beside it). `cargo test --release` runs the unit tests.

## What "the same outputs" means

For the same input and options the port writes, byte for byte, what `locusmasker.py` writes:

- the masked alignment and `masked.csv`, as files;
- `tests.tsv.gz` (`--all-tests`), decompressed;
- every member of `report.npz` (the arrays, `taxa.json`, `loci.json`, the trees and pair records) and `manifest.json`
  apart from the fields listed below;
- every member of the report's core (inline in `report.html`) and of every `report_files/` chunk, decoded, and the page
  around the core.

Not the same, by nature: the compressed containers themselves (`report.npz`, the base64 ZIPs, the gzip stream: their
bytes hold timestamps and depend on the zlib build, as they do between two runs of the Python tool), `run.log` and the
progress lines (timings, and each tool describes its own machinery), and these manifest fields: `date`, `timings`,
`command`, `system`, `commit`, `python` and `numpy` (the port writes the Rust compiler there, so the viewer's run line
stays as it is). And the tool's version: the port is 0.15.0 and writes that where the script writes its 0.14.0 (the
manifest's `version`, `report_files/index.js` and the report's loading line); the method and every result are those of
locusmasker.py 0.14.0.

`validation/compare_ports.py` runs both tools on the fixture datasets (every input format, the option combinations, and
the size thresholds forced low so that the sparse report, the no-full-matrix mode and the extended candidate lists run)
and compares all of the above; `validation/fuzz_ports.py` does the same on random alignments with random options.

    python validation/compare_ports.py            # 42 cases
    python validation/fuzz_ports.py --cases 150 --seed 3

The test suite runs against the binary with `LOCUSMASKER_BIN=rust/target/release/locusmasker pytest tests/e2e` (see
`tests/README.md`).

## How the numbers are held to the Python

The port mirrors the order of floating-point operations of the numpy expressions it replaces: numpy's pairwise
summation, its medians, percentiles and interpolation, float32 where the Python computes in float32, Python's
`round()`, `repr()` of floats, `json.dumps` and `csv.writer`, the compensated `sum()` of Python 3.12 and later (the
neighbour-joining trees use it), numpy's PCG64 generator and `choice` (the report's fixed sample of cells), and the
tool's own Lanczos log-gamma and exact tails. `exp`, `log`, `pow` and `erfc` come from the C library, as they do in
Python and in numpy on a machine without AVX-512: so the two tools agree to the bit where they share a C library, and
a different libm (or numpy's AVX-512 loops) can move a last digit in either of them. The pinned hashes in
`tests/data/golden` show this for the Python tool itself: `loci.json` (the dispersions at full precision) and a few
trees hash differently on a machine other than the one that pinned them.

## Differences in behaviour

- **The pair sweep.** The Python walks pairs of taxon blocks with matrix products and merges block results (its
  reductions are order-independent); the port walks the rows of the pair matrix against every column, with the
  arithmetic of the Python's compiled kernel, counting shared and identical residues on bit planes with popcounts
  (AVX2 where the machine has it). The rows go in chunks of about half a megabyte of planes, the columns in the outer
  loop, so a column is read from memory once per chunk: the sweep is bound by memory traffic otherwise. Per row it
  keeps only the candidate list so far and the closest match at each locus so far, never a row of counts.
  `--block-taxa` and `--no-numba` are accepted and do nothing. Threads replace the worker processes.
- **The broad-sample test** is computed row by row in three passes over a locus (the draws and the curve; the
  too-close side; the too-divergent side and the direction) instead of on rows x members matrices: the same numbers
  from a fraction of the memory.
- **NEXUS files with comments or quoted names.** The Python hands these to Biopython; the port has its own reader for
  them (comments, quoted names, interleaved matrices, TAXA / DATA / CHARACTERS blocks, a SETS block of charset, taxset
  and charpartition statements; other blocks are ignored) and writes such a file back the way Biopython does. Its error
  messages are its own.
- **Restart points** are written and used as in the Python (`restart/rates`, `restart/pairs`, `--restart`,
  `LOCUSMASKER_STOP_AFTER`), but hold the port's own structures: a point written by one tool is not used by the other.
- **`--memory`** sets the budget the run's estimate is checked against before it starts; the port has no chunk sizes
  to derive from it.
- **Test hooks.** The constants the test suite changes in the Python module are read from the environment:
  `LOCUSMASKER_DENSE_CELLS_MAX`, `LOCUSMASKER_BEST_MATCH_ALL_MAX_TAXA`, `LOCUSMASKER_FULL_M_MAX_TAXA`,
  `LOCUSMASKER_SWEEP_BLOCK_MIN`, `LOCUSMASKER_REPORT_SAMPLE_CELLS`, `LOCUSMASKER_REPORT_SEQ_ALL_MAX`.
  `LOCUSMASKER_NO_AVX2` turns the AVX2 counting kernel off (the same counts from the scalar one).

## Layout

    src/main.rs      the command line and the run
    src/input.rs     FASTA, PHYLIP, NEXUS, partitions, the character checks, the masked alignment writers
    src/prepare.rs   residue codes, plateaus, the two pair sweeps, the locus rates, the broad-sample pools
    src/runz.rs      comparators, dispersion, the tests, the decisions, masked.csv, tests.tsv.gz
    src/report.rs    report.npz, report.html, report_files/
    src/tree.rs      neighbour joining and newick
    src/npmath.rs    numpy's reductions, the special functions and exact tails, the hashes
    src/pyfmt.rs     Python's float formatting, round(), json.dumps, csv.writer
    src/npy.rs       .npy members, ZIP, base64          src/rng.rs      numpy's default_rng(0).choice
    src/rep.rs       the report's data                  src/restart.rs  the restart points
    src/log.rs       the log and the progress lines     src/tests.rs    unit tests
