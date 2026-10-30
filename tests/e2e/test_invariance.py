"""Reproducibility: the same input gives byte-identical outputs whatever the input format, the thread count, the
chunking of the work, or the memory the machine reports. These guard the promise that a run is the same anywhere."""
import os
import numpy as np
import pytest
from tests.conftest import lo, dataset, run_tool, npz_hashes, sha

pytestmark = pytest.mark.e2e


def outputs(run):
    return open(run.file("masked.csv")).read(), npz_hashes(run.file("report.npz"))


class TestFormats:
    @pytest.mark.parametrize("fmt", ["nexus", "nexus_quirky", "phylip", "phylip_interleaved"])
    def test_every_input_format_gives_the_same_result(self, runs, tmp_out, fmt):
        ds = dataset("dna_small"); base = outputs(runs("dna_small"))
        path = ds.file(fmt); r = run_tool(path, tmp_out, "dna", None if "nexus" in fmt else ds.parts)
        assert r.rc == 0, r.log
        got = outputs(r); assert got[0] == base[0]
        diff = {k for k in base[1] if base[1][k] != got[1].get(k)}
        assert not diff, diff

    def test_nexus_output_keeps_its_charsets_and_the_masks(self, runs, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.file("nexus"), tmp_out, "dna"); out = r.file("alignment.masked.nex")
        names, rows, charsets, info = lo.read_nexus_fast(out)
        base = lo.read_fasta(open(runs("dna_small").file("alignment.masked.fasta")).read())
        assert names == base[0] and [bytes(x).decode() for x in rows] == base[1] and list(charsets) == [f"locus{i + 1}" for i in range(ds.nloci)]

    def test_quirky_nexus_output_is_written_by_biopython(self, tmp_out):
        pytest.importorskip("Bio.Nexus"); from Bio.Nexus import Nexus
        ds = dataset("dna_small"); r = run_tool(ds.file("nexus_quirky"), tmp_out, "dna"); out = r.file("alignment_quirky.masked.nex")
        nex = Nexus.Nexus(out); names, seqs = lo.read_fasta(open(ds.aln).read())
        assert list(nex.taxlabels) == names and len(nex.charsets) == ds.nloci
        for row in r.masked:
            g, (a, b) = names.index(row["taxon"]), ds.truth["locus_bounds"][int(row["locus"][5:]) - 1]
            assert str(nex.matrix[row["taxon"]])[a - 1:b] == "-" * (b - a + 1)

    def test_phylip_output_round_trips(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.file("phylip"), tmp_out, "dna", ds.parts)
        names, seqs = lo.read_phylip(open(r.file("alignment.masked.phy")).read())
        assert names == lo.read_fasta(open(ds.aln).read())[0] and all(len(s) == ds.truth["sites"] for s in seqs)


class TestThreads:
    @pytest.mark.parametrize("name", ["dna_small", "dna_dense", "aa_small"])
    @pytest.mark.parametrize("threads", ["2", "7"])
    def test_threads_do_not_change_the_outputs(self, runs, tmp_out, name, threads, monkeypatch):
        """with the thread pools forced on even for small loci"""
        monkeypatch.setattr(lo, "PASS1_FLOPS_PER_THREAD", 1.0); monkeypatch.setattr(lo, "DISPERSION_PAIRS_PER_WORKER", 1)
        monkeypatch.setattr(lo, "COMPARATOR_WORK_PER_WORKER", 1.0); monkeypatch.setattr(lo, "TEST_BATCH_COMPARISONS", 100); monkeypatch.setattr(lo, "TEST_BATCH_MIN", 100)
        ds = dataset(name); r = run_tool(ds.aln, tmp_out, ds.dtype, ds.parts, ["--threads", threads]); assert r.rc == 0
        assert outputs(r) == outputs(runs(name))

    def test_the_thread_count_is_recorded(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--threads", "3"]); assert r.report.man["options"]["threads"] == 3


class TestChunking:
    def test_pair_chunking_does_not_change_the_outputs(self, runs, tmp_out, monkeypatch):
        monkeypatch.setattr(lo, "LSET_CHUNK", 3)
        monkeypatch.setattr(lo, "_mem_available", lambda: 64 << 20)                      # tiny pair chunks, no cached comparator arrays
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--threads", "4"]); assert r.rc == 0
        assert outputs(r) == outputs(runs("dna_small"))

    def test_running_twice_is_identical(self, runs, tmp_out):
        ds = dataset("aa_small"); r = run_tool(ds.aln, tmp_out, "aa", ds.parts)
        assert outputs(r) == outputs(runs("aa_small"))
        assert sha(r.file("alignment.masked.fasta")) == sha(runs("aa_small").file("alignment.masked.fasta"))

    def test_the_sparse_report_carries_the_same_decisions(self, runs, tmp_out, monkeypatch):
        """above DENSE_CELLS_MAX cells the per-cell numbers are carried only for cells with a failing test"""
        monkeypatch.setattr(lo, "DENSE_CELLS_MAX", 10)
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); R = r.report; B = runs("dna_small").report
        assert not R.dense and R.man["carried"]["numeric_columns"] == "sparse" and open(r.file("masked.csv")).read() == open(runs("dna_small").file("masked.csv")).read()
        for k in ("n_comp", "n_fail", "n_similar", "side"): assert np.array_equal(R.arr(f"cells/{k}.npy"), B.arr(f"cells/{k}.npy")), k
        side = B.arr("cells/side.npy"); nf = B.arr("cells/n_fail.npy")
        rows = set(R.sparse); want = {(g, li) for g, li in zip(*np.nonzero((nf > 0) | ((side & (2 | 128)) > 0)))}
        assert rows == want
        for (g, li) in list(rows)[:30]:
            assert R.cell("excess_m", g, li) == B.cell("excess_m", g, li) and R.comps(g, li) == B.comps(g, li)
            assert np.array_equal(R.cell("comp_div_log10p", g, li), B.cell("comp_div_log10p", g, li), equal_nan=True)

    def test_best_matches_for_failing_cells_only_above_the_taxon_cap(self, runs, tmp_out, monkeypatch):
        monkeypatch.setattr(lo, "BEST_MATCH_ALL_MAX_TAXA", 5)
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); R = r.report; B = runs("dna_small").report
        assert R.man["carried"]["best_match_cells"] == "fails"
        bm, bb = R.arr("cells/best_match.npy"), B.arr("cells/best_match.npy"); side = B.arr("cells/side.npy"); nf = B.arr("cells/n_fail.npy"); nc = B.arr("cells/n_comp.npy")
        fail = (nf > 0) | ((side & (2 | 128)) > 0)
        assert np.array_equal(bm[fail], bb[fail]) and (bm[(nc > 0) & ~fail] == -1).all()

    def test_sequences_for_failing_cells_only_above_the_size_cap(self, tmp_out, monkeypatch):
        monkeypatch.setattr(lo, "REPORT_SEQ_ALL_MAX", 1)
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); R = r.report
        assert R.man["carried"]["sequences"] == "fails"
        side = R.arr("cells/side.npy"); nf = R.arr("cells/n_fail.npy"); pt = R.arr("cells/partner.npy")
        for li in range(R.nloc):
            want = set()
            for g in range(R.n):
                if nf[g, li] > 0 or side[g, li] & (2 | 128): want.add(g); want.update(R.comps(g, li)); want.add(int(pt[g, li])) if pt[g, li] >= 0 else None
            got = set(R.arr(f"seq/{li}_taxa.npy").tolist()) if R.has(f"seq/{li}.npy") else set()
            assert got == want, li


class TestKernel:
    """the compiled pair stage (numba) gives the outputs of the numpy one, whatever the threads and the slabs"""
    @pytest.mark.parametrize("name,threads", [("dna_small", 1), ("dna_small", 4), ("aa_random_trees", 2)])
    def test_the_compiled_pair_stage_gives_the_same_outputs(self, runs, tmp_out, name, threads, monkeypatch):
        if lo._pair_kernel(threads) is None: pytest.skip("numba is not installed")
        monkeypatch.setattr(lo, "NUMBA_PAIR_LOCI_MIN", 0); monkeypatch.setattr(lo, "NUMBA_PAIR_SLAB", 1000); monkeypatch.setattr(lo, "PAIRS_PER_THREAD", 1)
        ds = dataset(name); r = run_tool(ds.aln, tmp_out, ds.dtype, ds.parts, ["--threads", str(threads)]); assert r.rc == 0, r.log
        assert r.log.count("compiled kernel (numba)") == 2 and outputs(r) == outputs(runs(name))
