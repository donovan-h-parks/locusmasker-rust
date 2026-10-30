"""Alignments at the edge: so much missing data that nothing can be tested, every sequence identical, almost
identical, a single contaminant among identical sequences. The tool must behave sensibly and say what it did."""
import re
import numpy as np
import pytest
from tests.conftest import lo, run_tool

pytestmark = pytest.mark.e2e
DNA = np.frombuffer(b"ACGT", dtype=np.uint8)


def write(d, X, nloc, name="a"):
    """X: taxa x sites of uint8 letters (ord('-') for missing); loci of equal length"""
    n, L = X.shape; ll = L // nloc; aln = d / f"{name}.fasta"; parts = d / f"{name}.parts.txt"
    with open(aln, "w") as fh:
        for t in range(n): fh.write(f">t{t:03d}\n{bytes(X[t]).decode()}\n")
    parts.write_text("".join(f"DNA, locus{i + 1} = {i * ll + 1}-{(i + 1) * ll}\n" for i in range(nloc)))
    return str(aln), str(parts)


def simulated(rng, n, nloc, ll, d=0.1):
    """a star-like alignment: every taxon the ancestor with a fraction d of its sites changed"""
    anc = rng.integers(0, 4, nloc * ll); X = np.tile(anc, (n, 1))
    mut = rng.random(X.shape) < d; X[mut] = (X[mut] + rng.integers(1, 4, int(mut.sum()))) % 4
    return DNA[X]


def not_tested_summary(log):
    """the numbers the log must state: taxa, loci and sequences not tested"""
    m = re.search(r"not tested: (\d+) of (\d+) taxa[^\n]*?(\d+) of (\d+) loci[^\n]*?(\d+) of ([\d,]+) (?:present )?sequences", log, re.I)
    assert m, "no not-tested summary in the log:\n" + log
    return tuple(int(x.replace(",", "")) for x in m.groups())


class TestNothingWorks:
    def test_everything_missing_is_a_clear_message(self, tmp_path, tmp_out):
        """100 taxa, 100 loci, 95% of the cells empty: no taxon has enough loci, nothing can be tested"""
        rng = np.random.default_rng(1); X = simulated(rng, 100, 100, 40)
        for t in range(100):
            for li in range(100):
                if rng.random() < 0.95: X[t, li * 40:(li + 1) * 40] = ord("-")
        aln, parts = write(tmp_path, X, 100); r = run_tool(aln, tmp_out, "dna", parts)
        assert "Traceback" not in r.log and r.rc != 0, r.log
        assert re.search(r"100 taxa[^\n]*excluded|excluded[^\n]*100 taxa|every taxon", r.log, re.I) and re.search(r"nothing[^\n]*tested|no sequence[^\n]*tested", r.log, re.I), r.log

    def test_no_pair_overlaps_enough_loci(self, tmp_path, tmp_out):
        """every taxon has data at 12 loci, but hardly any pair shares 5: no comparators anywhere"""
        rng = np.random.default_rng(2); X = simulated(rng, 100, 100, 40)
        for t in range(100):
            keep = set(rng.choice(100, 12, replace=False).tolist())
            for li in range(100):
                if li not in keep: X[t, li * 40:(li + 1) * 40] = ord("-")
        aln, parts = write(tmp_path, X, 100); r = run_tool(aln, tmp_out, "dna", parts)
        assert "Traceback" not in r.log, r.log
        assert r.rc != 0 and re.search(r"nothing[^\n]*tested|no sequence[^\n]*tested", r.log, re.I), r.log
        assert re.search(r"overlap in <5 loci|fewer than 5 loci|never compared", r.log, re.I), r.log

    def test_half_missing_still_works(self, tmp_path, tmp_out):
        rng = np.random.default_rng(3); X = simulated(rng, 100, 100, 40, d=0.08)
        for t in range(100):
            for li in range(100):
                if rng.random() < 0.5: X[t, li * 40:(li + 1) * 40] = ord("-")
        aln, parts = write(tmp_path, X, 100); r = run_tool(aln, tmp_out, "dna", parts)
        assert r.rc == 0 and "Traceback" not in r.log; c = r.report.man["counts"]
        assert c["cells"] > 0.8 * c["present"] and c["masked"] < 0.05 * c["cells"]
        nt = not_tested_summary(r.log); assert nt[0] == c["excluded_taxa"] and nt[2] == c["skipped_loci"] and nt[4] == c["present"] - c["cells"] and nt[5] == c["present"]

    def test_the_not_tested_summary_matches_the_report(self, runs):
        run = runs("dna_small"); c = run.report.man["counts"]; nt = not_tested_summary(run.log)
        assert nt == (c["excluded_taxa"], c["taxa"], c["skipped_loci"], c["loci"], c["present"] - c["cells"], c["present"])
        assert "untested_sequences" in c and c["untested_sequences"] == c["present"] - c["cells"]


class TestIdentical:
    def test_all_identical_sequences_cannot_be_tested_and_say_so(self, tmp_path, tmp_out):
        """with no pair resolved beyond its resolution floor no locus has a rate, so nothing can be tested: the run must say
        exactly that, not run over nothing"""
        rng = np.random.default_rng(4); X = simulated(rng, 20, 12, 100, d=0.0); aln, parts = write(tmp_path, X, 12)
        r = run_tool(aln, tmp_out, "dna", parts); assert "Traceback" not in r.log, r.log
        assert r.rc != 0 and re.search(r"nothing could be tested", r.log, re.I) and re.search(r"resolved beyond|too similar", r.log, re.I), r.log
        m = lo.divergence_matrix(aln, parts, "dna")[1]; off = ~np.eye(20, dtype=bool); assert np.isnan(m[off]).all()     # no locus has a rate, so no pair has a divergence

    def test_almost_identical_sequences_cannot_be_tested_and_say_so(self, tmp_path, tmp_out):
        rng = np.random.default_rng(5); X = simulated(rng, 20, 12, 100, d=0.004); aln, parts = write(tmp_path, X, 12)
        r = run_tool(aln, tmp_out, "dna", parts); assert "Traceback" not in r.log, r.log
        assert r.rc != 0 and re.search(r"nothing could be tested", r.log, re.I), r.log

    def test_slightly_divergent_sequences_are_tested_and_mask_nothing(self, tmp_path, tmp_out):
        rng = np.random.default_rng(10); X = simulated(rng, 20, 12, 100, d=0.03); aln, parts = write(tmp_path, X, 12)
        r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0 and "Traceback" not in r.log, r.log
        c = r.report.man["counts"]; assert c["masked"] == 0 and c["cells"] == 20 * 12

    @pytest.mark.xfail(strict=True, reason="LIMIT: a contaminant among identical relatives cannot be found: no locus has the 8 resolved comparator pairs a rate needs")
    def test_one_contaminant_among_identical_sequences_is_found(self, tmp_path, tmp_out):
        rng = np.random.default_rng(6); X = simulated(rng, 20, 12, 100, d=0.0)
        X[7, 300:400] = DNA[rng.integers(0, 4, 100)]                                   # taxon 7 at locus 4
        aln, parts = write(tmp_path, X, 12); r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0
        assert [(x["taxon"], x["locus"], x["reason"]) for x in r.masked] == [("t007", "locus4", "divergent")]

    def test_one_contaminant_among_slightly_divergent_sequences_is_found(self, tmp_path, tmp_out):
        rng = np.random.default_rng(11); X = simulated(rng, 20, 12, 100, d=0.03)
        X[7, 300:400] = DNA[rng.integers(0, 4, 100)]                                   # taxon 7 at locus 4
        aln, parts = write(tmp_path, X, 12); r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0
        assert [(x["taxon"], x["locus"], x["reason"]) for x in r.masked] == [("t007", "locus4", "divergent")]

    def test_a_blank_taxon_is_excluded_not_crashed_on(self, tmp_path, tmp_out):
        rng = np.random.default_rng(7); X = simulated(rng, 20, 12, 100); X[3, :] = ord("-"); aln, parts = write(tmp_path, X, 12)
        r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0 and "t003" in r.log and r.report.taxa["status"][3].startswith("excluded")

    def test_a_blank_locus_is_skipped_and_named(self, tmp_path, tmp_out):
        rng = np.random.default_rng(8); X = simulated(rng, 20, 12, 100); X[:, 500:600] = ord("-"); aln, parts = write(tmp_path, X, 12)
        r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0 and r.report.loci["status"][5] != "tested" and "locus6" in r.log

    def test_a_taxon_with_a_residue_or_two_per_locus_is_excluded(self, tmp_path, tmp_out):
        rng = np.random.default_rng(9); X = simulated(rng, 20, 12, 100)
        for li in range(12): X[4, li * 100 + 2:(li + 1) * 100] = ord("-")
        aln, parts = write(tmp_path, X, 12); r = run_tool(aln, tmp_out, "dna", parts); assert r.rc == 0 and r.report.taxa["status"][4].startswith("excluded")


class TestRandom:
    """random sequences: no structure, every pair at the plateau; nothing should be masked"""
    def random_alignment(self, d, dtype, seed, n=30, nloc=12, ll=150):
        rng = np.random.default_rng(seed); A = "ACGT" if dtype == "dna" else "ACDEFGHIKLMNPQRSTVWY"
        with open(d / "r.fasta", "w") as fh:
            for t in range(n): fh.write(f">t{t:02d}\n" + "".join(rng.choice(list(A), nloc * ll)) + "\n")
        (d / "p.txt").write_text("".join(f"l{i} = {i * ll + 1}-{(i + 1) * ll}\n" for i in range(nloc))); return str(d / "r.fasta"), str(d / "p.txt")

    @pytest.mark.parametrize("dtype", ["dna", "aa"])
    def test_random_sequences_run_and_test_every_cell(self, tmp_path, tmp_out, dtype):
        aln, parts = self.random_alignment(tmp_path, dtype, 1); r = run_tool(aln, tmp_out, dtype, parts)
        assert r.rc == 0 and "Traceback" not in r.log, r.log; c = r.report.man["counts"]
        assert c["cells"] == 30 * 12 and c["masked_divergent"] == 0 and c["excluded_taxa"] == 0 and c["skipped_loci"] == 0

    @pytest.mark.xfail(strict=True, reason="BUG: random sequences give false too-similar pairs (3 of 6 alignments): the closest match is chosen as the minimum mismatch over every taxon, which the closest-match arm's expectation ignores")
    def test_random_sequences_mask_nothing(self, tmp_path):
        found = []
        for dtype in ("dna", "aa"):
            for seed in (1, 2, 3):
                d = tmp_path / f"{dtype}{seed}"; d.mkdir(); aln, parts = self.random_alignment(d, dtype, seed); r = run_tool(aln, str(d / "out"), dtype, parts); assert r.rc == 0
                found += [(dtype, seed, x["taxon"], x["locus"], x["reason"], x["similar_role"]) for x in r.masked]
        assert found == [], found
