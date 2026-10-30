"""The command line: argument checks and their messages, the options' effects, the files a run leaves behind, and the
log. Most runs are in-process (coverage sees them); a few are real processes, the way a user runs the tool."""
import csv, gzip, json, os, re, zipfile
import numpy as np
import pytest
from tests.conftest import lo, dataset, run_tool, run_cli, Report


@pytest.mark.e2e
class TestArguments:
    def test_help_lists_the_options(self):
        r = run_cli(["--help"]); assert r.returncode == 0
        for opt in ("--aln", "--partitions", "--out", "--datatype", "--n-comp", "--fit-overlap", "--alpha", "--donor-threshold", "--threads", "--no-report", "--all-tests", "--name"):
            assert opt in r.stdout, opt

    def test_required_arguments(self, tmp_out):
        ds = dataset("dna_small")
        r = run_cli(["--aln", ds.aln, "--out", tmp_out, "--partitions", ds.parts]); assert r.returncode == 2 and "--datatype" in r.stderr
        r = run_cli(["--aln", ds.aln, "--datatype", "dna", "--partitions", ds.parts]); assert r.returncode == 2 and "--out" in r.stderr
        r = run_cli(["--aln", ds.aln, "--out", tmp_out, "--datatype", "rna", "--partitions", ds.parts]); assert r.returncode == 2

    @pytest.mark.parametrize("ncomp", ["2", "4", "1", "0"])
    def test_n_comp_must_be_odd_and_at_least_three(self, tmp_out, ncomp):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--n-comp", ncomp])
        assert r.rc != 0 and "odd" in r.log

    def test_both_mask_restrictions_together(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--mask-similar-only", "--mask-different-only"])
        assert r.rc != 0 and "mutually exclusive" in r.log

    def test_datatype_must_match_a_nexus_declaration(self, tmp_out):
        ds = dataset("cannon_aa"); r = run_tool(ds.aln, tmp_out, "dna")
        assert r.rc != 0 and "declares" in r.log

    def test_too_few_loci(self, tmp_path, tmp_out):
        ds = dataset("dna_small"); p = tmp_path / "p.txt"
        p.write_text("".join(f"l{i} = {a}-{b}\n" for i, (a, b) in enumerate(ds.truth["locus_bounds"][:9])))
        r = run_tool(ds.aln, tmp_out, "dna", str(p)); assert r.rc != 0 and "at least 10" in r.log

    def test_too_few_taxa(self, tmp_path, tmp_out):
        ds = dataset("dna_small"); names, seqs = lo.read_fasta(open(ds.aln).read()); f = tmp_path / "four.fasta"
        f.write_text("".join(f">{n}\n{s}\n" for n, s in list(zip(names, seqs))[:4]))
        r = run_tool(str(f), tmp_out, "dna", ds.parts); assert r.rc != 0 and "at least 5" in r.log

    def test_fasta_needs_partitions(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna"); assert r.rc != 0 and "partitions" in r.log

    def test_the_output_directory_is_created(self, tmp_path):
        ds = dataset("dna_clean"); out = str(tmp_path / "a" / "b" / "c"); r = run_tool(ds.aln, out, "dna", ds.parts)
        assert r.rc == 0 and os.path.isdir(out)


@pytest.mark.e2e
class TestOptions:
    def test_n_comp_sets_k_and_the_majority(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--n-comp", "3"]); o = r.report.man["options"]
        assert o["n_comp"] == 3 and o["K_eff"] == 3 and o["M_eff"] == 2 and r.report.arr("cells/comparators.npy").shape[2] == 3

    def test_n_comp_shrinks_to_the_largest_odd_below_the_taxon_count(self, tmp_out):
        ds = dataset("aa_clean"); r = run_tool(ds.aln, tmp_out, "aa", ds.parts, ["--n-comp", "99"]); o = r.report.man["options"]
        assert o["n_comp"] == 99 and o["K_eff"] == 19 and o["M_eff"] == 10 and "shrunk" in r.log

    def test_alpha_moves_the_bar(self, tmp_out, tmp_path):
        ds = dataset("dna_clean"); r1 = run_tool(ds.aln, tmp_out, "dna", ds.parts); r2 = run_tool(ds.aln, str(tmp_path / "b"), "dna", ds.parts, ["--alpha", "0.001"])
        t1, t2 = r1.report.man["test"], r2.report.man["test"]
        assert t1["tests"] == t2["tests"] and np.isclose(t2["log10_bar"] - t1["log10_bar"], np.log10(0.001 / 0.05), atol=1e-3)
        assert r2.report.man["options"]["alpha"] == 0.001

    def test_mask_similar_only_leaves_divergent_sequences_unmasked(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--mask-similar-only"]); side = r.report.arr("cells/side.npy")
        assert (side & 1).sum() == 0 and (side & 16).sum() > 0 and (side & 2).sum() > 0
        assert all(row["reason"] == "similar" for row in r.masked) and r.report.man["options"]["mask_similar_only"]
        for e in ds.ev("contam"): assert r.masked_at(e["taxon"], e["locus"]) is None

    def test_mask_different_only_leaves_copies_unmasked(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--mask-different-only"]); side = r.report.arr("cells/side.npy")
        assert (side & 2).sum() == 0 and (side & 32).sum() > 0 and all(row["reason"] == "divergent" for row in r.masked)

    def test_donor_threshold_one_always_calls_a_direction(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--donor-threshold", "1"])
        assert r.report.man["counts"]["similar_uncertain"] == 0 and r.report.man["options"]["donor_threshold"] == 1.0
        assert all(row["similar_role"] != "uncertain" for row in r.masked)

    def test_huge_donor_threshold_masks_both_members_of_every_pair(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--donor-threshold", "1e9"]); R = r.report; side = R.arr("cells/side.npy")
        assert R.man["counts"]["similar_pairs"] > 0 and R.man["counts"]["similar_sources"] == 0
        for li in range(R.nloc):
            for e in R.similar(li):
                if e["fail"]: assert e["uncertain"] and side[e["taxon"], li] & 2 and side[e["partner"], li] & 2

    def test_fit_overlap_changes_which_taxa_can_be_compared(self, tmp_out, tmp_path):
        ds = dataset("dna_small"); r1 = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--fit-overlap", "20"]); r2 = run_tool(ds.aln, str(tmp_path / "b"), "dna", ds.parts, ["--fit-overlap", "120"])
        assert r2.report.man["options"]["fit_overlap"] == 120 and r2.report.man["counts"]["cells"] <= r1.report.man["counts"]["cells"]
        assert r2.report.arr("cells/n_comp.npy").sum() < r1.report.arr("cells/n_comp.npy").sum()

    def test_name_sets_the_dataset_name(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--name", "My Data"])
        assert r.report.man["dataset"] == "My Data"

    def test_dataset_name_is_the_folder_when_the_file_is_called_alignment(self, tmp_out, tmp_path):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); assert r.report.man["dataset"] == "dna_clean"
        f = tmp_path / "primates.fasta"; f.write_text(open(ds.aln).read()); r = run_tool(str(f), str(tmp_path / "o"), "dna", ds.parts)
        assert r.report.man["dataset"] == "primates"

    def test_no_report_writes_no_report_files(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--no-report"])
        assert r.rc == 0 and not r.exists("report.html") and not r.exists("report.npz") and r.exists("masked.csv") and r.exists("alignment.masked.fasta")

    def test_stale_files_from_an_earlier_run_are_named(self, tmp_out):
        ds = dataset("dna_clean"); run_tool(ds.aln, tmp_out, "dna", ds.parts); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--no-report"])
        assert "report.html" in r.log and "earlier run" in r.log and r.exists("report.html")

    def test_all_tests_writes_one_row_per_test(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--all-tests"]); R = r.report
        rows = list(csv.DictReader(gzip.open(r.file("tests.tsv.gz"), "rt"), delimiter="\t"))
        assert list(rows[0].keys()) == lo.TESTS_COLUMNS
        kinds = {k: sum(1 for x in rows if x["test"] == k) for k in ("divergent", "similar_comparator", "similar_closest")}
        nslots = int(R.arr("cells/n_comp.npy").sum()); assert kinds["divergent"] == nslots and kinds["similar_comparator"] == nslots
        assert kinds["similar_closest"] == R.man["test"]["similar_closest_tests"]
        bar = R.man["test"]["log10_bar"]
        for x in rows: assert (float(x["log10p"]) < bar) == (x["fail"] == "1"), x
        names = set(R.taxa["names"]); assert all(x["taxon"] in names and x["comparator"] in names for x in rows)

    def test_report_codes_policy(self, tmp_out, tmp_path):
        ds = dataset("dna_small")
        r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--report-codes", "none"]); assert not any(k.startswith("codes/") for k in r.report.names)
        r = run_tool(ds.aln, str(tmp_path / "m"), "dna", ds.parts, ["--report-codes", "masked"]); R = r.report; side = R.arr("cells/side.npy")
        for li in range(R.nloc):
            want = [g for g in range(R.n) if side[g, li] & 3]
            got = R.arr(f"codes/{li}_cells.npy").tolist() if R.has(f"codes/{li}.npy") else []
            assert got == want, li

    def test_report_sequences_policy(self, tmp_out, tmp_path):
        ds = dataset("dna_small")
        r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--report-sequences", "none"]); assert not any(k.startswith("seq/") for k in r.report.names) and r.report.man["carried"]["sequences"] == "none"
        r = run_tool(ds.aln, str(tmp_path / "m"), "dna", ds.parts, ["--report-sequences", "masked"]); R = r.report; side = R.arr("cells/side.npy"); pt = R.arr("cells/partner.npy")
        for li in range(R.nloc):
            want = set()
            for g in range(R.n):
                if side[g, li] & 3: want.add(g); want.update(R.comps(g, li)); want.add(int(pt[g, li])) if pt[g, li] >= 0 else None
            got = set(R.arr(f"seq/{li}_taxa.npy").tolist()) if R.has(f"seq/{li}.npy") else set()
            assert got == want, li

    def test_threads_zero_means_every_cpu_and_more_than_available_is_clamped(self, tmp_out, tmp_path):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--threads", "0"]); assert r.rc == 0 and r.report.man["options"]["threads"] == lo._cpus()
        r = run_tool(ds.aln, str(tmp_path / "b"), "dna", ds.parts, ["--threads", "100000"]); assert r.rc == 0 and "asked" in r.log and r.report.man["options"]["threads"] == lo._cpus()

    def test_progress_auto_without_a_terminal_prints_sparse_lines(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, progress="auto")      # stderr is captured: not a terminal
        lines = [l for l in r.log.split("\n") if l.startswith("progress\t")]
        stages = [l.split("\t")[1] for l in lines]
        for st in ("read", "counts", "pairs0", "rates", "pairs", "comparators", "dispersion", "test", "write"): assert st in stages, st
        assert "\r" not in r.log and not any(l.startswith("event\t") for l in r.log.split("\n"))

    def test_progress_lines_is_no_longer_a_mode(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, progress="lines"); assert r.rc != 0 and "invalid choice" in r.log

    def test_progress_bar_mode_without_a_terminal(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, progress="bar")
        assert r.rc == 0 and "\r" in r.log and "test" in r.log


@pytest.mark.e2e
class TestRunRecord:
    def test_log_starts_with_the_version_and_the_command(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts)
        head = r.log.split("\n")[:2]
        assert head[0].startswith(f"LocusMasker {lo.__version__}") and head[1].startswith("command: ") and "--aln" in head[1]
        assert f"{ds.ntaxa} taxa" in r.log and f"{ds.nloci} loci" in r.log

    def test_outputs_are_listed_at_the_end(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts)
        assert "outputs:" in r.log and "masked.csv" in r.log and "report.html" in r.log

    def test_manifest_records_the_input_and_options(self, tmp_out):
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--alpha", "0.01", "--threads", "2"]); m = r.report.man
        assert m["input"]["file"] == "alignment.fasta" and m["input"]["sha256"] == lo._sha256(ds.aln) and m["input"]["partitions"] == "partitions.txt" and m["input"]["format"] == "fasta"
        assert m["options"]["alpha"] == 0.01 and m["options"]["threads"] == 2 and m["options"]["datatype"] == "dna" and m["tool"] == "LocusMasker" and m["version"] == lo.__version__
        assert "--alpha 0.01" in m["command"] and set(m["timings"]) >= {"total", "test", "counts"}

    def test_a_real_process_run(self, tmp_out):
        ds = dataset("dna_clean"); r = run_cli(["--aln", ds.aln, "--partitions", ds.parts, "--datatype", "dna", "--out", tmp_out, "--progress", "none"])
        assert r.returncode == 0 and os.path.exists(os.path.join(tmp_out, "masked.csv")) and "done" in r.stderr and r.stdout == ""

    def test_a_run_emits_no_warnings(self, tmp_out):
        import warnings
        ds = dataset("dna_small")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts)
        assert r.rc == 0 and [str(x.message) for x in w if "locus_outliers" in (x.filename or "") and not issubclass(x.category, ResourceWarning)] == []
