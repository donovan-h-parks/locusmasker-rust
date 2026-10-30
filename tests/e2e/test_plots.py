"""The plotting scripts on a run's outputs: run in this process (so coverage sees them), plus one real process."""
import contextlib, io, os, runpy, subprocess, sys
import pytest
from tests.conftest import ROOT, PY, dataset, run_tool

pytestmark = pytest.mark.e2e


def script(name, args, cwd=None):
    """a plot script in this process: sys.argv set, stdout and stderr captured, SystemExit caught"""
    old_argv, old_cwd = sys.argv, os.getcwd(); out, err = io.StringIO(), io.StringIO(); rc = 0
    sys.argv = [os.path.join(ROOT, name)] + [str(a) for a in args]
    if cwd: os.chdir(cwd)
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try: runpy.run_path(os.path.join(ROOT, name), run_name="__main__")
            except SystemExit as e:
                rc = e.code if isinstance(e.code, int) else 1
                if e.code and not isinstance(e.code, int): err.write(str(e.code) + "\n")
    finally:
        sys.argv = old_argv; os.chdir(old_cwd)
        import matplotlib.pyplot as plt; plt.close("all")
    return subprocess.CompletedProcess(sys.argv, rc, out.getvalue(), err.getvalue())


def with_fails_only(run, d):
    """a copy of a run whose report says best matches were carried for failing cells only"""
    import json, zipfile, shutil
    shutil.copytree(run.out, d)
    src = zipfile.ZipFile(run.file("report.npz")); dst = zipfile.ZipFile(d / "report.npz", "w", zipfile.ZIP_STORED)
    for it in src.infolist():
        b = src.read(it.filename)
        if it.filename == "manifest.json": m = json.loads(b); m["carried"]["best_match_cells"] = "fails"; b = json.dumps(m).encode()
        dst.writestr(it, b)
    dst.close(); return str(d)


class TestPlotPanels:
    def test_every_page_renders(self, runs, tmp_path):
        run = runs("dna_small"); r = script("plot_panels.py", [run.out, "--png", "--out", str(tmp_path / "plots.pdf")], str(tmp_path))
        assert r.returncode == 0, r.stderr
        assert os.path.getsize(tmp_path / "plots.pdf") > 1000
        for page in ("taxa", "loci", "bars", "matches", "residuals_taxa", "residuals_loci"): assert (tmp_path / f"plots_{page}.png").exists(), page

    def test_a_subset_of_pages_and_the_title(self, runs, tmp_path):
        run = runs("aa_small"); r = script("plot_panels.py", [run.out, "--plots", "taxa,bars", "--title", "A title", "--out", str(tmp_path / "p.pdf"), "--png"], str(tmp_path))
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "p_taxa.png").exists() and not (tmp_path / "p_matches.png").exists()

    def test_all_units_and_columns(self, runs, tmp_path):
        run = runs("dna_dense"); r = script("plot_panels.py", [run.out, "--plots", "taxa,loci", "--all", "--cols", "5", "--out", str(tmp_path / "p.pdf")], str(tmp_path))
        assert r.returncode == 0, r.stderr

    def test_unknown_page_is_an_error(self, runs, tmp_path):
        r = script("plot_panels.py", [runs("dna_small").out, "--plots", "nope"], str(tmp_path)); assert r.returncode != 0 and "unknown" in r.stderr

    def test_matches_page_from_the_alignment(self, runs, tmp_path):
        """when the report carries best matches for failing cells only, the page recomputes them from the alignment"""
        d = with_fails_only(runs("dna_small"), tmp_path / "run")
        r = script("plot_panels.py", [d, "--plots", "matches", "--png", "--out", str(tmp_path / "m.pdf")], str(tmp_path)); assert r.returncode == 0 and "skipped" in r.stdout
        ds = dataset("dna_small")
        r = script("plot_panels.py", [d, "--plots", "matches", "--png", "--out", str(tmp_path / "m2.pdf"), "--alignment", ds.aln, "--partitions", ds.parts], str(tmp_path))
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "m2_matches.png").exists()

    def test_the_alignment_must_be_the_runs(self, runs, tmp_path):
        d = with_fails_only(runs("dna_small"), tmp_path / "run"); ds = dataset("aa_small")
        r = script("plot_panels.py", [d, "--plots", "matches", "--out", str(tmp_path / "m.pdf"), "--alignment", ds.aln, "--partitions", ds.parts], str(tmp_path))
        assert r.returncode != 0 and "run's taxa" in r.stderr

    def test_old_outputs_are_refused_with_advice(self, tmp_path):
        (tmp_path / "cells.tsv").write_text("x\n")
        r = script("plot_panels.py", [str(tmp_path)], str(tmp_path)); assert r.returncode != 0 and "re-run" in r.stderr
        r = script("plot_panels.py", [str(tmp_path / "nothing")], str(tmp_path)); assert r.returncode != 0 and "report.npz" in r.stderr

    @pytest.mark.parametrize("name", ["cannon_aa", "near_dna", "dna_dense", "aa_clean"])
    def test_pages_on_other_datasets(self, runs, tmp_path, name):
        run = runs(name); r = script("plot_panels.py", [run.out, "--plots", "bars,taxa,loci,matches", "--out", str(tmp_path / "p.pdf")], str(tmp_path))
        assert r.returncode == 0, r.stderr

    def test_a_sparse_report(self, tmp_path, monkeypatch):
        from tests.conftest import lo
        monkeypatch.setattr(lo, "DENSE_CELLS_MAX", 10); ds = dataset("dna_small"); run = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts); assert not run.report.dense
        r = script("plot_panels.py", [run.out, "--png", "--out", str(tmp_path / "p.pdf")], str(tmp_path)); assert r.returncode == 0, r.stderr
        assert (tmp_path / "p_taxa.png").exists()

    @pytest.mark.filterwarnings("ignore:Tight layout not applied")      # matplotlib cannot fit 80-character titles on the matches page: cosmetic
    def test_long_names_in_the_titles(self, tmp_path):
        ds = dataset("dna_small"); from tests.conftest import lo
        names, seqs = lo.read_fasta(open(ds.aln).read()); names = [n + "_" + "longname" * 6 + "-" + "x" * 30 for n in names]; f = tmp_path / "a.fasta"
        with open(f, "w") as fh:
            for n, sq in zip(names, seqs): fh.write(f">{n}\n{sq}\n")
        run = run_tool(str(f), str(tmp_path / "o"), "dna", ds.parts); r = script("plot_panels.py", [run.out, "--plots", "taxa,matches", "--out", str(tmp_path / "p.pdf")], str(tmp_path))
        assert r.returncode == 0, r.stderr

    def test_an_older_report_format_is_refused(self, runs, tmp_path):
        import json, zipfile, shutil
        run = runs("dna_small"); d = tmp_path / "run"; shutil.copytree(run.out, d)
        src = zipfile.ZipFile(run.file("report.npz")); dst = zipfile.ZipFile(d / "report.npz", "w", zipfile.ZIP_STORED)
        for it in src.infolist():
            b = src.read(it.filename)
            if it.filename == "manifest.json": m = json.loads(b); m["format_version"] = 3; b = json.dumps(m).encode()
            dst.writestr(it, b)
        dst.close(); r = script("plot_panels.py", [str(d), "--plots", "bars", "--out", str(tmp_path / "p.pdf")], str(tmp_path)); assert r.returncode != 0 and "format 3" in r.stderr

    def test_seven_comparators_ramp_the_colours(self, tmp_path):
        """with K=7 there are four failing levels (4..7) and four non-failing: the colour ramps are interpolated"""
        ds = dataset("dna_dense"); run = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts, ["--n-comp", "7"]); assert run.report.K == 7
        r = script("plot_panels.py", [run.out, "--plots", "bars,taxa,loci", "--out", str(tmp_path / "p.pdf")], str(tmp_path)); assert r.returncode == 0, r.stderr

    def test_a_real_process_run(self, runs, tmp_path):
        run = runs("dna_clean"); r = subprocess.run([PY, os.path.join(ROOT, "plot_panels.py"), run.out, "--plots", "bars", "--out", str(tmp_path / "b.pdf")], capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, r.stderr


class TestOtherScripts:
    def test_plot_taxon(self, runs, tmp_path):
        run = runs("dna_small"); taxon = run.masked[0]["taxon"]
        r = script("plot_taxon.py", [run.out, taxon, "--out", str(tmp_path / "t.png")], str(tmp_path)); assert r.returncode == 0, r.stderr
        assert (tmp_path / "t.png").exists()

    def test_plot_taxon_with_labels(self, runs, tmp_path):
        run = runs("dna_dense"); taxon = run.masked[0]["taxon"]
        r = script("plot_taxon.py", [run.out, taxon, "--label", "4", "--out", str(tmp_path / "t.png")], str(tmp_path)); assert r.returncode == 0, r.stderr

    def test_plot_taxon_refuses_old_outputs(self, tmp_path):
        (tmp_path / "cells.tsv").write_text("x\n")
        r = script("plot_taxon.py", [str(tmp_path), "t01"], str(tmp_path)); assert r.returncode != 0 and "re-run" in r.stderr
        r = script("plot_taxon.py", [str(tmp_path / "nothing"), "t01"], str(tmp_path)); assert r.returncode != 0 and "report.npz" in r.stderr

    def test_plot_taxon_on_a_sparse_report(self, tmp_path, monkeypatch):
        """a sparse report has numbers only for cells with a failing test: a taxon with such cells plots, one without is refused"""
        from tests.conftest import lo
        monkeypatch.setattr(lo, "DENSE_CELLS_MAX", 10); ds = dataset("dna_small"); run = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts); R = run.report; assert not R.dense
        st = set(R.arr("cells/sparse_taxon.npy").tolist()); nc = R.arr("cells/n_comp.npy"); names = R.taxa["names"]
        with_rows = next(g for g in range(R.n) if g in st); r = script("plot_taxon.py", [run.out, names[with_rows], "--out", str(tmp_path / "t.png")], str(tmp_path))
        assert r.returncode == 0 and (tmp_path / "t.png").exists(), r.stderr
        without = [g for g in range(R.n) if g not in st and (nc[g] > 0).any()]
        if without:
            r = script("plot_taxon.py", [run.out, names[without[0]], "--out", str(tmp_path / "u.png")], str(tmp_path)); assert r.returncode != 0 and "sparse" in r.stderr

    def test_plot_taxon_unknown_taxon(self, runs, tmp_path):
        r = script("plot_taxon.py", [runs("dna_small").out, "nobody", "--out", str(tmp_path / "t.png")], str(tmp_path)); assert r.returncode != 0

    def test_plot_pair(self, tmp_path):
        ds = dataset("dna_small"); run = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts, ["--all-tests"]); row = run.masked[0]; other = row["comparators"].split("; ")[0]
        r = script("plot_pair.py", [run.out, row["taxon"], other, "--out", str(tmp_path / "p.png")], str(tmp_path)); assert r.returncode == 0, r.stderr
        assert (tmp_path / "p.png").exists()
        r = script("plot_pair.py", [run.out, row["taxon"], row["taxon"]], str(tmp_path)); assert r.returncode != 0 and "two different" in r.stderr
        r = script("plot_pair.py", [run.out, row["taxon"], "nobody"], str(tmp_path)); assert r.returncode != 0 and "not a taxon" in r.stderr
        excluded = next(e["taxon"] for e in ds.events if e["kind"] == "excluded")
        r = script("plot_pair.py", [run.out, row["taxon"], excluded], str(tmp_path)); assert r.returncode != 0 and "never tested" in r.stderr

    def test_plot_pair_on_a_too_similar_pair(self, tmp_path):
        ds = dataset("dna_small"); run = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts, ["--all-tests"])
        row = next(x for x in run.masked if x["similar_partner"])
        r = script("plot_pair.py", [run.out, row["taxon"], row["similar_partner"], "--out", str(tmp_path / "p.png")], str(tmp_path)); assert r.returncode == 0, r.stderr
        assert (tmp_path / "p.png").exists()

    def test_plot_pair_needs_the_tests_file(self, runs, tmp_path):
        run = runs("dna_small"); row = run.masked[0]
        r = script("plot_pair.py", [run.out, row["taxon"], row["comparators"].split("; ")[0], "--out", str(tmp_path / "p.png")], str(tmp_path))
        assert r.returncode != 0 and "--all-tests" in r.stderr
