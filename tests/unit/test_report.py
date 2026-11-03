"""The report files: report.npz as numpy opens it, the manifest, and report.html as one self-contained page."""
import base64, json, os, re, zipfile
import numpy as np
import pytest
from tests.conftest import lo, dataset, run_tool


@pytest.mark.e2e
class TestNpz:
    def test_numpy_opens_it_like_an_npz(self, runs):
        r = runs("dna_small"); z = np.load(r.file("report.npz"))
        assert "cells/n_comp" in z.files and z["cells/n_comp"].shape == (30, 14) and z["taxa/closest"].shape[1] == 5

    def test_members_are_version_1_npy_or_json(self, runs):
        z = zipfile.ZipFile(runs("dna_small").file("report.npz"))
        for k in z.namelist():
            b = z.read(k)
            if k.endswith(".npy"): assert b[:6] == b"\x93NUMPY" and b[6] == 1
            else: assert k.endswith(".json") and json.loads(b) is not None

    def test_manifest_side_bits_and_code_map_are_the_documented_ones(self, runs):
        m = runs("dna_small").report.man
        assert m["side_bits"] == dict(masked_divergent=1, masked_similar=2, ambiguous=4, similar_uncertain=8, flagged_not_masked=16, similar_not_masked=32, similar_tested=64, similar_source=128)
        assert m["code_map"] == dict(mismatch_count_max=5, similar_base=128, tick=254, absent=255)
        assert m["carried"]["numeric_columns"] == "dense" and m["carried"]["best_match_cells"] == "all"

    def test_manifest_records_the_python_and_numpy_and_no_platform(self, runs):
        m = runs("dna_small").report.man; assert "platform" not in m and m["python"] and m["numpy"] == np.__version__


@pytest.mark.e2e
class TestHtml:
    def test_the_page_is_one_self_contained_file_with_the_data_inside(self, runs):
        r = runs("dna_small"); html = open(r.file("report.html"), encoding="utf-8").read(); data = open(r.file("report.npz"), "rb").read()
        assert base64.b64encode(data).decode("ascii") in html
        assert 'id="lm-bundle"' in html and not r.exists("report_data.js") and '<script src="report_data.js">' not in html
        assert not re.search(r'<script src="[^"]+"></script>', html) and not re.search(r'<link rel="stylesheet" href=', html)
        assert "@@BUNDLE@@" not in html and "@@SIDECAR@@" not in html

    def test_icons_are_inlined(self, runs):
        html = open(runs("dna_small").file("report.html"), encoding="utf-8").read()
        assert 'rel="icon"' in html and html.count('href="data:image/png;base64,') >= 2 and 'src="data:image/svg+xml;base64,' in html
        assert 'icons/' not in html

    def test_viewer_and_style_are_inlined(self, runs):
        html = open(runs("dna_small").file("report.html"), encoding="utf-8").read()
        assert "LocusMasker report viewer" in html and "<style>" in html and "d3" in html and "<title>LocusMasker report</title>" in html

    def test_the_data_always_goes_inside_the_page(self, runs):
        """the report is one file whatever its size (Rob, 29 Oct 2026): no size cap, no sidecar"""
        assert not hasattr(lo, "REPORT_INLINE_MAX")
        r = runs("dna_clean"); html = open(r.file("report.html"), encoding="utf-8").read()
        assert base64.b64encode(open(r.file("report.npz"), "rb").read()).decode("ascii") in html and not r.exists("report_data.js")

    def test_the_page_size_is_the_data_plus_a_small_viewer(self, runs):
        r = runs("dna_small"); page = os.path.getsize(r.file("report.html")); data = os.path.getsize(r.file("report.npz"))
        assert data * 4 / 3 < page < data * 4 / 3 + 1_300_000                 # the viewer, d3, fuse and the embedded Inter face


class TestOneFileCeiling:
    """report.html carries the bundle up to REPORT_ONE_FILE_MAX; above it the run says so at the start and the end, and
    the tables stand in; the files are whole or absent; every run leaves run.log"""
    def test_above_the_ceiling_the_tables_stand_in_for_the_page(self, runs, tmp_out, monkeypatch):
        monkeypatch.setattr(lo, "REPORT_ONE_FILE_MAX", 1)
        ds = dataset("dna_small"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); assert r.rc == 0, r.log
        assert not r.exists("report.html") and r.exists("report.npz") and "report.html not written" in r.log and "above the one-file ceiling" in r.log
        assert "likely to exceed" in r.log and r.log.index("likely to exceed") < r.log.index("pair sweep:")
        taxa = open(r.file("taxa.tsv")).read().splitlines(); loci = open(r.file("loci.tsv")).read().splitlines()
        with zipfile.ZipFile(r.file("report.npz")) as Z: man = json.loads(Z.read("manifest.json")); tj = json.loads(Z.read("taxa.json")); lj = json.loads(Z.read("loci.json"))
        assert len(taxa) == ds.ntaxa + 1 and len(loci) == ds.nloci + 1 and taxa[0].startswith("taxon\tstatus\tloci_present") and loci[0].startswith("locus\tlength\tstatus")
        assert [l.split("\t")[0] for l in taxa[1:]] == tj["names"] and [l.split("\t")[0] for l in loci[1:]] == lj["names"]
        assert json.load(open(r.file("summary.json"))) == man
        cols = taxa[0].split("\t"); row = dict(zip(cols, taxa[1].split("\t"))); assert row["status"] == tj["status"][0] and int(row["loci_present"]) == tj["loci_present"][0]
        assert all(c in cols for c in ("masked_divergent", "masked_similar", "short_of_comparators", "closest_1"))
        assert not [f for f in os.listdir(tmp_out) if f.endswith(".part")] and "taxa.tsv" in r.log.split("outputs:")[1] and "report.html" not in r.log.split("outputs:")[1]
        base = runs("dna_small"); assert open(r.file("masked.csv")).read() == open(base.file("masked.csv")).read()

    def test_a_normal_run_writes_whole_files_and_a_log(self, runs):
        r = runs("dna_small")
        assert r.exists("report.html") and r.exists("run.log") and not r.exists("taxa.tsv") and not r.exists("summary.json")
        assert not [f for f in os.listdir(r.out) if f.endswith(".part")]
        log = open(r.file("run.log"), encoding="utf-8").read(); assert "LocusMasker" in log and "outputs:" in log and "run.log" in log.split("outputs:")[1]

    def test_the_tables_of_an_earlier_run_are_removed(self, tmp_out, monkeypatch):
        ds = dataset("dna_small")
        with monkeypatch.context() as mp:
            mp.setattr(lo, "REPORT_ONE_FILE_MAX", 1); r = run_tool(ds.aln, tmp_out, "dna", ds.parts); assert r.rc == 0 and r.exists("taxa.tsv")
        r = run_tool(ds.aln, tmp_out, "dna", ds.parts); assert r.rc == 0 and r.exists("report.html") and not r.exists("taxa.tsv") and not r.exists("loci.tsv") and not r.exists("summary.json")
