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
        assert data * 4 / 3 < page < data * 4 / 3 + 600_000
