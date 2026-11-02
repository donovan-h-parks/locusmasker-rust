"""Loading a report: the page comes up with the right header, refuses what it cannot read, and copes with reports that
carry less (older formats, sparse numbers, no sequences, nothing masked)."""
import base64, os
import pytest
from tests.conftest import dataset, run_tool
from tests.viewer.conftest import open_report, inline_viewer, READY

pytestmark = pytest.mark.viewer


class TestLoads:
    def test_header_shows_the_dataset_and_the_counts(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); c = R.man["counts"]
        assert page.title() == "LocusMasker report"
        assert page.text_content(".topbar h1") == R.man["dataset"] == "dna_small"
        sub = page.text_content(".topbar .sub"); assert f"{c['masked']} of {c['cells']} sequences masked" in sub and f"{c['taxa']} taxa × {c['loci']} loci" in sub
        meta = page.text_content(".topbar .meta summary"); assert R.man["version"] in meta and R.man["date"] in meta

    def test_logo_and_favicons_are_in_the_page(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        assert page.evaluate("() => document.querySelector('.topbar img.logo').src.startsWith('data:image/svg+xml')")
        assert page.evaluate("() => Array.from(document.querySelectorAll('link[rel=icon]')).map(l => l.sizes.value + ':' + l.href.slice(0, 14))") == ["16x16:data:image/png", "32x32:data:image/png"]

    def test_the_api_the_tests_use_is_exposed(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        assert page.evaluate("() => Object.keys(window.LM).sort()") == ["App", "openCell", "openUnit", "setView", "start"]
        assert page.evaluate("() => [LM.App.D.n, LM.App.D.nloc, LM.App.D.format, LM.App.D.K, LM.App.D.M]") == [R.n, R.nloc, 8, R.K, R.M]

    def test_the_loading_placeholder_is_gone(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        assert page.query_selector(".loading") is None and page.query_selector(".topbar") is not None

    @pytest.mark.parametrize("name", ["aa_small", "dna_dense", "dna_clean", "aa_clean", "cannon_aa", "near_dna"])
    def test_every_fixture_report_opens(self, page, report_of, name):
        html, R = report_of(name); open_report(page, html)
        assert page.text_content(".topbar h1") == R.man["dataset"]

    def test_a_report_with_nothing_masked(self, page, report_of):
        html, R = report_of("dna_clean"); open_report(page, html)
        tiles = page.eval_on_selector_all(".tile .v", "els => els.map(e => e.textContent)")
        assert tiles[:3] == ["0", "0", "0"] and page.query_selector_all(".worst .col .row") == [] and "donors" not in page.text_content(".worst")

    def test_sparse_numbers_and_no_sequences(self, page, tmp_path, monkeypatch):
        from tests import conftest as C
        ds = dataset("dna_small"); monkeypatch.setattr(C.lo, "DENSE_CELLS_MAX", 10)
        r = run_tool(ds.aln, str(tmp_path / "o"), "dna", ds.parts, ["--report-sequences", "none", "--report-codes", "none"]); R = r.report
        assert not R.dense
        open_report(page, r.file("report.html"))
        assert page.text_content(".topbar h1") == "dna_small"
        g, li = [(int(a), int(b)) for a, b in zip(*(R.arr("cells/side.npy") & 3).nonzero())][0]
        page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_timeout(600)
        body = page.text_content(".panel[data-key=sequence]"); assert "not carried" in body or "neither sequences nor site codes" in body


class TestRefusals:
    @pytest.mark.allow_errors
    def test_a_newer_format_is_refused_with_the_version(self, page, report_of):
        html, R = report_of("dna_small", "format9"); page.goto("file://" + html); page.wait_for_selector(".error", timeout=30000)
        txt = page.text_content(".error"); assert "format 9" in txt and "9.9.9" in txt and "formats 1 to 8" in txt

    @pytest.mark.allow_errors
    def test_corrupted_data_shows_an_error_not_a_blank_page(self, page, report_of, tmp_path):
        html, R = report_of("dna_small"); b64 = base64.b64encode(b"this is not a zip file" * 100).decode()
        out = inline_viewer(b64, str(tmp_path / "bad.html")); page.goto("file://" + out); page.wait_for_selector(".error", timeout=30000)
        assert "incomplete" in page.text_content(".error")                    # data that does not end like the tool's ZIP: the same message as a cut file

    @pytest.mark.allow_errors
    def test_a_truncated_file_says_it_is_incomplete(self, page, report_of, tmp_path):
        """a page cut inside the bundle still boots (the viewer now sits before the data) and says the file is
        incomplete, with how far its data reaches"""
        html, R = report_of("dna_small"); txt = open(html, encoding="utf-8").read()
        b = txt.index('id="lm-bundle"'); cut = str(tmp_path / "cut_bundle.html")
        open(cut, "w", encoding="utf-8").write(txt[:b + (len(txt) - b) // 2]); page.goto("file://" + cut)
        page.wait_for_selector(".error", timeout=30000)
        body = page.text_content("body"); assert "incomplete" in body and "stops after" in body and "Loading the report" not in body

    @pytest.mark.allow_errors
    def test_a_page_cut_inside_the_viewer_falls_to_the_watchdog(self, page, report_of, tmp_path):
        """cut before the viewer script ends, nothing can boot: the template's head watchdog replaces the silent
        loading line and says to open the file in a regular browser"""
        html, R = report_of("dna_small"); txt = open(html, encoding="utf-8").read()
        v = txt.index("LocusMasker report viewer. One classic"); cut = str(tmp_path / "cut_viewer.html")
        open(cut, "w", encoding="utf-8").write(txt[:v + 2000]); page.goto("file://" + cut)
        page.wait_for_function("() => document.querySelector('#app') && document.querySelector('#app').textContent.includes('regular browser')", timeout=30000)

    def test_nothing_in_the_page_reads_undefined_or_nan(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.wait_for_timeout(300)
        shown = page.evaluate("() => document.getElementById('app').innerText")                      # rendered text only, not the inlined script
        assert "undefined" not in shown and "NaN" not in shown and "null" not in shown.split()

    def test_a_page_without_data_says_so(self, page, tmp_path):
        """the viewer template without a report: a message, no file picker (the browser version is gone)"""
        out = inline_viewer("", str(tmp_path / "empty.html")); page.goto("file://" + out); page.wait_for_selector(".error", timeout=30000)
        assert "carries no data" in page.text_content(".error") and page.query_selector("input[type=file]") is None


class TestOlderFormats:
    @pytest.mark.parametrize("variant", ["format5", "format3"])
    def test_older_reports_open_and_summarise(self, page, report_of, variant):
        html, R = report_of("dna_small", variant); open_report(page, html)
        assert page.evaluate("() => LM.App.D.format") == (5 if variant == "format5" else 3)
        tiles = page.eval_on_selector_all(".tile .v", "els => els.map(e => e.textContent)"); assert tiles[0] == str(R.man["counts"]["masked"])
        for view in ("matrix", "loci", "taxa", "table", "summary"): page.evaluate(f"() => LM.setView('{view}')"); page.wait_for_timeout(150)
