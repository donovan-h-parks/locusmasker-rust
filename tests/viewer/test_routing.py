"""The URL hash as the page's state, in both directions, and the panel chrome."""
import pytest
from urllib.parse import unquote
from tests.viewer.conftest import open_report, cells_of_kind

pytestmark = pytest.mark.viewer


class TestHash:
    def test_opens_a_taxon_with_a_view(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "t=3&dv=matches"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(300)
        assert R.taxa["names"][3] in page.text_content(".panel[data-key=detail] .crumb") and page.eval_on_selector(".panel[data-key=detail] .crumb select", "e => e.value") == "matches"

    def test_opens_a_locus(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "l=2"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(300)
        assert R.loci["names"][2] in page.text_content(".panel[data-key=detail] .crumb")

    def test_opens_a_cell(self, page, report_of):
        html, R = report_of("dna_small"); g, li = cells_of_kind(R)["divergent"]; open_report(page, html, f"c={g},{li}"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot"); page.wait_for_timeout(200)
        crumb = page.text_content(".panel[data-key=sequence] .crumb"); assert R.taxa["names"][g] in crumb and R.loci["names"][li] in crumb

    def test_opens_a_view_and_an_order(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "view=table&o=name,length"); page.wait_for_timeout(400)
        assert page.evaluate("() => LM.App.view") == "table" and page.query_selector(".tfoot") is not None
        assert page.evaluate("() => [LM.App.order.taxaMode, LM.App.order.lociMode]") == ["name", "length"]

    def test_a_hash_from_another_dataset_is_ignored(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "d=0000000000&t=3"); page.wait_for_timeout(300)
        assert page.query_selector(".panel[data-key=detail]") is None

    def test_the_hash_follows_the_actions(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        page.evaluate("() => LM.setView('matrix')"); page.wait_for_timeout(300); h = page.evaluate("() => location.hash")
        assert "view=matrix" in h and f"d={R.man['input']['sha256'][:10]}" in h
        page.evaluate("() => LM.openUnit('locus', 1)"); page.wait_for_timeout(400); assert "l=1" in page.evaluate("() => location.hash")
        g, li = cells_of_kind(R)["copy"]; page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_timeout(400); assert f"c={g},{li}" in unquote(page.evaluate("() => location.hash"))

    def test_copies_hash(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "cp=recipient"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(300)
        assert page.eval_on_selector(".panel[data-key=detail] .crumb select", "e => e.value") == "recipient" and "copies" in page.text_content(".panel[data-key=detail] .crumb")


class TestPanelChrome:
    def test_close_collapse_and_maximise(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.openUnit('taxon', 1)"); page.wait_for_selector(".panel[data-key=detail]")
        page.click(".panel[data-key=detail] .tools button[title=collapse]"); page.wait_for_timeout(100); assert page.evaluate("() => document.querySelector('.panel[data-key=detail]').classList.contains('collapsed')")
        page.click(".panel[data-key=detail] .tools button[title=collapse]"); page.click(".panel[data-key=detail] .tools button[title=maximise]"); page.wait_for_timeout(100)
        assert page.evaluate("() => document.querySelector('.panel[data-key=detail]').classList.contains('maximised')")
        page.click(".panel[data-key=detail] .tools button[title=maximise]"); page.click(".panel[data-key=detail] .tools button[title=close]"); page.wait_for_timeout(200)
        assert page.query_selector(".panel[data-key=detail]") is None

    def test_the_sequence_panel_closes_and_reopens(self, page, report_of):
        html, R = report_of("dna_small"); g, li = cells_of_kind(R)["divergent"]; open_report(page, html)
        page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot")
        page.click(".panel[data-key=sequence] .tools button[title=close]"); page.wait_for_timeout(200); assert page.query_selector(".panel[data-key=sequence]") is None
        page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot")


class TestOlderFormatPanels:
    def test_format5_alignment_has_no_p_column_and_says_not_recorded(self, page, report_of):
        html, R = report_of("dna_small", "format5"); kinds = cells_of_kind(report_of("dna_small")[1])
        open_report(page, html); g, li = kinds["clean"]
        page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot"); page.wait_for_timeout(200)
        assert page.query_selector(".seq-p") is None
        names = page.eval_on_selector_all(".seq-names .seqnm .nm", "els => els.map(e => e.textContent)"); assert names[1] == "not recorded"
        g, li = kinds["copy"]; page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_timeout(400)
        names = page.eval_on_selector_all(".seq-names .seqnm .nm", "els => els.map(e => [e.textContent, e.className])"); assert "st-blue" in names[1][1]
        page.click(".seq-ctl button:has-text('trees')"); page.wait_for_selector(".trees-pane svg", timeout=15000)

    def test_format3_opens_the_panels(self, page, report_of):
        html, R = report_of("dna_small", "format3"); kinds = cells_of_kind(report_of("dna_small")[1]); open_report(page, html)
        page.evaluate("() => LM.openUnit('taxon', 2)"); page.wait_for_selector(".panel[data-key=detail]")
        g, li = kinds["divergent"]; page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot")
        assert page.text_content(".seq-line .state").startswith("Masked")
