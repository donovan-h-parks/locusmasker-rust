"""The taxon and locus panel: the crumb, the three views and the copies blocks, against the report."""
import numpy as np
import pytest
from tests.viewer.conftest import open_report, hover_text

pytestmark = pytest.mark.viewer


def most_masked_taxon(R):
    side = R.arr("cells/side.npy"); return int(np.argmax(((side & 3) > 0).sum(1)))


class TestTaxonPanel:
    def test_crumb_and_views(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); g = most_masked_taxon(R)
        page.evaluate(f"() => LM.openUnit('taxon', {g})"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(300)
        crumb = page.text_content(".panel[data-key=detail] .crumb"); side = R.arr("cells/side.npy"); nc = R.arr("cells/n_comp.npy")
        assert R.taxa["names"][g] in crumb and f"{int(((side[g] & 3) > 0).sum())} of {int((nc[g] > 0).sum())} sequences masked" in crumb
        opts = page.eval_on_selector_all(".panel[data-key=detail] .crumb select option", "els => els.map(e => e.value)"); assert opts == ["scatter", "matches", "ribbons"]

    def test_scatter_has_axes_diagonal_and_points(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); g = most_masked_taxon(R)
        page.evaluate(f"() => LM.openUnit('taxon', {g})"); page.wait_for_selector(".panel[data-key=detail] svg"); page.wait_for_timeout(300)
        texts = page.eval_on_selector_all(".panel[data-key=detail] svg text", "els => els.map(e => e.textContent)")
        assert any("expected divergence" in t for t in texts) and any("observed divergence" in t for t in texts)
        painted = page.evaluate("() => { const c = document.querySelector('.panel[data-key=detail] .scatter-wrap canvas'); const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data; let n = 0; for (let i = 3; i < d.length; i += 4) if (d[i]) n++; return n; }")
        assert painted > 50
        nc = R.arr("cells/n_comp.npy")
        page.hover(".panel[data-key=detail] .help-q"); page.wait_for_timeout(200); assert "Observed vs expected" in (hover_text(page) or "")

    def test_matches_bars_count_the_best_matches(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); g = most_masked_taxon(R)
        page.evaluate(f"() => LM.openUnit('taxon', {g})"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(200)
        page.select_option(".panel[data-key=detail] .crumb select", "matches"); page.wait_for_timeout(400)
        side = R.arr("cells/side.npy"); bm = R.arr("cells/best_match.npy")
        masked_bm = [int(bm[g, li]) for li in range(R.nloc) if side[g, li] & 3 and bm[g, li] >= 0]
        labels = [int(t) for t in page.eval_on_selector_all(".panel[data-key=detail] svg > text", "els => els.map(e => e.textContent)") if t.isdigit()]
        assert sum(labels) == len(masked_bm) and labels == sorted(labels, reverse=True)
        names = page.eval_on_selector_all(".panel[data-key=detail] svg .axis text", "els => els.map(e => e.textContent)")
        assert set(R.taxa["names"][b] for b in masked_bm) <= set(names) | {n for n in names if n.startswith("other")}
        assert page.query_selector_all(".panel[data-key=detail] svg > line") == []           # no expected-count ticks

    def test_ribbons_list_the_masked_cells(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); g = most_masked_taxon(R)
        page.evaluate(f"() => LM.openUnit('taxon', {g})"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(200)
        page.select_option(".panel[data-key=detail] .crumb select", "ribbons"); page.wait_for_timeout(500)
        rows = page.query_selector_all(".panel[data-key=detail] .rib-row2"); side = R.arr("cells/side.npy")
        assert len(rows) == int(((side[g] & 3) > 0).sum())
        rows[0].click(); page.wait_for_selector(".panel[data-key=sequence]"); assert R.taxa["names"][g] in page.text_content(".panel[data-key=sequence] .crumb")

    def test_locus_panel_has_no_matches_view(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); side = R.arr("cells/side.npy"); li = int(np.argmax(((side & 3) > 0).sum(0)))
        page.evaluate(f"() => LM.openUnit('locus', {li})"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(300)
        assert R.loci["names"][li] in page.text_content(".panel[data-key=detail] .crumb")
        assert page.eval_on_selector_all(".panel[data-key=detail] .crumb select option", "els => els.map(e => e.value)") == ["scatter", "ribbons"]


class TestCopies:
    def test_blocks_by_donor_recipient_and_pair(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html, "cp=donor"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(400)
        side = R.arr("cells/side.npy"); ncop = int((((side & 2) > 0) & ((side & 8) == 0)).sum())
        assert f"{ncop} copies" in page.text_content(".panel[data-key=detail] .crumb")
        rules = page.eval_on_selector_all(".panel[data-key=detail] .rule", "els => els.map(e => e.textContent)"); assert rules and all("cop" in r for r in rules)
        counts = [int(r.split("·")[1].split()[0]) for r in rules]; assert counts == sorted(counts, reverse=True) and sum(counts) == ncop
        rows = page.query_selector_all(".panel[data-key=detail] .rib-row2"); assert len(rows) == ncop
        for mode in ("recipient", "pair"):
            page.select_option(".panel[data-key=detail] .crumb select", mode); page.wait_for_timeout(400)
            assert len(page.query_selector_all(".panel[data-key=detail] .rib-row2")) == ncop and page.evaluate("() => location.hash").find(f"cp={mode}") >= 0

    def test_a_single_block_from_the_summary(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        first = page.query_selector(".worst .worst-row:nth-child(2) .col:nth-child(1) .row"); name = first.query_selector(".nm").text_content(); n = int(first.query_selector(".num").text_content())
        first.click(); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(400)
        crumb = page.text_content(".panel[data-key=detail] .crumb"); assert name in crumb and f"{n} of" in crumb
        assert len(page.query_selector_all(".panel[data-key=detail] .rib-row2")) == n
