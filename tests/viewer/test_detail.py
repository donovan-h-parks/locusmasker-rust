"""The taxon and locus panel: the crumb, the three views and the copies blocks, against the report."""
import numpy as np
import pytest
from tests.viewer.conftest import open_report, hover_text, lm_str, lm_fmt

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
        page.hover(".panel[data-key=detail] .help-q"); page.wait_for_timeout(200); assert lm_str(page, "unit.note_scatter_head") in (hover_text(page) or "")

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
        assert lm_fmt(lm_str(page, "unit.copies_n"), n=ncop) in page.text_content(".panel[data-key=detail] .crumb")
        rules = page.eval_on_selector_all(".panel[data-key=detail] .rule", "els => els.map(e => e.textContent)"); assert rules and all(r.strip() for r in rules)
        counts = [int(x) for x in page.eval_on_selector_all(".panel[data-key=detail] .rule", "els => els.map(e => e.dataset.n)")]; assert counts == sorted(counts, reverse=True) and sum(counts) == ncop
        rows = page.query_selector_all(".panel[data-key=detail] .rib-row2"); assert len(rows) == ncop
        for mode in ("recipient", "pair"):
            page.select_option(".panel[data-key=detail] .crumb select", mode); page.wait_for_timeout(400)
            assert len(page.query_selector_all(".panel[data-key=detail] .rib-row2")) == ncop and page.evaluate("() => location.hash").find(f"cp={mode}") >= 0

    def test_a_single_donor_block_by_the_api(self, page, report_of):
        """the copies view holds (reachable by the API and the chord view to come; the summary's bar lists are gone)"""
        html, R = report_of("dna_small"); open_report(page, html)
        side = R.arr("cells/side.npy"); pt = R.arr("cells/partner.npy")
        donors = {}
        for g in range(R.n):
            for li in range(R.nloc):
                if side[g, li] & 2 and not side[g, li] & 8: donors[int(pt[g, li])] = donors.get(int(pt[g, li]), 0) + 1
        assert donors, "the fixture has no called copy"
        d, n = max(donors.items(), key=lambda kv: kv[1])
        page.evaluate(f"() => LM.openCopies('donor', '{d}')"); page.wait_for_selector(".panel[data-key=detail]"); page.wait_for_timeout(400)
        crumb = page.text_content(".panel[data-key=detail] .crumb"); assert R.taxa["names"][d] in crumb and lm_fmt(lm_str(page, "unit.copies_of"), k=n, n=sum(donors.values())) in crumb
        assert len(page.query_selector_all(".panel[data-key=detail] .rib-row2")) == n
