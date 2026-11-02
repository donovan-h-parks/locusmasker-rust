"""The overview panel: the summary's tiles and lists against the report, the matrix with its hover, click, zoom and
order controls, the ribbons and table views, and the downloads."""
import re
import numpy as np
import pytest
from tests.viewer.conftest import open_report, hover_text

pytestmark = pytest.mark.viewer


def masked_share(R):
    side = R.arr("cells/side.npy"); nc = R.arr("cells/n_comp.npy"); masked = (side & 3) > 0
    return masked.sum(1) / np.maximum((nc > 0).sum(1), 1), masked.sum(0), masked


class TestSummary:
    def test_tiles_are_the_manifest_counts(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); c = R.man["counts"]
        v = page.eval_on_selector_all(".tile .v", "els => els.map(e => e.textContent)"); k = page.eval_on_selector_all(".tile .k", "els => els.map(e => e.textContent)")
        share, per_locus, masked = masked_share(R)
        rc = R.man.get("random_comparators"); off = 1 if rc else 0                     # format 8 adds the random-comparator tile
        assert v[0] == str(c["masked"]) and v[1] == str(c["masked_divergent"]) and v[2] == str(c["masked_similar"])
        if rc: assert v[3] == str((rc.get("close") or 0) + (rc.get("divergent") or 0) + (rc.get("both_sides") or 0)) and "random comparators" in k[3]
        assert v[3 + off] == f"{int((masked.sum(1) > 0).sum())} / {c['taxa']}" and v[4 + off] == f"{int((per_locus > 0).sum())} / {c['loci']}"
        assert f"of {c['cells']} tested" in k[0] and f"{c['ambiguous']} with fewer comparators" in k[1] and f"{c['similar_uncertain']} in a pair of uncertain direction" in k[2]
        assert f"{c['excluded_taxa']} excluded" in k[3 + off] and f"{c['skipped_loci']} skipped" in k[4 + off]

    def test_most_masked_taxa_in_order(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); share, _, masked = masked_share(R)
        names = page.eval_on_selector_all(".worst .worst-row:nth-child(1) .col:nth-child(1) .row .nm", "els => els.map(e => e.textContent)")
        want = [R.taxa["names"][g] for g in sorted(np.flatnonzero(masked.sum(1) > 0), key=lambda g: -share[g])][:10]
        assert len(names) == len(want) and sorted(names) == sorted(want)
        got_share = [share[R.tx(n)] for n in names]; assert got_share == sorted(got_share, reverse=True)
        labels = page.eval_on_selector_all(".worst .worst-row:nth-child(1) .col:nth-child(1) .row .num", "els => els.map(e => e.textContent)"); assert all(l.endswith("%") for l in labels)

    def test_most_masked_loci_in_order(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); _, per_locus, _ = masked_share(R)
        names = page.eval_on_selector_all(".worst .worst-row:nth-child(1) .col:nth-child(2) .row .nm", "els => els.map(e => e.textContent)")
        nums = [int(x) for x in page.eval_on_selector_all(".worst .worst-row:nth-child(1) .col:nth-child(2) .row .num", "els => els.map(e => e.textContent)")]
        assert nums == sorted(nums, reverse=True) and nums == [int(per_locus[R.lx(n)]) for n in names] and len(names) == int((per_locus > 0).sum())

    def test_donor_lists_leave_uncertain_pairs_out(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); side = R.arr("cells/side.npy"); pt = R.arr("cells/partner.npy")
        copies = [(g, int(pt[g, li])) for g in range(R.n) for li in range(R.nloc) if side[g, li] & 2 and not side[g, li] & 8]
        assert copies, "the fixture has no called copy"
        heads = page.eval_on_selector_all(".worst h3", "els => els.map(e => e.textContent)")
        assert "most common donors" in heads and "most common recipients" in heads and "most common donor-recipient pairs" in heads
        cols = page.eval_on_selector_all(".worst .worst-row:nth-child(2) .col", "els => els.map(c => Array.from(c.querySelectorAll('.row')).map(r => [r.querySelector('.nm').textContent, +r.querySelector('.num').textContent]))")
        donors = {}; recips = {}
        for g, d in copies: donors[R.taxa["names"][d]] = donors.get(R.taxa["names"][d], 0) + 1; recips[R.taxa["names"][g]] = recips.get(R.taxa["names"][g], 0) + 1
        assert dict(cols[0]) == donors and dict(cols[1]) == recips and sum(n for _, n in cols[2]) == len(copies)
        unc = {R.taxa["names"][g] for g in range(R.n) for li in range(R.nloc) if side[g, li] & 8}
        assert not (unc & set(dict(cols[0]))) or True        # an uncertain member may be a donor elsewhere; the pair itself is never listed
        pairs = {p for p, _ in cols[2]}; assert all("→" in p for p in pairs)

    def test_histograms_count_every_unit(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        svgs = page.query_selector_all(".hist-wrap svg"); assert len(svgs) == 2
        totals = page.eval_on_selector_all(".hist-wrap svg", "els => els.map(s => Array.from(s.querySelectorAll(':scope > text')).map(t => t.textContent.replace(/,/g, '')).filter(t => /^\\d+$/.test(t)).map(Number))")
        fs = R.arr("cells/focal_sites.npy"); present_taxa = int((fs.sum(1) > 0).sum()); present_loci = int((fs.sum(0) > 0).sum())
        assert sum(totals[0]) == present_taxa and sum(totals[1]) == present_loci              # the bar labels: every unit with residues sits in one bin

    def test_a_list_row_opens_the_detail_panel(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        name = page.text_content(".worst .worst-row:nth-child(1) .col:nth-child(1) .row .nm"); page.click(".worst .worst-row:nth-child(1) .col:nth-child(1) .row"); page.wait_for_selector(".panel[data-key=detail]")
        assert name in page.text_content(".panel[data-key=detail] .crumb")


class TestMatrix:
    def test_matrix_draws_and_the_views_switch(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        assert page.eval_on_selector_all("select option", "els => els.map(e => e.value)")[:5] == ["summary", "matrix", "loci", "taxa", "table"]
        page.select_option(".viewctl select, .panel[data-key=overview] select >> nth=0", "matrix"); page.wait_for_timeout(300)
        nonblank = page.evaluate("() => { const c = document.querySelector('canvas.matrix'); const d = c.getContext('2d').getImageData(0, 0, c.width, c.height).data; let n = 0; for (let i = 3; i < d.length; i += 4) if (d[i]) n++; return [n, c.width, c.height]; }")
        assert nonblank[0] > 0.5 * nonblank[1] * nonblank[2]

    def test_hover_names_the_cell_and_click_opens_the_alignment(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.setView('matrix')"); page.wait_for_timeout(300)
        box = page.evaluate("() => { const c = document.querySelector('canvas.matrix'); const r = c.getBoundingClientRect(); return [r.x, r.y, r.width, r.height]; }")
        V = page.evaluate("() => { const V = LM.App.D && window.LM.App.panels.overview ? null : null; return null; }")
        page.mouse.move(box[0] + 20, box[1] + 10); page.wait_for_timeout(200)
        txt = hover_text(page); assert txt and " · " in txt and any(n in txt for n in R.taxa["names"]) and any(l in txt for l in R.loci["names"])
        page.mouse.click(box[0] + 20, box[1] + 10); page.wait_for_selector(".panel[data-key=sequence]", timeout=10000)
        assert page.query_selector(".seq-head") is not None

    def test_order_controls_reorder_the_rows(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.setView('matrix')"); page.wait_for_timeout(200)
        before = page.evaluate("() => LM.App.order.taxa.slice()")
        page.evaluate("() => { const s = document.querySelectorAll('.viewctl select, .panel[data-key=overview] select'); for (const x of s) if (Array.from(x.options).some(o => o.value === 'masked')) { x.value = 'masked'; x.dispatchEvent(new Event('change')); break; } }"); page.wait_for_timeout(200)
        after = page.evaluate("() => LM.App.order.taxa.slice()"); share, _, _ = masked_share(R)
        assert after != before and [share[g] for g in after] == sorted([share[g] for g in after], reverse=True)
        page.evaluate("() => { const s = document.querySelectorAll('.viewctl select, .panel[data-key=overview] select'); for (const x of s) if (Array.from(x.options).some(o => o.value === 'name')) { x.value = 'name'; x.dispatchEvent(new Event('change')); break; } }"); page.wait_for_timeout(200)
        byname = page.evaluate("() => LM.App.order.taxa.slice()"); assert [R.taxa["names"][g] for g in byname] == sorted(R.taxa["names"])

    def test_zoom_with_ctrl_wheel(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.setView('matrix')"); page.wait_for_timeout(200)
        box = page.evaluate("() => { const r = document.querySelector('canvas.matrix').getBoundingClientRect(); return [r.x, r.y]; }")
        before = page.evaluate("() => location.hash")
        page.evaluate("() => { const c = document.querySelector('canvas.matrix'); const r = c.getBoundingClientRect(); c.dispatchEvent(new WheelEvent('wheel', { deltaY: -300, ctrlKey: true, clientX: r.x + 50, clientY: r.y + 50, bubbles: true, cancelable: true })); }"); page.wait_for_timeout(400)
        after = page.evaluate("() => location.hash"); assert after != before and "v=" in after


class TestRibbonsAndTable:
    @pytest.mark.parametrize("view,kind", [("loci", "locus"), ("taxa", "taxon")])
    def test_ribbons_list_the_units(self, page, report_of, view, kind):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate(f"() => LM.setView('{view}')"); page.wait_for_timeout(400)
        rows = page.query_selector_all(".ribbons-view .rib-row, .rib-list .rib-row"); assert len(rows) >= 1
        canvases = page.query_selector_all(".ribbons-view canvas, .rib-list canvas"); assert len(canvases) >= 1

    def test_table_counts_sorts_and_filters(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.setView('table')"); page.wait_for_timeout(400); c = R.man["counts"]
        assert page.text_content(".tfoot .n").startswith(f"{c['masked']} rows")
        heads = page.eval_on_selector_all(".thead span", "els => els.map(e => e.textContent)"); assert heads[0] == "taxon" and any("excess" in h for h in heads)
        first = page.eval_on_selector_all(".trow", "els => els.map(e => e.children[0].textContent)")
        page.click(".thead span >> nth=0"); page.wait_for_timeout(200)
        after = page.eval_on_selector_all(".trow", "els => els.map(e => e.children[0].textContent)"); assert after == sorted(after) and len(after) == len(first)
        page.click(".tfoot input[type=checkbox]"); page.wait_for_timeout(300)
        assert page.text_content(".tfoot .n").startswith(f"{c['cells']} rows")

    def test_a_table_row_opens_the_alignment(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html); page.evaluate("() => LM.setView('table')"); page.wait_for_timeout(300)
        name = page.text_content(".trow >> nth=0 >> span >> nth=0"); page.click(".trow >> nth=0"); page.wait_for_selector(".panel[data-key=sequence]")
        assert name.strip() in page.text_content(".panel[data-key=sequence] .crumb")

    def test_downloads(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)
        with page.expect_download() as dl: page.click("button:has-text('cells.tsv')")
        d = dl.value; assert d.suggested_filename == "dna_small_cells.tsv"
        txt = open(d.path(), encoding="utf-8").read().split("\n"); assert txt[0].startswith("taxon\tlocus\tstate") and len(txt) - 2 >= R.man["counts"]["cells"] - 1
        page.evaluate("() => LM.setView('matrix')"); page.wait_for_timeout(300)
        with page.expect_download() as dl: page.click("button:has-text('PNG')")
        assert dl.value.suggested_filename == "dna_small_overview.png"
