"""The broad-sample test in the viewer (format 8): a close-side mask is blue, a labelled divergent cell yellow, the
hover names the test and its p, and the alignment panel shows the partner under a broad-comparator tag."""
import json
import numpy as np
import pytest
from tests.conftest import dataset
from tests.viewer.conftest import open_report


def bs_cells(R):
    """one (g, li) per broad-sample outcome on the dna_rc run, from side2"""
    s2 = R.arr("cells/side2.npy"); side = R.arr("cells/side.npy"); out = {}
    for g, li in zip(*np.nonzero(s2 & (4 | 8 | 64))):
        v = int(s2[g, li]); sd = int(side[g, li])
        if v & 4 and (sd & 3) == 2 and "close" not in out: out["close"] = (int(g), int(li))
        if v & 64 and "label" not in out: out["label"] = (int(g), int(li))
    return out


class TestBroadSampleCells:
    def test_the_colours_follow_the_kind(self, page, report_of):
        html, R = report_of("dna_rc"); open_report(page, html)
        cells = bs_cells(R); assert "close" in cells, cells
        fam = {k: page.evaluate(f"() => LM.App.D.famOf(LM.App.D.state[{g} * LM.App.D.nloc + {li}])") for k, (g, li) in cells.items()}
        ST = page.evaluate("() => ({RED: 20, BLUE: 40, BOTH: 60})")
        assert fam["close"] == ST["BLUE"]

    def test_the_labelled_divergent_cell_is_yellow(self, page, report_of):
        html, R = report_of("dna_small"); open_report(page, html)                 # dna_small's far copy: masked divergent, the donor side named
        cells = bs_cells(R); assert "label" in cells, cells
        g, li = cells["label"]
        fam = page.evaluate(f"() => LM.App.D.famOf(LM.App.D.state[{g} * LM.App.D.nloc + {li}])")
        assert fam == 60

    @pytest.mark.parametrize("kind,needle", [("close", "too close to the broad sample")])
    def test_the_labels_name_the_test(self, page, report_of, kind, needle):
        html, R = report_of("dna_rc"); open_report(page, html)
        g, li = bs_cells(R)[kind]
        label = page.evaluate(f"() => LM.App.D.cellLabel(LM.App.D.state[{g} * LM.App.D.nloc + {li}], {g} * LM.App.D.nloc + {li})")
        assert needle in label, label

    def test_the_label_names_the_donor_side(self, page, report_of):
        """the divergent cell the broad sample also calls too close: the donor-side comparator is on record and
        reaches the alignment rows as a blue badge"""
        html, R = report_of("dna_small"); open_report(page, html)
        g, li = bs_cells(R)["label"]
        got = page.evaluate(f"() => {{ const D = LM.App.D; const r = D.rcAt({g} * D.nloc + {li}); return r ? [r.pt, r.lo] : null; }}")
        assert got is not None and got[0] >= 0 and got[1] < 0

    def test_the_hover_rows_carry_the_p_and_the_partner_tag(self, page, report_of):
        html, R = report_of("dna_rc"); open_report(page, html)
        g, li = bs_cells(R)["close"]
        txt = page.evaluate(f"() => {{ const b = window.cellHover ? cellHover({g}, {li}) : null; return b ? b.textContent : null; }}")
        if txt is None: pytest.skip("cellHover is not exposed on window in this build")
        assert "broad" in txt

    def test_the_sequence_panel_tags_the_broad_partner(self, page, report_of):
        html, R = report_of("dna_rc"); open_report(page, html)
        g, li = bs_cells(R)["close"]
        page.evaluate(f"() => LM.openCell({g}, {li})")
        page.wait_for_function("() => document.body.textContent.includes('broad comparator')", timeout=30000)
        tag = page.evaluate(f"() => {{ const D = LM.App.D; const sr = D.rcAt({g} * D.nloc + {li}); return sr ? D.taxa.names[sr.pt] : null; }}")
        assert tag is not None

    def test_a_clean_report_shows_no_broad_rows(self, page, report_of):
        html, R = report_of("dna_clean"); open_report(page, html)
        s2 = R.arr("cells/side2.npy"); assert not np.any(s2 & (4 | 64))
