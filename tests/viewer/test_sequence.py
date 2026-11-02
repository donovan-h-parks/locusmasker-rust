"""The alignment panel: the header verdict, the five columns, the names' colours, the p column's stars, the site colours
on the canvas (checked pixel by pixel against the rules computed here from the carried sequences), the hovers, the
controls, the trees. Every kind of cell the fixture holds."""
import math
import numpy as np
import pytest
from urllib.parse import unquote
from tests.conftest import lo
from tests.viewer.conftest import open_report, cells_of_kind, stars, hover_text, ORANGE_INK, RED, BLUE, INK, FAINT

pytestmark = pytest.mark.viewer
RULER_H, ROW_H, RULE_H = 18, 18, 7
GAP = set(b"-?.* ")


def top_of(ri): return RULER_H + 4 + ri * ROW_H + (RULE_H if ri >= 2 else 0)


def open_cell(page, g, li):
    page.evaluate(f"() => LM.openCell({g}, {li})"); page.wait_for_selector(".panel[data-key=sequence] .seq-foot", timeout=15000); page.wait_for_timeout(300)


class Expect:
    """what the panel should show for a cell, from the report"""
    def __init__(self, R, g, li):
        self.R, self.g, self.li = R, g, li; side = R.side(g, li); self.side = side; self.bar = R.man["test"]["log10_bar"]
        self.comps = R.comps(g, li); self.Y = int(R.cell("locus_closest", g, li)); self.Yp = float(R.cell("locus_closest_log10p", g, li))
        self.pairs = {}
        for e in R.similar(li):
            if e["fail"] and g in (e["taxon"], e["partner"]):
                o = e["partner"] if e["taxon"] == g else e["taxon"]; self.pairs[o] = min(self.pairs.get(o, 0), e["log10p"])
        cd = R.cell("comp_div_log10p", g, li); self.div = {c: float(cd[k]) for k, c in enumerate(self.comps)}
        self.red = {c for c, p in self.div.items() if p < self.bar and c not in self.pairs}
        self.divP = R.cell("div_log10p", g, li); self.simP = R.cell("sim_log10p", g, li)
        self.dv, self.sm = bool(side & (1 | 16)), bool(side & (2 | 32 | 128))
        seq = R.arr(f"seq/{li}.npy"); taxa = R.arr(f"seq/{li}_taxa.npy").tolist(); self.rows = {int(t): bytes(seq[k]) for k, t in enumerate(taxa)}

    def focal_stars(self):
        if self.dv and self.sm: p = min(x for x in (self.divP, self.simP) if x is not None and np.isfinite(x))
        elif self.dv: p = self.divP
        elif self.sm: p = self.simP
        else: return "–"
        return stars(p, self.bar)

    def site_colours(self):
        """per site: the focal's cell (red, blue, plain), the candidate's, and each relative's"""
        f = self.rows[self.g]; L = len(f); counting = [c for c in self.comps if c != self.Y and c in self.rows]; cand = self.rows.get(self.Y)
        out = []
        for s in range(L):
            if f[s] in GAP: out.append(None); continue
            n = d = 0
            for c in counting:
                b = self.rows[c][s]
                if b in GAP: continue
                n += 1; d += b != f[s]
            red = n >= 3 and d >= math.ceil(n / 2); blue = red and cand is not None and cand[s] == f[s]
            out.append(("blue" if blue else "red" if red else "plain", n, d))
        return out


def rgb(hexs): return [int(hexs[i:i + 2], 16) for i in (1, 3, 5)]


@pytest.fixture(scope="module")
def ctx(report_of):
    html, R = report_of("dna_small"); kinds = cells_of_kind(R)
    cols = None
    return dict(html=html, R=R, kinds=kinds)


class TestLayout:
    def test_columns_roles_and_rows(self, page, ctx):
        R = ctx["R"]; g, li = ctx["kinds"]["divergent"]; open_report(page, ctx["html"]); open_cell(page, g, li)
        assert page.eval_on_selector_all(".seq-col .trk", "els => els.map(e => e.textContent)") == ["test type", "sequence", "Taxon name", "p<α"]
        roles = page.eval_on_selector_all(".seq-role .seqnm", "els => els.map(e => e.textContent)"); assert roles == ["focal", "closest"] + [f"relative {i + 1}" for i in range(len(R.comps(g, li)))]
        assert page.eval_on_selector_all(".seq-gut .seqnm", "els => els.map(e => e.textContent).filter(Boolean)") == ["similarity", "divergence"]
        names = page.eval_on_selector_all(".seq-names .seqnm .nm", "els => els.map(e => e.textContent)")
        assert names[0] == R.taxa["names"][g] and names[1] == R.taxa["names"][int(R.cell("locus_closest", g, li))] and names[2:] == [R.taxa["names"][c] for c in R.comps(g, li)]
        assert len(page.query_selector_all(".seq-col .rule")) == 4                      # one rule per column, between the two groups
        assert page.text_content(".seq-key").startswith("focal") and "half or more of its relatives differ" in page.text_content(".seq-key")

    def test_crumb_and_hash(self, page, ctx):
        R = ctx["R"]; g, li = ctx["kinds"]["copy"]; open_report(page, ctx["html"]); open_cell(page, g, li)
        crumb = page.text_content(".panel[data-key=sequence] .crumb"); assert R.taxa["names"][g] in crumb and R.loci["names"][li] in crumb
        page.wait_for_timeout(300); assert f"c={g},{li}" in unquote(page.evaluate("() => location.hash"))

    def test_column_resize_sticks(self, page, ctx):
        g, li = ctx["kinds"]["divergent"]; open_report(page, ctx["html"]); open_cell(page, g, li)
        w0 = page.evaluate("() => getComputedStyle(document.querySelector('.seq-names')).width")
        h = page.query_selector(".seq-names .trk .hdl"); bb = h.bounding_box(); x, y = bb["x"] + 4, bb["y"] + 8
        page.mouse.move(x, y); page.mouse.down(); page.mouse.move(x + 60, y, steps=6); page.mouse.up(); page.wait_for_timeout(200)
        w1 = page.evaluate("() => getComputedStyle(document.querySelector('.seq-names')).width"); assert float(w1[:-2]) - float(w0[:-2]) >= 50
        assert "name" in page.evaluate("() => localStorage.getItem('lm.seqcols') || ''")
        open_cell(page, *ctx["kinds"]["clean"]); assert page.evaluate("() => getComputedStyle(document.querySelector('.seq-names')).width") == w1


class TestHeaderAndColours:
    @pytest.mark.parametrize("kind", ["clean", "near_miss", "divergent", "copy", "source", "uncertain", "both", "ambiguous", "closest_tested"])
    def test_verdict_names_and_stars(self, page, ctx, kind):
        if kind not in ctx["kinds"]: pytest.skip(f"the fixture has no {kind} cell")
        R = ctx["R"]; g, li = ctx["kinds"][kind]; E = Expect(R, g, li); open_report(page, ctx["html"]); open_cell(page, g, li)
        verdict = page.text_content(".seq-line .state")
        assert verdict.startswith("Masked" if E.side & 3 else "Not masked")
        assert ("too divergent" in verdict) == E.dv and ("too similar to" in verdict) == E.sm
        if E.side & 8: assert "donor uncertain: both masked" in verdict
        elif E.side & 2: assert f"({R.taxa['names'][int(R.cell('partner', g, li))]} looks more likely to be the donor)" in verdict
        elif E.side & 128: assert f"({R.taxa['names'][g]} looks more likely to be the donor)" in verdict
        if E.side & 4: assert "fewer comparators" in verdict
        # names: focal orange; the closest row blue when significantly too similar; a relative red when significantly too divergent
        cls = page.eval_on_selector_all(".seq-names .seqnm .nm", "els => els.map(e => e.className)")
        s2 = int(R.arr("cells/side2.npy")[g, li]) if R.has("cells/side2.npy") else 0
        rc_y = bool(s2 & (4 | 64)) and R.has("cells/rc_partner.npy") and E.Y is not None and any(
            int(rl) == li and int(rt) == g and int(rp) == E.Y for rl, rt, rp in zip(R.arr("cells/rc_locus.npy"), R.arr("cells/rc_taxon.npy"), R.arr("cells/rc_partner.npy")))
        off = 1 if any("random comparator" in t for t in page.eval_on_selector_all(".seq-role .seqnm", "els => els.map(e => e.textContent)")) else 0
        assert "st-focal" in cls[0] and (("st-blue" in cls[1]) == ((E.Y in E.pairs) or rc_y))   # the closest row is also blue when it is the labelled donor side
        for k, c in enumerate(E.comps): assert ("st-red" in cls[2 + off + k]) == (c in E.red), (kind, c)
        # the p column
        p = page.eval_on_selector_all(".seq-p .pcell span", "els => els.map(e => [e.textContent, e.style.background])")
        assert p[0][0] == E.focal_stars()
        bg = p[0][1].replace(" ", ""); want = ("rgb(233,202,0)" if E.dv and E.sm else "rgb(181,23,46)" if E.dv else "rgb(47,107,216)" if E.sm else "rgb(222,220,215)")
        assert bg == want, (kind, bg, want)
        if not rc_y: assert p[1][0] == (stars(E.pairs[E.Y], E.bar) if E.Y in E.pairs else "–")
        for k, c in enumerate(E.comps): assert p[2 + off + k][0] == (stars(E.div[c], E.bar) if c in E.red else "–"), (kind, c)

    def test_hovers_carry_the_numbers(self, page, ctx):
        R = ctx["R"]; g, li = ctx["kinds"]["divergent"]; E = Expect(R, g, li); open_report(page, ctx["html"]); open_cell(page, g, li)
        titles = page.eval_on_selector_all(".seq-names .seqnm", "els => els.map(e => e.title)")
        assert f"log10 p {E.divP:.2f}".replace("-", "−") in titles[0] and f"bar {E.bar:.2f}".replace("-", "−") in titles[0]
        for k, c in enumerate(E.comps): assert "too divergent from the focal: log10 p" in titles[2 + k] and "too similar to the focal" in titles[2 + k]
        assert "closest match at this locus" in titles[1]
        assert page.get_attribute(".seq-p .trk", "title").startswith("significant")


class TestSites:
    def test_letters_mode_colours_every_visible_site_by_the_rules(self, page, ctx):
        R = ctx["R"]
        for kind in ("copy", "divergent", "near_miss", "uncertain"):
            if kind not in ctx["kinds"]: continue
            g, li = ctx["kinds"][kind]; E = Expect(R, g, li); open_report(page, ctx["html"]); open_cell(page, g, li)
            page.click(".seq-ctl button:has-text('letters')"); page.wait_for_timeout(300)
            col = page.evaluate("() => { const D = LM.App.D; return { red: D.col.redTint, blue: D.col.blueTint, cell: D.col.cell }; }"); col = {k: rgb(v) for k, v in col.items()}
            W = page.evaluate("() => document.querySelector('.seq-canvas-wrap canvas').clientWidth"); nvis = min(int(W // 12) - 1, 80)
            px = page.evaluate(f"""() => {{ const c = document.querySelector('.seq-canvas-wrap canvas'); const g = c.getContext('2d'); const out = [];
                for (let ri = 0; ri < {2 + len(E.comps)}; ri++) {{ const row = []; for (let s = 0; s < {nvis}; s++) {{ const d = g.getImageData(s * 12 + 1, {RULER_H + 4} + ri * {ROW_H} + (ri >= 2 ? {RULE_H} : 0) + 2, 1, 1).data; row.push([d[0], d[1], d[2]]); }} out.push(row); }} return out; }}""")
            sites = E.site_colours(); f = E.rows[g]; cand = E.rows.get(E.Y)
            for s in range(nvis):
                if sites[s] is None: continue
                kind_s = sites[s][0]; assert px[0][s] == (col["blue"] if kind_s == "blue" else col["red"] if kind_s == "red" else col["cell"]), (kind, "focal", s, sites[s])
                if cand is not None and cand[s] not in GAP:
                    assert px[1][s] == (col["blue"] if (kind_s != "plain" and cand[s] == f[s]) else col["cell"]), (kind, "closest", s)
                for k, c in enumerate(E.comps):
                    b = E.rows[c][s]
                    if b in GAP: continue
                    assert px[2 + k][s] == (col["red"] if b != f[s] else col["cell"]), (kind, "relative", k, s)

    def test_site_hover_counts_the_relatives(self, page, ctx):
        R = ctx["R"]; g, li = ctx["kinds"]["divergent"]; E = Expect(R, g, li); open_report(page, ctx["html"]); open_cell(page, g, li)
        page.click(".seq-ctl button:has-text('letters')"); page.wait_for_timeout(200)
        box = page.evaluate("() => { const r = document.querySelector('.seq-canvas-wrap canvas').getBoundingClientRect(); return [r.x, r.y]; }")
        sites = E.site_colours(); s = next(i for i, v in enumerate(sites[:60]) if v and v[0] != "plain")
        page.mouse.move(box[0] + s * 12 + 6, box[1] + top_of(0) + 9); page.wait_for_timeout(200)
        txt = hover_text(page); assert txt and f"{sites[s][2]} of {sites[s][1]} relatives differ from the focal" in txt and f"site{s + 1}" in txt.replace(" ", "")

    def test_dots_and_block_mode_change_the_drawing(self, page, ctx):
        g, li = ctx["kinds"]["divergent"]; open_report(page, ctx["html"]); open_cell(page, g, li)
        snap = "() => document.querySelector('.seq-canvas-wrap canvas').toDataURL()"
        page.click(".seq-ctl button:has-text('letters')"); page.wait_for_timeout(200); a = page.evaluate(snap)
        page.click(".seq-ctl input[type=checkbox]"); page.wait_for_timeout(200); b = page.evaluate(snap); assert a != b
        page.click(".seq-ctl button:has-text('whole locus')"); page.wait_for_timeout(200); c = page.evaluate(snap); assert c != b
        labels = page.eval_on_selector_all(".seq-foot span", "els => els.map(e => e.textContent)"); assert any("sites" in l for l in labels)


class TestTrees:
    def test_tips_take_the_rows_colours(self, page, ctx):
        R = ctx["R"]
        for kind in ("copy", "divergent", "uncertain"):
            if kind not in ctx["kinds"]: continue
            g, li = ctx["kinds"][kind]; E = Expect(R, g, li); open_report(page, ctx["html"]); open_cell(page, g, li)
            page.click(".seq-ctl button:has-text('trees')"); page.wait_for_selector(".trees-card svg", timeout=15000); page.wait_for_timeout(300)
            tips = page.evaluate("() => Array.from(document.querySelectorAll('.trees-card svg text')).filter(t => t.getAttribute('font-size') === '11.5').map(t => [t.textContent.split(':')[0], t.getAttribute('fill'), t.getAttribute('font-weight')])")
            assert len(tips) >= 2 * (1 + len(E.comps))
            for name, fill, weight in tips:
                t = R.tx(name) if name in R.taxa["names"] else None
                if t == g: assert fill == ORANGE_INK and weight == "700"
                elif t in E.pairs: assert fill == BLUE and weight == "700"
                elif t in E.red: assert fill == RED and weight == "700"
                elif t in E.comps: assert fill == INK
            assert page.eval_on_selector_all(".trees-card svg text", "els => els.filter(e => e.textContent.includes('expected')).length") >= 1
            page.click(".trees-card .close"); page.wait_for_timeout(100); assert page.query_selector(".trees-card") is None

    def test_tree_hover_states_the_role(self, page, ctx):
        R = ctx["R"]; g, li = ctx["kinds"]["copy"]; open_report(page, ctx["html"]); open_cell(page, g, li)
        page.click(".seq-ctl button:has-text('trees')"); page.wait_for_selector(".trees-card svg"); page.wait_for_timeout(200)
        titles = page.evaluate("() => Array.from(document.querySelectorAll('.trees-card svg text title')).map(t => t.textContent)")
        assert any(t.endswith("focal sequence") for t in titles) and any("closest match at locus; significantly closer than expected" in t for t in titles)
        assert all(("inside focal group" in t) or ("outside focal group" in t) or t.endswith("focal sequence") for t in titles)
