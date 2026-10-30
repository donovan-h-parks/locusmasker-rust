"""The method against the planted truths of the synthetic datasets: what it must find, what it must not."""
import numpy as np
import pytest
from tests.conftest import lo, dataset, SYNTH

pytestmark = pytest.mark.e2e


def failing_pairs(R):
    out = {}
    for li in range(R.nloc):
        for e in R.similar(li):
            if e["fail"]: out[(R.loci["names"][li], frozenset((R.taxa["names"][e["taxon"]], R.taxa["names"][e["partner"]])))] = e
    return out


@pytest.mark.parametrize("name", [n for n in SYNTH if "clean" not in n])
class TestPlanted:
    def test_contaminants_are_masked_as_too_divergent_by_every_comparator(self, runs, name):
        run = runs(name); ds = dataset(name)
        for e in ds.ev("contam"):
            r = run.masked_at(e["taxon"], e["locus"]); assert r is not None and r["reason"] in ("divergent", "both"), e
            assert int(r["divergent_comparators"]) == run.report.K and int(r["divergent_fails"]) >= run.report.M and float(r["divergent_log10p"]) < run.report.man["test"]["log10_bar"]

    def test_far_copies_are_masked_by_both_tests_with_the_donor_named(self, runs, name):
        run = runs(name); ds = dataset(name)
        for e in ds.ev("farcopy"):
            r = run.masked_at(e["recipient"], e["locus"]); assert r is not None and r["reason"] in ("both", "similar"), e
            assert r["similar_partner"] == e["donor"] and r["similar_role"] == "copy" and r["similar_test"] == "closest", e
            assert run.masked_at(e["donor"], e["locus"]) is None, e

    def test_near_copies_are_reported_as_too_similar_pairs(self, runs, name):
        """a copy from a close relative is a failing pair; the recipient is masked as the copy, or both are masked when
        the direction is uncertain (the two fit their relatives almost equally)"""
        run = runs(name); ds = dataset(name); pairs = failing_pairs(run.report)
        for e in ds.ev("nearcopy"):
            key = (e["locus"], frozenset((e["donor"], e["recipient"]))); assert key in pairs, e
            r = run.masked_at(e["recipient"], e["locus"]); assert r is not None and r["similar_partner"] == e["donor"], e
            assert r["similar_role"] in ("copy", "uncertain"), e

    def test_sister_copies_are_reported_as_too_similar_pairs(self, runs, name):
        run = runs(name); ds = dataset(name); pairs = failing_pairs(run.report)
        for e in ds.ev("sistercopy"):
            key = (e["locus"], frozenset((e["donor"], e["recipient"])))
            assert key in pairs, e
            masked = [t for t in (e["donor"], e["recipient"]) if run.masked_at(t, e["locus"]) is not None]
            assert masked, e

    def test_the_excluded_taxon_is_untested_unmasked_and_named(self, runs, name):
        run = runs(name); ds = dataset(name); R = run.report
        for e in ds.ev("excluded"):
            g = R.tx(e["taxon"]); assert R.taxa["status"][g].startswith("excluded") and (R.arr("cells/n_comp.npy")[g] == 0).all()
            assert not any(r["taxon"] == e["taxon"] for r in run.masked) and e["taxon"] in run.log and "excluded" in run.log

    def test_the_skipped_locus_is_untested_and_named(self, runs, name):
        run = runs(name); ds = dataset(name); R = run.report
        for e in ds.ev("skipped"):
            li = R.lx(e["locus"]); assert R.loci["status"][li] != "tested" and (R.arr("cells/n_comp.npy")[:, li] == 0).all()
            assert e["locus"] in run.log and "not tested" in run.log

    def test_the_ambiguous_cell_fails_with_fewer_comparators(self, runs, name):
        run = runs(name); ds = dataset(name); R = run.report
        for e in ds.ev("ambiguous"):
            g, li = R.tx(e["taxon"]), R.lx(e["locus"]); sd = R.side(g, li)
            assert sd & 1 and sd & 4 and R.cell("n_comp", g, li) == 3 == len(e["comparators"]) and set(R.comps(g, li)) == {R.tx(c) for c in e["comparators"]}
            r = run.masked_at(e["taxon"], e["locus"]); assert r["reason"] in ("divergent", "both") and int(r["divergent_comparators"]) == 3

    def test_few_masks_beyond_the_planted_ones(self, runs, name):
        run = runs(name); ds = dataset(name)
        planted = {(e.get("taxon") or e["recipient"], e.get("locus")) for e in ds.events if e["kind"] not in ("excluded", "skipped")}
        planted |= {(e["donor"], e["locus"]) for e in ds.events if "donor" in e}
        extra = [(r["taxon"], r["locus"], r["reason"]) for r in run.masked if (r["taxon"], r["locus"]) not in planted]
        assert len(extra) <= 1, extra


@pytest.mark.parametrize("name", [n for n in SYNTH if "clean" in n])
class TestClean:
    def test_nothing_is_masked(self, runs, name):
        run = runs(name); assert run.masked == [] and run.report.man["counts"]["masked"] == 0
        assert (run.report.arr("cells/side.npy") & (1 | 2 | 16 | 32)).sum() == 0

    def test_everything_present_is_tested(self, runs, name):
        R = runs(name).report; nc = R.arr("cells/n_comp.npy"); fs = R.arr("cells/focal_sites.npy")
        assert ((nc > 0) == (fs > 0)).all() and R.man["counts"]["excluded_taxa"] == 0 and R.man["counts"]["skipped_loci"] == 0
