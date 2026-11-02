"""The random-comparator test end to end, on the fixture built for it (dna_rc): a clade moved wholesale is caught by
the shortfall side (blue, the donor named), contamination hidden where the near side of the tree has no coverage is
caught by the excess side (red), clean data raises nothing, the family enlarges the run's single bar, and
--no-random-comparators switches all of it off."""
import csv
import numpy as np
import pytest
from tests.conftest import lo, dataset, run_tool


@pytest.fixture(scope="module")
def rc_run(tmp_path_factory):
    ds = dataset("dna_rc")
    r = run_tool(ds.aln, str(tmp_path_factory.mktemp("rc")), "dna", ds.parts, ["--threads", "4"])
    assert r.rc == 0, r.log
    return ds, r


def masked_rows(r): return list(csv.DictReader(open(r.file("masked.csv"))))


class TestOnThePlantedFixture:
    def test_cladecopy_recipients_are_caught_too_close(self, rc_run):
        """a distant donor's sequence over a whole clade: each member's closest relatives moved with it, so only the
        random comparators see the young coalescence; most recipients are masked blue with the event's locus"""
        ds, r = rc_run; rows = masked_rows(r)
        events = ds.ev("cladecopy"); recips = {(e["locus"], t) for e in events for t in e["recipients"]}
        close = [row for row in rows if row["random_side"] == "close"]
        hits = [row for row in close if (row["locus"], row["taxon"]) in recips]
        assert len(hits) >= 4, [(_r["locus"], _r["taxon"]) for _r in close]
        assert all(row["reason"] in ("similar", "both") and float(row["random_log10p"]) < -1.3 for row in hits)
        donors = {e["donor"] for e in events}
        assert any(row["random_partner"] in donors for row in hits)        # the strongest too-close comparator names the donor

    def test_hidden_contamination_is_the_designed_blind_spot(self, rc_run):
        """a contaminated stretch that NO close relative covers: the closest-relative tests compare on the clean
        remainder, and the broad-sample excess side was prototyped and not adopted (Rob, 3 Nov 2026) - the cell passes
        the variance stage's bar but nothing masks it. This pins the limitation deliberately"""
        ds, r = rc_run; rows = masked_rows(r)
        (e,) = ds.ev("hiddencontam")
        assert not [q for q in rows if (q["locus"], q["taxon"]) == (e["locus"], e["taxon"])]
        g = r.report.taxa["names"].index(e["taxon"]); li = r.report.loci["names"].index(e["locus"])
        s2 = int(r.report.arr("cells/side2.npy")[g, li])
        assert (s2 & 3) == 3 and not (s2 & 4)                             # tested and past the bar; no mask

    def test_plain_contaminants_still_fall_to_the_main_test(self, rc_run):
        ds, r = rc_run; rows = masked_rows(r)
        for e in ds.ev("contam"):
            (row,) = [q for q in rows if (q["locus"], q["taxon"]) == (e["locus"], e["taxon"])]
            assert row["reason"] == "divergent" and row["random_side"] == ""

    def test_unplanted_masks_are_collateral_at_event_loci(self, rc_run):
        """moving a clade orphans its sisters (their closest relatives all moved), so a few unplanted cells fail the
        main tests AT THE EVENT LOCI; the random-comparator test itself adds no unplanted mask"""
        ds, r = rc_run; rows = masked_rows(r)
        planted = {(e["locus"], t) for e in ds.events for t in (e.get("recipients") or e.get("taxa") or [e.get("taxon")]) if t}
        event_loci = {e["locus"] for e in ds.events if e["kind"] in ("cladecopy", "hiddencontam", "contam")}
        stray = [q for q in rows if (q["locus"], q["taxon"]) not in planted]
        assert all(q["locus"] in event_loci for q in stray) and len(stray) <= 8
        assert all(q["random_side"] == "" for q in stray)

    def test_the_family_includes_one_test_per_present_sequence(self, rc_run):
        ds, r = rc_run; man = r.report.man
        t = man["test"]
        assert t["random_tests"] > 0 and t["tests"] == 2 * t["divergent_tests"] + t["similar_closest_tests"] + t["random_tests"]
        assert man["random_comparators"]["tested"] > 0 and man["random_comparators"]["family"] == t["random_tests"]
        assert man["side2_bits"]["random_close"] == 4 and "random_divergent" not in man["side2_bits"]

    def test_report_carries_the_rc_members(self, rc_run):
        ds, r = rc_run
        for k in ("cells/side2.npy", "cells/rc_locus.npy", "cells/rc_taxon.npy", "cells/rc_lo_log10p.npy",
                  "cells/rc_partner.npy", "cells/rc_var_log10p.npy"):
            assert k in r.report.names, k
        rl = r.report.arr("cells/rc_locus.npy"); assert len(rl) >= 5


class TestCleanAndOff:
    def test_clean_data_raises_nothing(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, []); assert r.rc == 0, r.log
        s2 = r.report.arr("cells/side2.npy")
        assert not np.any(s2 & (4 | 8 | 64)) and (s2 & 1).any()
        assert all(row["random_side"] == "" for row in masked_rows(r))

    def test_the_switch_turns_the_test_off(self, tmp_out):
        ds = dataset("dna_rc"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--no-random-comparators"]); assert r.rc == 0, r.log
        assert r.report.man["test"]["random_tests"] == 0 and r.report.man["random_comparators"] is None
        assert "cells/rc_locus.npy" not in r.report.names
        assert not np.any(r.report.arr("cells/side2.npy"))
        assert all(row["random_side"] == "" for row in masked_rows(r))
        assert "random comparators:" not in r.log
