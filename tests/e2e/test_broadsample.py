"""The broad-sample test end to end, on the fixture built for it (dna_rc): a clade moved wholesale is caught by
the far shortfall (blue, the donor named), contamination hidden where the near side of the tree has no coverage is
the pinned blind spot (the excess arm was removed), clean data raises nothing, the family enlarges the run's single
bar, and --no-broad-sample switches all of it off."""
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


def broad_masked(row):
    """masked by the broad-sample test itself: broad columns filled, no other family claims the mask"""
    return row["broad_partner"] != "" and row["divergent_fails"] == "" and \
        row["comparator_role"] == "" and row["closest_role"] == ""


class TestOnThePlantedFixture:
    def test_cladecopy_recipients_are_caught_too_close(self, rc_run):
        """a distant donor's sequence over a whole clade: each member's closest relatives moved with it, so only the
        broad sample sees the young coalescence; most recipients are masked blue with the event's locus"""
        ds, r = rc_run; rows = masked_rows(r)
        events = ds.ev("cladecopy"); recips = {(e["locus"], t) for e in events for t in e["recipients"]}
        fired = [row for row in rows if row["broad_partner"] != ""]
        hits = [row for row in fired if (row["locus"], row["taxon"]) in recips]
        assert len(hits) >= 4, [(_r["locus"], _r["taxon"]) for _r in fired]
        assert all(row["reason"] in ("similar", "both") and float(row["broad_log10p"]) < -1.3 for row in hits)
        assert all(float(row["broad_observed"]) < float(row["broad_expected"]) for row in hits)
        donors = {e["donor"] for e in events}
        assert any(row["broad_partner"] in donors for row in hits)          # the strongest too-close comparator names the donor

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
            assert row["reason"] == "divergent" and not broad_masked(row)

    def test_unplanted_masks_are_collateral_at_event_loci(self, rc_run):
        """moving a clade orphans its sisters (their closest relatives all moved), so a few unplanted cells fail the
        main tests AT THE EVENT LOCI; the broad-sample test itself adds no unplanted mask"""
        ds, r = rc_run; rows = masked_rows(r)
        planted = {(e["locus"], t) for e in ds.events for t in (e.get("recipients") or e.get("taxa") or [e.get("taxon")]) if t}
        event_loci = {e["locus"] for e in ds.events if e["kind"] in ("cladecopy", "hiddencontam", "contam")}
        stray = [q for q in rows if (q["locus"], q["taxon"]) not in planted]
        assert all(q["locus"] in event_loci for q in stray) and len(stray) <= 8
        assert not any(broad_masked(q) for q in stray)

    def test_the_family_includes_one_test_per_present_sequence(self, rc_run):
        ds, r = rc_run; man = r.report.man
        t = man["test"]
        assert t["broad_sample_tests"] > 0 and t["tests"] == 2 * t["divergent_tests"] + t["closest_tests"] + t["broad_sample_tests"]
        assert man["broad_sample"]["tested"] > 0 and man["broad_sample"]["family"] == t["broad_sample_tests"]
        assert man["side2_bits"]["broad_close"] == 4 and "broad_divergent" not in man["side2_bits"]

    def test_report_carries_the_bs_members(self, rc_run):
        ds, r = rc_run
        for k in ("cells/side2.npy", "cells/bs_locus.npy", "cells/bs_taxon.npy", "cells/bs_lo_log10p.npy",
                  "cells/bs_partner.npy", "cells/bs_var_log10p.npy"):
            assert k in r.report.names, k
        rl = r.report.arr("cells/bs_locus.npy"); assert len(rl) >= 5


class TestCleanAndOff:
    def test_clean_data_raises_nothing(self, tmp_out):
        ds = dataset("dna_clean"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, []); assert r.rc == 0, r.log
        s2 = r.report.arr("cells/side2.npy")
        assert not np.any(s2 & (4 | 8 | 64)) and (s2 & 1).any()
        assert all(row["broad_partner"] == "" for row in masked_rows(r))

    def test_the_switch_turns_the_test_off(self, tmp_out):
        ds = dataset("dna_rc"); r = run_tool(ds.aln, tmp_out, "dna", ds.parts, ["--no-broad-sample"]); assert r.rc == 0, r.log
        assert r.report.man["test"]["broad_sample_tests"] == 0 and r.report.man["broad_sample"] is None
        assert "cells/bs_locus.npy" not in r.report.names
        assert not np.any(r.report.arr("cells/side2.npy"))
        assert all(row["broad_partner"] == "" for row in masked_rows(r))
        assert "broad sample:" not in r.log
