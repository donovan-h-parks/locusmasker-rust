"""The fixture datasets themselves: the committed files are what the generator writes, and each truth file is consistent
with its alignment."""
import filecmp, json, os, tempfile
import numpy as np
import pytest
from tests.conftest import lo, DATA, SYNTH, REAL, dataset
from tests.data import make_fixtures as MF


class TestGenerator:
    @pytest.mark.parametrize("name", SYNTH)
    def test_committed_files_match_the_generator(self, name, tmp_path):
        MF.write_dataset(str(tmp_path / name), name, MF.DATASETS[name])
        for f in sorted(os.listdir(tmp_path / name)):
            assert filecmp.cmp(tmp_path / name / f, os.path.join(DATA, "synth", name, f), shallow=False), f

    def test_every_dataset_is_in_the_registry(self):
        assert set(SYNTH) == set(MF.DATASETS) and set(REAL) == set(MF.REAL)


@pytest.mark.parametrize("name", SYNTH)
class TestTruth:
    def test_events_refer_to_real_taxa_loci_and_bounds(self, name):
        ds = dataset(name); names, rows, loci, *_ = lo.read_input(ds.aln, ds.parts); lnames = [k for k, _ in loci]
        assert ds.truth["taxa"] == len(names) and ds.truth["loci"] == len(loci) and ds.truth["sites"] == len(rows[0])
        assert [(a, b) for a, b in ds.truth["locus_bounds"]] == [(int(s[0]) + 1, int(s[-1]) + 1) for _, s in loci]
        for e in ds.events:
            for k in ("taxon", "donor", "recipient"):
                if k in e: assert e[k] in names, e
            if "locus" in e: assert e["locus"] in lnames, e

    def test_planted_taxa_are_disjoint(self, name):
        ds = dataset(name); used = []
        for e in ds.events: used += [e[k] for k in ("taxon", "donor", "recipient") if k in e]
        assert len(used) == len(set(used)), used

    def test_copies_are_copies_in_the_alignment(self, name):
        ds = dataset(name); names, rows, loci, *_ = lo.read_input(ds.aln, ds.parts); lx = {k: s for k, s in loci}
        for e in ds.events:
            if e["kind"] in ("nearcopy", "farcopy", "sistercopy"):
                s = lx[e["locus"]]; d, r = names.index(e["donor"]), names.index(e["recipient"])
                assert bytes(np.frombuffer(rows[d], dtype=np.uint8)[s]) == bytes(np.frombuffer(rows[r], dtype=np.uint8)[s]), e

    def test_patristic_matrix_is_a_metric_on_the_names(self, name):
        ds = dataset(name); P = ds.truth["patristic"]; names = list(P)
        assert names == ds.truth and False or set(names) == set(lo.read_input(ds.aln, ds.parts)[0])
        M = np.array([P[n] for n in names]); assert np.allclose(M, M.T) and (np.diag(M) == 0).all() and (M[~np.eye(len(names), dtype=bool)] > 0).all()
