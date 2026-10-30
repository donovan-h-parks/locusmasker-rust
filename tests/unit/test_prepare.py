"""prepare(): the per-locus pair counts, plateaus, the exclusion of thin taxa, the pair divergences m and the locus rates,
checked against brute force on the fixture alignments; and divergence_matrix(), the entry point scripts use."""
import argparse, time
import numpy as np
import pytest
from tests.conftest import lo, dataset


def encode(rows, dtype):
    alphabet = lo.DNA if dtype == "dna" else lo.PROTEIN; lut = np.full(256, -1, np.int8)
    for i, c in enumerate(alphabet): lut[ord(c)] = i; lut[ord(c.lower())] = i
    if dtype == "dna": lut[ord("U")] = lut[ord("u")] = 3
    return np.stack([lut[np.frombuffer(r, dtype=np.uint8)] for r in rows])


@pytest.fixture(scope="module")
def prepared():
    """prepare() on dna_small with one thread, plus everything needed to check it"""
    ds = dataset("dna_small"); names, rows, loci, fmt, declared, nexus = lo.read_input(ds.aln, ds.parts)
    a = argparse.Namespace(fit_overlap=20, threads=1); P = lo.prepare(a, names, rows, loci, "dna", 5, time.time())
    return dict(ds=ds, names=names, rows=rows, loci=loci, P=P, X=encode(rows, "dna"))


class TestCounts:
    def test_residue_codes(self, prepared):
        assert np.array_equal(prepared["P"].X, prepared["X"])

    def test_shared_and_same_counts_are_exact(self, prepared):
        X, loci, P = prepared["X"], prepared["loci"], prepared["P"]; n = X.shape[0]; iu = P.iu
        thin = set(P.thin.tolist())
        for li in (0, 3, len(loci) - 1):
            sites = loci[li][1]; Xl = X[:, sites]; pres = Xl >= 0
            for k in range(0, len(iu[0]), 7):
                i, j = int(iu[0][k]), int(iu[1][k]); both = pres[i] & pres[j]
                if i in thin or j in thin: assert P.W[li, k] == 0 and P.S[li, k] == 0; continue
                assert P.W[li, k] == both.sum() and P.S[li, k] == ((Xl[i] == Xl[j]) & both).sum()

    def test_focal_sites_counts_residues(self, prepared):
        X, loci, P = prepared["X"], prepared["loci"], prepared["P"]
        want = np.array([[(X[t, s] >= 0).sum() for _, s in loci] for t in range(X.shape[0])])
        assert np.array_equal(P.focal_sites, want)

    def test_plateaus_from_residue_frequencies(self, prepared):
        X, loci, P = prepared["X"], prepared["loci"], prepared["P"]
        for li in range(len(loci)):
            v = X[:, loci[li][1]]; v = v[v >= 0]; f = np.bincount(v, minlength=4) / len(v)
            assert np.isclose(P.plateaus[li], max(100 * (1 - (f ** 2).sum()), 1.0))

    def test_tri_index_covers_every_pair_symmetrically(self, prepared):
        P = prepared["P"]; n = prepared["X"].shape[0]
        assert np.array_equal(P.tri_idx, P.tri_idx.T) and len(set(P.tri_idx[np.triu_indices(n)].tolist())) == n * (n + 1) // 2


class TestExclusion:
    def test_thin_taxa_are_those_with_too_few_covered_loci(self, prepared):
        X, loci, P, ds = prepared["X"], prepared["loci"], prepared["P"], prepared["ds"]
        cov = np.array([[(X[t, s] >= 0).sum() >= 20 for _, s in loci] for t in range(X.shape[0])]).sum(1)
        assert set(P.thin.tolist()) == set(np.flatnonzero(cov < lo.MIN_LOCI_PER_TAXON).tolist())
        assert [prepared["names"][t] for t in P.thin] == [e["taxon"] for e in ds.ev("excluded")]
        assert np.array_equal(P.present, cov)

    def test_thin_taxa_have_no_divergence_to_anyone(self, prepared):
        P = prepared["P"]
        for t in P.thin: assert np.isnan(P.m[t]).all() and np.isnan(P.m[:, t]).all()


class TestDivergence:
    def test_m_is_symmetric_with_nan_diagonal(self, prepared):
        m = prepared["P"].m; n = m.shape[0]
        assert np.isnan(np.diag(m)).all()
        off = ~np.eye(n, dtype=bool); assert np.array_equal(np.isnan(m), np.isnan(m.T)) and np.allclose(m[off], m.T[off], equal_nan=True)

    def test_m_is_a_percentage_and_finite_where_pairs_overlap(self, prepared):
        P = prepared["P"]; m = P.m; n = m.shape[0]; keep = np.ones(n, bool); keep[P.thin] = False
        sub = m[np.ix_(keep, keep)]; off = ~np.eye(keep.sum(), dtype=bool)
        assert np.isfinite(sub[off]).all() and (sub[off] >= 0).all() and (sub[off] < 100).all()

    def test_floor_never_lowers_a_divergence(self, prepared):
        P = prepared["P"]; ok = np.isfinite(P.m)
        assert np.all(P.m_hi[ok] >= P.m[ok]) and np.all(P.m_xs[ok] >= 0) and np.isnan(np.diag(P.m_hi)).all()

    def test_identical_taxa_have_zero_divergence(self):
        ds = dataset("dna_small"); names, rows, loci, *_ = lo.read_input(ds.aln, ds.parts)
        rows = [bytearray(r) for r in rows]; rows[1][:] = rows[0]                     # taxon 1 becomes a copy of taxon 0 everywhere
        P = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), names, rows, loci, "dna", 5, time.time())
        assert P.m[0, 1] == 0.0 and P.m_hi[0, 1] > 0 and P.m_xs[0, 1] == 0.0

    def test_m_orders_taxa_like_the_simulated_tree(self, prepared):
        """a taxon's closest relative by m is its closest by patristic distance, for most taxa"""
        P, ds, names = prepared["P"], prepared["ds"], prepared["names"]; pat = ds.truth["patristic"]; hits = 0; tot = 0
        for t, nm in enumerate(names):
            if t in P.thin: continue
            row = np.where(np.isfinite(P.m[t]), P.m[t], np.inf); nearest_m = int(np.argmin(row))
            pr = np.array(pat[nm], float); pr[t] = np.inf; pr[P.thin] = np.inf; nearest_pat = int(np.argmin(pr))
            tot += 1; hits += nearest_m == nearest_pat
        assert hits >= 0.8 * tot, (hits, tot)

    def test_rates_median_one(self, prepared):
        P, loci = prepared["P"], prepared["loci"]; ok = np.isfinite(P.rate) & (P.rate > 0); ls = np.array([len(s) for _, s in loci], float)
        assert ok.sum() >= 3 and np.isclose(lo.wmedian_stack(P.rate[ok, None], ls[ok, None])[0], 1.0)

    def test_rates_are_positive_where_defined_and_vary_by_locus(self, prepared):
        r = prepared["P"].rate; ok = np.isfinite(r)
        assert (r[ok] > 0).all() and r[ok].std() > 0

    def test_ranking_key_is_a_strict_order(self, prepared):
        P = prepared["P"]
        for t in range(P.mj.shape[0]):
            row = P.mj[t]; fin = row[np.isfinite(row)]; assert len(np.unique(fin)) == len(fin)

    def test_threads_do_not_change_anything(self):
        ds = dataset("dna_small"); names, rows, loci, *_ = lo.read_input(ds.aln, ds.parts)
        P1 = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), names, rows, loci, "dna", 5, time.time())
        P4 = lo.prepare(argparse.Namespace(fit_overlap=20, threads=4), names, rows, loci, "dna", 5, time.time())
        for k in ("S", "W", "m", "m_hi", "m_xs", "mj", "plateaus", "rate"):
            assert np.array_equal(getattr(P1, k), getattr(P4, k), equal_nan=True), k


class TestAlphabets:
    def test_dna_lowercase_and_u_are_residues_ambiguity_is_missing(self):
        rows = [bytearray(b"acgtU"), bytearray(b"ACGTT"), bytearray(b"NRY-?")]
        X = encode(rows, "dna")
        assert X[0].tolist() == [0, 1, 2, 3, 3] and X[1].tolist() == [0, 1, 2, 3, 3] and X[2].tolist() == [-1] * 5

    def test_protein_non_standard_letters_are_missing(self):
        X = encode([bytearray(b"ACDEFxBZ-*")], "aa")
        assert X[0].tolist()[:5] == [0, 1, 2, 3, 4] and X[0].tolist()[5:] == [-1] * 5

    def test_prepare_uses_the_same_encoding(self):
        ds = dataset("aa_small"); names, rows, loci, *_ = lo.read_input(ds.aln, ds.parts)
        P = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), names, rows, loci, "aa", 5, time.time())
        assert np.array_equal(P.X, encode(rows, "aa"))


class TestDivergenceMatrix:
    def test_returns_names_and_matrices(self):
        ds = dataset("dna_small"); names, m, m_hi = lo.divergence_matrix(ds.aln, ds.parts, "dna")
        assert names == lo.read_input(ds.aln, ds.parts)[0] and m.shape == (ds.ntaxa, ds.ntaxa) and np.allclose(m, m.T, equal_nan=True)
        assert np.all(m_hi[np.isfinite(m)] >= m[np.isfinite(m)])

    def test_rejects_a_bad_datatype(self):
        ds = dataset("dna_small")
        with pytest.raises(ValueError): lo.divergence_matrix(ds.aln, ds.parts, "rna")

    def test_rejects_a_datatype_that_contradicts_the_nexus_file(self):
        ds = dataset("cannon_aa")
        with pytest.raises(ValueError): lo.divergence_matrix(ds.aln, None, "dna")

    def test_same_as_prepare(self):
        ds = dataset("dna_small"); names, m, m_hi = lo.divergence_matrix(ds.aln, ds.parts, "dna", threads=2)
        n2, rows, loci, *_ = lo.read_input(ds.aln, ds.parts); P = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), n2, rows, loci, "dna", 5, time.time())
        assert np.array_equal(m, P.m, equal_nan=True)


class TestUnresolvedPairs:
    """an identical genome is never a comparator and never votes on a rate; an unresolved pair never votes; the rates
    come from the nearest resolved, non-identical neighbours, so every taxon duplicated still gives every locus a rate"""
    def test_identical_twins_are_never_comparators_and_the_rates_come_from_the_rest(self):
        rng = np.random.default_rng(21); n, nloc, ll = 12, 12, 120
        anc = rng.integers(0, 4, nloc * ll); X = np.tile(anc, (n, 1)); mut = rng.random(X.shape) < 0.08
        X[mut] = (X[mut] + rng.integers(1, 4, int(mut.sum()))) % 4; X = np.repeat(X, 2, axis=0)                  # every taxon twice
        rows = [bytearray(bytes(np.frombuffer(b"ACGT", dtype=np.uint8)[r])) for r in X]; names = [f"t{i:02d}" for i in range(2 * n)]
        loci = [(f"l{i}", np.arange(i * ll, (i + 1) * ll)) for i in range(nloc)]
        P = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), names, rows, loci, "dna", 5, time.time())
        for g in range(2 * n):
            assert P.m[g, g ^ 1] == 0 and P.m_xs[g, g ^ 1] == 0 and not np.isfinite(P.mj[g, g ^ 1])                 # the twin: identical, never a comparator
            others = [h for h in range(2 * n) if h not in (g, g ^ 1)]; assert np.isfinite(P.mj[g, others]).all()
        assert np.isfinite(P.rate).all() and (P.rate > 0).all()

    def test_two_clades_of_identical_sequences_get_rates_from_the_other_clade(self):
        rng = np.random.default_rng(22); nloc, ll = 12, 100; anc = rng.integers(0, 4, nloc * ll); X = np.tile(anc, (20, 1))
        mut = rng.random(nloc * ll) < 0.05; X[10:, mut] = (X[10:, mut] + rng.integers(1, 4, int(mut.sum()))) % 4
        rows = [bytearray(bytes(np.frombuffer(b"ACGT", dtype=np.uint8)[r])) for r in X]; names = [f"t{i:02d}" for i in range(20)]
        loci = [(f"l{i}", np.arange(i * ll, (i + 1) * ll)) for i in range(nloc)]
        P = lo.prepare(argparse.Namespace(fit_overlap=20, threads=1), names, rows, loci, "dna", 5, time.time())
        assert np.isfinite(P.rate).all() and (P.rate > 0).all()
        assert (P.m_xs[:10, :10][~np.eye(10, dtype=bool)] == 0).all() and (P.m_xs[:10, 10:] > 0).all()
