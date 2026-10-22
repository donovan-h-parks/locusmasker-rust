# locus_outliers: masking phylogenetic outliers from a multi-locus alignment

# Why?

In a large multi-locus alignment, it's possilbe that some loci seem very different from others in the same taxon. This can occur for a lot of reasons, including: contamination, mis-labelling, bioinformatic errors (like BUSCO pulling out a different locus), hybridsation or HGT, etc. If the aim is to produce a single estimate of a species tree from a concatenated or ILS-based analysis, these loci can only mislead the analysis. This method attempts to find and mask them from a single multiple sequence alignment, before the expense of building any trees. 

# Input/output

The input is the concatenated alignment and the locus boundaries. This can be supplied in one of two ways (details here: https://iqtree.github.io/doc/Advanced-Tutorial):

1. A single NEXUS format alignment file with a SETS block
2. A FASTA or PHYLIP format alignment with a RAxML-style partitions.txt file.

The output consists of the following:

1. The alignment in its input format with every outlier locus (locus × taxon cell) masked as gaps: `<name>.masked.<ext>`
2. `masked.csv`, one row per masked cell with the evidence for masking it: locus, taxon, sites masked, number of comparators, number of comparators failed, median residual, whether the cell was ambiguous (fewer than five comparators), and the closest comparator with its alignment-wide identity, overlap, locus identity, expected identity and residual
3. `plots.pdf`, the diagnostic figures (see below)
4. The full results behind the masking:
   - `cells.tsv`: every tested cell, with the same columns as `masked.csv` plus the flag, so unmasked cells can be inspected too
   - `pairs.tsv`: every individual comparison (locus, taxon, comparator rank, comparator, alignment-wide identity, overlap, locus identity, expected identity, residual, pass/fail)
   - `per_locus.tsv`: per locus, the number of sites, taxa, tested cells and pairs, the locus rate, the robust spread of residuals, the band, and the counts flagged by number of comparators failed
   - `trends.tsv`: the locus trend on a grid, as plotted
   - `ani.npy` and `taxa.txt`: the taxon × taxon alignment-wide identity matrix (%) and its row order

## Method overview

The basic idea of the method is to compare a genome to its closest relative in the alignment, and then look at the distribution of relatedness for each locus in turn. An outlier loci is one that is far _less_ related than expected, based on this comparison. We improve the power of the method by developing the expectation of the relatedness of each locus across all genomes, using the five closest relatives to each genome. We consider an outlier to be any locus that is 25% less related than the expected median relatedness, based on the full-alignment to single-locus comparisons for that locus.

A single locus in a single taxon is marked as an outlier if it is found to be an outlier when compared to three or more of the five closest relatives.

In a bit more detail, the method does the following things:

**Alignment-wide identity.** Calculate the pairwise identity (ignoring amiguities and gaps) of all taxa in the alignment. For taxa *g* and *h*, *A(g,h)* is the fraction of sites at which both
have a residue and the residues are identical, over the whole alignment.

**Comparators.** For locus *ℓ* and genome *g* with residues in *ℓ*, the comparator set *C(g,ℓ)*
is the five genomes with the highest *A(g,·)* among those that overlap *g* in *ℓ* by at least
min(20, *n_g,ℓ*) sites, where *n_g,ℓ* is the number of residues *g* has in *ℓ*. Ties are broken
at random. Fewer than five may exist.

**Locus identity.** *I_ℓ(g,h)* is the identity of *g* and *h* on the sites they share within *ℓ*.

**Trend.** Locus identity rises with genome-wide identity along a curve specific to the locus.
*T_ℓ(a)* is the running median of *I_ℓ(g,h)* over all comparator pairs of the locus with at
least 20 shared sites, as a function of *x* = log₁₀(100 − *A(g,h)*): windows of 0.1 decade,
each widened until it holds 40 pairs, joined by linear interpolation.

**Pair test.** For each pairing of a single locus *ℓ* in genomes *g*, *h* the residual is

  *r_ℓ(g,h)* = *I_ℓ(g,h)* − *T_ℓ(A(g,h))*,

This is the pair's identity in the locus minus the identity the locus trend predicts at the pair's
alignment-wide identity. 

The pair is considered an outlier when *r_ℓ(g,h)* < *m_ℓ* − *β_ℓ*, where *m_ℓ* is the median
residual of the locus (which will usually be close to zero) and *β_ℓ* = max(25, 5 *s_ℓ*), with
*s_ℓ* = 1.4826 × MAD (Median Absolute Deviation) of the residuals, which is a robust estimate of the standard deviation. In words: a pair fails when its
locus identity is below the trend by more than 25 percentage points or five robust standard
deviations, whichever is larger. The test is one-sided. 

**Cell test.** Cell (*g,ℓ*) is flagged when

  median{ *r_ℓ(g,h)* : *h* ∈ *C(g,ℓ)* } < *m_ℓ* − *β_ℓ*,

which for five comparators is the same as failing against at least three. A cell with fewer than
five comparators that satisfies the inequality is reported as *ambiguous* and treated as flagged. We use five comparators because if you use a single comparator you can't know which genome in the focal pair has the outlier locus. You only know that one of them is. Using five helps to make this a directional test. For example, a contaminant locus should appear very different from all five of the most closely related genomes being compared. While an unlucky inclusion of a contaminant locus in the _comparison set_ of a good locus will affect only a single one of the five comparisons.

## Usage

    python locus_outliers.py --aln alignment.nex --out scan
    python locus_outliers.py --aln alignment.fasta --partitions partitions.txt --out scan
    python locus_outliers.py --aln alignment.phy --partitions partitions.txt --out scan
    python plot_panels.py scan                          # scan/plots.pdf, all pages
    python plot_panels.py scan --plots log,per_locus     # any subset of log, linear, per_locus
    python plot_panels.py scan --png                     # also write each page as a PNG

A NEXUS alignment is read with Biopython's `Bio.Nexus`; each `charset` in the SETS block is one
locus. A FASTA, PHYLIP or NEXUS alignment without charsets takes a RAxML-style partition file
(`DNA, locus1 = 1-100, 250-300`, `WAG, locus2 = 101-249\3`; the model prefix is optional). A
NEXUS file with charsets and a partition file together is an error: supply one or the other.
Data type is read from the NEXUS file when it states one, otherwise inferred from the sequences
(DNA when at least 90% of the non-missing characters are A, C, G, T or U). DNA is analysed on
A, C, G, T only, with U read as T; protein on the 20 standard amino acids; gaps and every other
character are treated as missing.

Method options:

- `--n-comp 5`: comparators per cell, the closest taxa by alignment-wide identity.
- `--fit-overlap 20`: sites a pair must share in a locus to be a valid comparison and to enter
  the trend fit; a taxon with fewer residues than this in the locus needs an overlap of all its
  residues instead.
- `--min-band 25`: a pair fails when its locus identity is more than this many percentage
  points below the trend.
- `--k 5`: the band is the larger of `--min-band` and `k` robust standard deviations of the
  locus residuals.
- `--bin 0.1`: running-median window width, in decades of (100 − alignment-wide identity).
- `--bin-min 40`: pairs a window must hold; a window widens until it does.
- `--seed 1`: random seed for breaking ties between equally close comparators.

Requirements: Python 3, numpy, Biopython (NEXUS input only); matplotlib for the figures.

## Reading the panels

The cutoff is to be judged, not trusted. Flagged points hugging the dashed line in a locus with
no separate population below the cloud are the marginal cases. A trend that bends away from the
cloud at the sparse low-identity end reflects few comparators there. A cloud that is wide
everywhere usually means alignment problems rather than foreign sequences.
