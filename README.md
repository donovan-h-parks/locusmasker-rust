# locus_outliers: masking phylogenetic outliers from a multi-locus alignment

# What's this?

A script which attempts to highlight and mask potentially misleading seuqences (individual loci from individual taxa) from large phylogenomic analyses. Read on if you have a large (>10 loci) dataset, and are intending to use it to estimate a single phylogenetic tree with a method that doesn't account for things like contamination, introgression, horizontal gene transfer, or accidental use of paralogs. (NB: this includes methods in software like IQ-TREE, RAxML, BEAST, PhyloBayes, MCMCTREE, MrBayes, ASTRAL, etc.)

[add summary figure - alignment of 10 loci from 7 taxa, all present to start with, big arrow, alignment of the same loci, but with a handful of cells obviously gapped out; maybe one in the middle with the outlier loci highlited in red; even better if the first alignment shows the contaminated sequences as clearly coming from a different taxon, and then shows a tree grouping these taxa; then the masked alignment shows the tree fixed]

# Why?

A lot of people (me included) are interested in estimating phylogenies from large datasets. Although we all know that some data doesn't fit on a tree, we are often willing to intentionally ignore this and as a first pass estimate a single tree topology that we hope is a good explanation of the entire dataset (e.g. a species tree). 

However, in a large multi-locus alignment not all sequenced loci from all taxa will *fit* on a single tree, and these loci can be problematic for inferring that tree. Loci can mislead tree inference for a lot of reasons, including but not limited to: contamination, mis-labelling, bioinformatic errors (like BUSCO pulling out a different locus from different genomes), hybridsation and introgression, horizontal gene transfer, accidental labelling of a paralog as an ortholog, etc. If the aim is to produce a single estimate of a tree from a concatenated (e.g. IQ-TREE), ILS-aware (e.g. BEAST), or ILS-robust (ASTRAL) analysis, these loci can mislead the analysis. This repo contains scripts which look at the alignment and attempt to find and mask such loci.

# Pros and cons of masking data

The intention of this repo is that one would use this script to mask data before building a tree. An obvious use case would be to detect and remove contaminants and mislabelled sequences from large datasets. These are _obviously_ misleading because they're _wrong_, and nobody wants them in their datasets. Another use case is more nuanced - removing loci that have moved horizontally, due to processes like introgression or HGT. These loci are _not_ wrong, but should still be interested in masking them if the intention is to use a method that doesn't account for these biological processes when building the tree. 

So, one should be aware that the individual loci from individual taxa highlighted by this method could be highlighted due to any of the following three reasons: (i) the loci are strictly wrong (e.g. contamination during sequencing); (ii) the loci are biologically correct but not accounted for by a downstream method (e.g. HGT in a standard BEAST analysis); (iii) the loci have highlighted but it's a false positive. The method below seeks to maximise identification of the first two categories, and minimise the third. 

One needs to be aware though that because categories i and ii above have identical signals in the data, if one is interested in the biological side, one can consider this method as highlighting exactly the loci that might have been invovled in the interesting biology. In this case, the sensible use of the script would be to do it in two stages:

i. use the script to mask data to build a single tree
ii. use un-masked data to then look for the interesting biological signals w.r.t to the tree estimated from the masked data (of course, the data on which loci were masked is also telling you something interesting, but the methods are intentionally rather crude and tree-free, because I wanted them to be a pre-tree step, and very fast)

# Input/output

## Input

Input either of the following:

1. A single NEXUS format alignment file with a SETS block defining the loci
2. A FASTA, PHYLIP, or NEXUS format alignment with a RAxML-style partitions.txt file defining the loci.

NB: your alignment must have at least 10 loci for this method to work.

More details on formats can be found here: https://iqtree.github.io/doc/Advanced-Tutorial

## Output

The output consists of the following:


1. `results_report.html`, an interactive HTML report showing all results [claude please note that we need to build this later, once we've finalised the plots]
2. The alignment in its input format with every outlier sequence (individual locus in an individual taxon) masked as gaps: `<name>.masked.<ext>`
3. `masked.csv`, one row per masked cell with the evidence for masking it: locus, taxon, sites masked, number of comparators, number of comparators failed, median residual, whether the cell was ambiguous (fewer than five comparators), and the closest comparator with its alignment-wide identity, overlap, locus identity, expected identity and residual
4. `plots.pdf`, the diagnostic figures (see below)
5. The full results behind the masking:
   - `cells.tsv`: every tested cell, with the same columns as `masked.csv` plus the flag, so unmasked cells can be inspected too
   - `pairs.tsv`: every individual comparison (locus, taxon, comparator rank, comparator, alignment-wide identity, overlap, locus identity, expected identity, residual, pass/fail)
   - `per_locus.tsv`: per locus, the number of sites, taxa, tested cells and pairs, the locus rate, the robust spread of residuals, the band, and the counts flagged by number of comparators failed
   - `trends.tsv`: the locus trend on a grid, as plotted
   - `ani.npy` and `taxa.txt`: the taxon × taxon alignment-wide identity matrix (%) and its row order

# Method overview

The basic idea of the method is to compare a genome to its closest relatives in the alignment, and then look at the distribution of relatedness for each locus in turn. An outlier sequence is a sequence in the focal genome that is far _more_ or _less_ related than expected based on this comparison. 

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
each widened until it holds 40 pairs, joined by linear interpolation. Note that the trend is binned on log divergence (not linear divergence) because linear bins flatten above 99% identity, and in many datasets closest taxa in the alignment are often >99% similar.

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

# Usage

## Installation

The simplest way to run this script is to use Conda as follows. 

[claude: write a conda env file and conda instructions here: insall conda; install env; start env;]

If you prefer to do it another way, the requirements for the script are in the `environment.yml` file.

## Running the script

If you have a nexus file with a SETS block (example filename here is `alignment.nex`): 

    python locus_outliers.py --aln alignment.nex --out scan

If you have a nexus file with no SETS block, or a FASTA or PHYLIP alignment, and a partitions file:

    python locus_outliers.py --aln alignment.fasta --partitions partitions.txt --out scan

[claude - fold the call to the plots right into the main script please]

Method options:

[claude - update these please]

- `--n-comp 5`: comparators per cell, the closest taxa by alignment-wide identity. Must be odd and at least 5.
- `--fit-overlap 20`: sites a pair must share in a locus to be a valid comparison for the trend fit.
- `--seed 1`: random seed for breaking ties between equally close comparators.
