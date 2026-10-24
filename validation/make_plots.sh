#!/bin/bash
# Run from a directory that holds one locus_outliers.py output folder per dataset, named as the dataset
# (and test_data/<dataset>/alignment.nex in the repo for the matches page). Writes one PDF per dataset into
# the repo's plots/datasets (cleared first) and keeps the page PNGs under pngs/ here for inspection.
set -e
REPO=$(cd "$(dirname "$0")/.." && pwd)
KP=$REPO/plots/datasets
rm -f "$KP"/*.pdf; mkdir -p "$KP" pngs
for d in $(ls -d */ | tr -d /); do
  [ -f $d/cells.tsv ] || continue
  n=$(($(wc -l < $d/per_locus.tsv) - 1))
  aln=$REPO/test_data/$d/alignment.nex; [ -f $aln ] && alnopt="--alignment $aln" || alnopt=""
  python $REPO/plot_panels.py $d --png --out $d/plots.pdf $alnopt > /dev/null 2>&1 \
    && cp $d/plots.pdf "$KP/${d}.pdf" && cp $d/plots_*.png pngs/ 2>/dev/null \
    && echo "$d: $n loci" || echo "$d: PLOT FAILED"
done
ls "$KP" | wc -l
