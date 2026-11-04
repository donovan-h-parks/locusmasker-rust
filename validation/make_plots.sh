#!/bin/bash
# Run from a directory that holds one locusmasker.py output folder per dataset, named as the dataset
# (and test_data/<dataset>/alignment.nex in the repo for the matches page). A folder counts when it holds a run's
# masked.csv and report.npz (a run made with --no-report has nothing to plot). Writes one PDF per dataset into
# the repo's plots/datasets (cleared first; PLOTS_DIR=dir writes there instead) and keeps the page PNGs under pngs/
# here for inspection.
set -e
REPO=$(cd "$(dirname "$0")/.." && pwd)
KP=${PLOTS_DIR:-$REPO/plots/datasets}
rm -f "$KP"/*.pdf; mkdir -p "$KP" pngs
for d in $(ls -d */ | tr -d /); do
  [ -f $d/masked.csv ] && [ -f $d/report.npz ] || continue
  n=$(python -c 'import json, sys, zipfile; print(len(json.loads(zipfile.ZipFile(sys.argv[1]).read("loci.json"))["names"]))' $d/report.npz 2>/dev/null) || n="?"
  aln=$REPO/test_data/$d/alignment.nex; [ -f $aln ] && alnopt="--alignment $aln" || alnopt=""
  python $REPO/plot_panels.py $d --png --out $d/plots.pdf $alnopt > /dev/null 2>&1 \
    && cp $d/plots.pdf "$KP/${d}.pdf" && cp $d/plots_*.png pngs/ 2>/dev/null \
    && echo "$d: $n loci" || echo "$d: PLOT FAILED"
done
ls "$KP" | wc -l
