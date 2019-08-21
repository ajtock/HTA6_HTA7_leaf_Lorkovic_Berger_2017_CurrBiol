#!/bin/bash

# Convert PeakRanger peak loci details file to narrowPeak format:
# https://genome.ucsc.edu/FAQ/FAQformat.html

# Usage:
# ./2_peaks_detailsTOnarrowPeak_Treads_Creads.sh HTA7_ChIP_SRR5298546

ChIP_prefix=$1

# Output line 24 to last line of file (skipping commented-out lines)
tail -n +24 \
  ${ChIP_prefix}_rangerPeaks_details > ${ChIP_prefix}_rangerPeaks_noHead
# Exclude peaks in chloroplast and mitochondria
grep -v 'chloroplast\|mitochondria' \
  ${ChIP_prefix}_rangerPeaks_noHead > ${ChIP_prefix}_rangerPeaks_noHead_tmp1
# Add "chr" to chromosome names in column 1
awk 'BEGIN {OFS="\t"}; {$1 = "chr"$1; print}' \
  ${ChIP_prefix}_rangerPeaks_noHead_tmp1  > ${ChIP_prefix}_rangerPeaks_noHead_tmp2
# Output chromosome name (column 1), peak start and end coordinates (columns 2 & 3),
# strand (column 8), P-value and FDR (columns 6 and 7), peak summit (column 5),
# peak ChIP reads and control reads (columns 9 and 10)
awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $8, $8, $8, $8, $6, $7, $5, $9, $10}' ${ChIP_prefix}_rangerPeaks_noHead_tmp2 > ${ChIP_prefix}_rangerPeaks_noHead_tmp3
# As strand is not known, replace strand information with "."
# Make peak summit relative to peak start coordinate (0-based)
# Note that the output file contains untransformed P-values and FDR values,
# unlike standard narrowPeak format files, which contain
# -log10(P) and -log10(FDR) values
awk 'BEGIN {OFS="\t"}; {$4 = "."; $5 = "."; $6 = "."; $7 = "."; $10 = $10-$2; print}' ${ChIP_prefix}_rangerPeaks_noHead_tmp3 > ${ChIP_prefix}_rangerPeaks_Treads_Creads.narrowPeak.UntransformedPQ
# Remove intermediate files
rm ${ChIP_prefix}_rangerPeaks_noHead*
