#!/bin/bash

## Use peakranger ranger (v1.18) to call narrow peaks in each REC8 replicate
## PeakRanger manual at http://ranger.sourceforge.net/manual1.18.html

# Usage:
# csmit -m 10G -c 48 "bash ./1_peaks_peakranger_ranger.sh '/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol' '/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol' HTA6_ChIP_SRR5298545 HTA_input_SRR5298544 0.05 0.05 150 48"

ChIP_bamDir=$1
input_bamDir=$2
ChIP_prefix=$3
input_prefix=$4
pval=$5
qval=$6
ext_length=$7
threads=$8

/home/ajt200/tools/PeakRanger-1.18/bin/peakranger ranger \
  --data ${ChIP_bamDir}/${ChIP_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --control ${input_bamDir}/${input_prefix}_RmDup_k10_bt2_mapped_lowmiss_unique_both_sort.bam \
  --format bam \
  --output ${ChIP_prefix}"_rangerPeaks" \
  --pval ${pval} \
  --FDR ${qval} \
  --ext_length ${ext_length} \
  --pad \
  --thread ${threads} \
  --verbose
