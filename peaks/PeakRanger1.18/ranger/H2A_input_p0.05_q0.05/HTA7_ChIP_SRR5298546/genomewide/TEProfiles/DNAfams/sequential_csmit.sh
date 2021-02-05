#!/bin/bash

csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI1 HTA7_HTA6 genomewide" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/log2wtNucNakedDNAuntrimmed_norm_allchrs_coverage_coord_tab_noZscore.bed MNase HTA7_HTA6 genomewide" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI3NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI3 HTA7_HTA6 genomewide" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI8NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI8 HTA7_HTA6 genomewide" & sleep 10;
#csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs.R 2000 2kb 20 /projects/ajt200/BAM_masters/H3K9me2/WT/coverage/log2ChIPinput/noZscore/log2_WT_H3K9me2_ChIP_WT_H3K9me2_input_norm_allchrs_coverage_coord_tab_noZscore.bed H3K9me2 HTA7_HTA6 genomewide" & sleep 10;
wait
