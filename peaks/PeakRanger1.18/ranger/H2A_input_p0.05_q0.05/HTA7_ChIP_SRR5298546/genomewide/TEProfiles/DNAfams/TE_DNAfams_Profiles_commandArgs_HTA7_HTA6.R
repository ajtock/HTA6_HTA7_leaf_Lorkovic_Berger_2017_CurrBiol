#!/applications/R/R-4.0.0/bin/Rscript

# Profile mean coverage around DNA TEs overlapped by HTA7_HTA6 peaks and random loci

# Usage via Condor submission system on node7:
#csmit -m 20G -c 1 "/applications/R/R-4.0.0/bin/Rscript ./TE_DNAfams_Profiles_commandArgs_HTA7_HTA6.R 2000 2kb 20 /projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed SPO11_1_oligos_RPI1 HTA7_HTA6 genomewide"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

#flankSize <- 2000
#flankName <- "2kb"
#winSize <- 20
#covDatPath <- "/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/noZscore/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab_noZscore.bed"
#libName <- "SPO11_1_oligos_RPI1"
#peaks <- "HTA7_HTA6"
#region <- "genomewide"

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])
peaks <- as.character(args[6])
region <- as.character(args[7])

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

#peakName1 <- "HTA7_overlapping_HTA6_genomewide_peaks"

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
#mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))
#mask <- toGRanges(data.frame(rep(chrs, 2),
#                             c(chrStart, pericenEnd),
#                             c(pericenStart, chrLens)))

# Test peaks
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData")
HTA7 <- sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))
rangerPeaksGR_arm_mergedOverlaps <- NULL
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(HTA7) <- "*"
print("***********HTA7 peaks***********")
print(HTA7)

load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA6_ChIP_SRR5298545_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA6_ChIP_SRR5298545_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData")
HTA6 <- sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))
rangerPeaksGR_arm_mergedOverlaps <- NULL
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(HTA6) <- "*"
print("***********HTA6 peaks***********")
print(HTA6)

HTA7_HTA6_overlaps <- findOverlaps(query = HTA6,
                                   subject = HTA7,
                                   ignore.strand = TRUE,
                                   select = "all")
peaksGR <- HTA7[unique(subjectHits(HTA7_HTA6_overlaps))]

DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

# Import DNA TEs as GRanges object
for(h in seq_along(DNAfamNames)) {
  # Import TEs as GRanges object
  TEs <- system(paste0("ls ", DNAdir, "TAIR10_Buisine_TEs_strand_tab_ann_", DNAfamNames[h], ".txt"),
                intern = T)
  TEsGR <- readGeneric(TEs,
                       header = T,
                       strand = 4,
                       meta.col = list(Transposon_Name = 5))
  # Get TEs that overlap peaks
  TEs_peaks_overlaps <- findOverlaps(query = peaksGR,
                                     subject = TEsGR,
                                     ignore.strand = TRUE,
                                     select = "all")
  TEsGR <- TEsGR[unique(subjectHits(TEs_peaks_overlaps))]
#  # Get TEs that do not overlap mask
#  TEs_mask_overlaps <- findOverlaps(query = mask,
#                                    subject = TEsGR,
#                                    ignore.strand = TRUE,
#                                    select = "all")
#  TEsGR <- TEsGR[-subjectHits(TEs_mask_overlaps)]
  print("***********TEs***********")
  print(TEsGR)
  # Generate GRanges object containing random loci of same number and
  # size distribution as TEsGR
  set.seed(374592)
  ranLocGR <- randomizeRegions(TEsGR,
                               genome = genome,
#                               mask = mask,
                               per.chromosome = TRUE,
                               allow.overlaps = TRUE)
  
  # Specify locations of normalised per base coverage files
  libPath <- system(paste0("ls ", covDatPath), intern = T)
  
  # Import coverage files as GRanges objects and assign to library names
  covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
  assign(paste0(libName), covGR)
  
  # Define matrix and column mean coverage outfile (mean profiles)
  outDF <- list(paste0(matDir, libName,
                       "_norm_cov_", DNAfamNames[h], "_TEs_overlapping_", peaks, "_", region, "_peaks_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
                paste0(matDir, libName,
                       "_norm_cov_", DNAfamNames[h], "_TEs_overlapping_", peaks, "_", region, "_peaks_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
  outDFcolMeans <- list(paste0(matDir, libName,
                               "_norm_cov_", DNAfamNames[h], "_TEs_overlapping_", peaks, "_", region, "_peaks_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                        paste0(matDir, libName,
                               "_norm_cov_", DNAfamNames[h], "_TEs_overlapping_", peaks, "_", region, "_peaks_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))
  
  # Run covMatrix() function on each coverage GRanges object to obtain matrices
  ## containing normalised coverage values around target and random loci
  covMatrix(signal = covGR,
            feature = TEsGR,
            ranLoc = ranLocGR,
            featureSize = mean(width(TEsGR)),
            flankSize = flankSize,
            winSize = winSize,
            outDF = outDF,
            outDFcolMeans = outDFcolMeans)
  print(paste0(libName, " ", DNAfamNames[h], " profile calculation complete"))
}
