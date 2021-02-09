#!/applications/R/R-4.0.0/bin/Rscript

# Get number of TEs in each superfamily overlapping HTA6 and HTA7 peaks

# /applications/R/R-4.0.0/bin/Rscript TEs_overlapping_HTA6_HTA7_peaks.R genomewide 

args <- commandArgs(trailingOnly = T)
region <- args[1]

options(stringsAsFactors = F)
library(regioneR)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
if(region == "peri") {
  mask <- toGRanges(data.frame(rep(chrs, 2),
                               c(chrStart, pericenEnd),
                               c(pericenStart, chrLens)))
} else if(region == "arm") {
  mask <- toGRanges(data.frame(chrs,
                               pericenStart,
                               pericenEnd))
} else if(region == "genomewide") {
  mask <- GRanges()
}

# Test peaks
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA6_ChIP_SRR5298545_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA6_ChIP_SRR5298545_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData")
HTA6 <- sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))
rangerPeaksGR_arm_mergedOverlaps <- NULL
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(HTA6) <- "*"
print("***********HTA6 peaks***********")
print(HTA6)

load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData")
load("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData")
HTA7 <- sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))
rangerPeaksGR_arm_mergedOverlaps <- NULL
rangerPeaksGR_peri_mergedOverlaps <- NULL
strand(HTA7) <- "*"
print("***********HTA7 peaks***********")
print(HTA7)

HTA6_HTA7_overlaps <- findOverlaps(query = HTA7,
                                   subject = HTA6,
                                   ignore.strand = TRUE,
                                   select = "all")
HTA6_HTA7 <- HTA6[unique(subjectHits(HTA6_HTA7_overlaps))]
HTA6_not_HTA7 <- HTA6[-subjectHits(HTA6_HTA7_overlaps)]

HTA7_HTA6_overlaps <- findOverlaps(query = HTA6,
                                   subject = HTA7,
                                   ignore.strand = TRUE,
                                   select = "all")
HTA7_HTA6 <- HTA7[unique(subjectHits(HTA7_HTA6_overlaps))]
HTA7_not_HTA6 <- HTA7[-subjectHits(HTA7_HTA6_overlaps)]

# TEs
DNAfamNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger")
RNAfamNames <- c("rna", "gypsy", "copia", "linel1", "sine")
DNAfamNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/Mariner", "MuDR", "EnSpm/CACTA", "hAT", "Harbinger")
RNAfamNamesPlot <- c("RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")
DNAdir <- "/projects/ajt200/TAIR10/TE_classes/DNA/"
RNAdir <- "/projects/ajt200/TAIR10/TE_classes/RNA/"

### DNA TEs
TEsDNAGR <- lapply(seq_along(DNAfamNames), function(x) {
  TEsDNA <- read.table(paste0(DNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              DNAfamNames[x], ".txt"),
                       header = T)
  TEsDNAGR_all <- GRanges(seqnames = TEsDNA$chr,
                          ranges = IRanges(start = TEsDNA$start,
                                           end = TEsDNA$end),
                          strand = "*")
  maskTEsDNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsDNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  if(length(subjectHits(maskTEsDNAoverlaps)) > 0) {
    TEsDNAGR_all[-subjectHits(maskTEsDNAoverlaps)]
  } else {
    TEsDNAGR_all
  }
})
### RNA TEs
TEsRNAGR <- lapply(seq_along(RNAfamNames), function(x) {
  TEsRNA <- read.table(paste0(RNAdir,
                              "TAIR10_Buisine_TEs_strand_tab_ann_",
                              RNAfamNames[x], ".txt"),
                       header = T)
  TEsRNAGR_all <- GRanges(seqnames = TEsRNA$chr,
                          ranges = IRanges(start = TEsRNA$start,
                                           end = TEsRNA$end),
                          strand = "*")
  maskTEsRNAoverlaps <- findOverlaps(query = mask,
                                     subject = TEsRNAGR_all,
                                     ignore.strand = TRUE,
                                     select = "all")
  if(length(subjectHits(maskTEsRNAoverlaps)) > 0) {
    TEsRNAGR_all[-subjectHits(maskTEsRNAoverlaps)]
  } else {
    TEsRNAGR_all
  }
})

otherNames <- c(
                DNAfamNames,
                RNAfamNames
               )
otherNamesPlot <- c(
                    DNAfamNamesPlot,
                    RNAfamNamesPlot
                   )
grl <- c(
         TEsDNAGR,
         TEsRNAGR
        )

# Count number of TEs in each superfamily overlapping given HTA6 or HTA7 peak sets
HTA6_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA6,
                    ignore.strand = TRUE) > 0)
})
HTA6_HTA7_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA6_HTA7,
                    ignore.strand = TRUE) > 0)
})
HTA6_not_HTA7_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA6_not_HTA7,
                    ignore.strand = TRUE) > 0)
})
HTA7_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA7,
                    ignore.strand = TRUE) > 0)
})
HTA7_HTA6_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA7_HTA6,
                    ignore.strand = TRUE) > 0)
})
HTA7_not_HTA6_TEs_overlaps <- sapply(seq_along(grl), function(x) {
  sum(countOverlaps(query = grl[[x]],
                    subject = HTA7_not_HTA6,
                    ignore.strand = TRUE) > 0)
})

TEs_overlapping_HTA6_HTA7 <- data.frame(feature = otherNames,
                                        total = sapply(grl, function(x) length(x)),
                                        featuresOverlapping_HTA6 = HTA6_TEs_overlaps,
                                        featuresOverlapping_HTA6_HTA7 = HTA6_HTA7_TEs_overlaps,
                                        featuresOverlapping_HTA6_not_HTA7 = HTA6_not_HTA7_TEs_overlaps,
                                        featuresOverlapping_HTA7 = HTA7_TEs_overlaps,
                                        featuresOverlapping_HTA7_HTA6 = HTA7_HTA6_TEs_overlaps,
                                        featuresOverlapping_HTA7_not_HTA6 = HTA7_not_HTA6_TEs_overlaps)
write.table(TEs_overlapping_HTA6_HTA7,
            file = paste0("TEs_overlapping_HTA6_and_HTA7_peak_sets_", region , ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
