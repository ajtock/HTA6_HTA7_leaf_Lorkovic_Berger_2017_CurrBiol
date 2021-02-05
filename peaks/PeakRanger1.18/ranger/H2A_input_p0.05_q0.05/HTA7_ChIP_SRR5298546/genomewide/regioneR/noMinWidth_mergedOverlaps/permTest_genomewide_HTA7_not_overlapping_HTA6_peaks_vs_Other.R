#!/applications/R/R-4.0.0/bin/Rscript

# Use permutation test function in regioneR to determine
# if peaks overlap features of interest more or less
# than expected by chance

# Usage on hydrogen node7:
# csmit -m 100G -c 32 "/applications/R/R-4.0.0/bin/Rscript permTest_genomewide_HTA7_not_overlapping_HTA6_peaks_vs_Other.R 10000"

args <- commandArgs(trailingOnly = T)
perms <- as.numeric(args[1])

peakName1 <- "HTA7_not_overlapping_HTA6_genomewide_peaks"

outDir <- "./"

library(regioneR)

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

plotDir <- "./histograms/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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
peaksGR <- HTA7[-subjectHits(HTA7_HTA6_overlaps)]

### Other
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/armrangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
load("/projects/ajt200/BAM_masters/SPO11-oligo/WT/peaks/PeakRanger1.18/ranger/p0.2_q0.2/idr/qValRank/minuslog10_pval_qval/perirangerPeaksGRmerge_RPI1_RPI8_idr0.05_noMinWidth.RData")
SPO11GR <- sort(c(armrangerPeaksGRmerge, perirangerPeaksGRmerge))
armrangerPeaksGRmerge <- NULL
perirangerPeaksGRmerge <- NULL
strand(SPO11GR) <- "*"
print(length(SPO11GR))
#[1] 5914

load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/armPeaksSH99GRmerge.RData")
load("/projects/ajt200/REC8_MSH4/nuc_peaks/log2ChIPinput/nucleR/trim/analysis_01/periPeaksSH99GRmerge.RData")
nucleRnucsGR <- sort(c(armPeaksSH99GRmerge, periPeaksSH99GRmerge))
armPeaksSH99GRmerge <- NULL
periPeaksSH99GRmerge <- NULL
strand(nucleRnucsGR) <- "*"
print(length(nucleRnucsGR))
#[1] 57734

OtherNames <- c(
                "SPO11GR",
                "nucleRnucsGR"
               )
OtherGR <- c(
             "SPO11GR" = SPO11GR,
             "nucleRnucsGR" = nucleRnucsGR
            )

# Perform permutation tests with randomized regions generated on a per chromosome basis;
# same per-chromosome number and size of regions as in A
set.seed(38402)
ptPeaksOtherPerChrom <- lapply(seq_along(OtherGR), function(x) {
  permTest(A = peaksGR,
           B = OtherGR[[x]],
           genome = genome,
#           mask = mask,
           randomize.function = randomizeRegions,
           allow.overlaps = TRUE,
           per.chromosome = TRUE,
           evaluate.function = numOverlaps,
           count.once = TRUE,
           ntimes = perms,
           mc.set.seed = FALSE,
           mc.cores = detectCores())
})

for(i in 1:length(ptPeaksOtherPerChrom)) {
  assign(paste0(OtherNames[i]), ptPeaksOtherPerChrom[[i]])
}
save(ptPeaksOtherPerChrom,
     file = paste0(outDir,
                   "permTest_", as.character(perms), "perms_", peakName1, "_vs_Other.RData"))

# Summarise results in a table
featureName <- NULL
noOfFeatures <- NULL
expected <- NULL
observed <- NULL
pval <- NULL
zscore <- NULL
for(i in 1:length(ptPeaksOtherPerChrom)) {
  featureNamei <- print(OtherNames[i])
  featureName <- c(featureName, featureNamei)
  noOfFeaturesi <- print(length(OtherGR[[i]]))
  noOfFeatures <- c(noOfFeatures, noOfFeaturesi)
  expectedi <- print(round(mean(ptPeaksOtherPerChrom[[i]]$numOverlaps$permuted)))
  expected <- c(expected, expectedi)
  observedi <- print(ptPeaksOtherPerChrom[[i]]$numOverlaps$observed)
  observed <- c(observed, observedi)
  pvali <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$pval, 4))
  pval <- c(pval, pvali)
  zscorei <- print(round(ptPeaksOtherPerChrom[[i]]$numOverlaps$zscore, 4))
  zscore <- c(zscore, zscorei)
}
ptPeaksOtherPerChromDataFrame <- cbind(featureName, noOfFeatures, expected, observed, pval, zscore)
colnames(ptPeaksOtherPerChromDataFrame) <- c("feature", "n", "expected", "observed", "pval", "zscore")
write.table(ptPeaksOtherPerChromDataFrame,
            file = paste0(outDir,
                          "permTest_", as.character(perms), "perms_", peakName1, "_vs_Other_DataFrame.txt"),
            sep = "\t", row.names = F)

# plot graphical summaries of results
for(i in 1:length(ptPeaksOtherPerChrom)) {
  pdf(file = paste0(plotDir, OtherNames[i], "_permTest_", as.character(perms), "perms_", peakName1, "_perChrom.pdf"), width = 10, height = 7)
  plot(ptPeaksOtherPerChrom[[i]], main = paste0(peakName1, " vs ", OtherNames[i]), xlab = "Number of overlaps", ylab = "Density")
  dev.off()

  # Using the localZScore() function, evaluate whether the association between peaks and Other is highly dependent on their exact position
  lz_1kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = OtherGR[[i]],
                        window = 1000, step = 50, count.once = TRUE)
  lz_10kb <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = OtherGR[[i]],
                         window = 10000, step = 500, count.once = TRUE)
  lz_custom <- localZScore(pt = ptPeaksOtherPerChrom[[i]], A = peaksGR, B = OtherGR[[i]],
                           window = 10*mean(width(peaksGR)), step = mean(width(peaksGR))/2, count.once = TRUE)
  win <- as.character(round((10*mean(width(peaksGR)))/1000))
  step <- as.character(round(mean(width(peaksGR))/2))
  pdf(file = paste0(plotDir, OtherNames[i], "_localZscore_permTest_", as.character(perms), "perms_", peakName1, "_w1kb_s50bp_w10kb_s500bp_w", win ,"kb_s", step, "bp_perChrom.pdf"))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(lz_1kb, main = paste0(peakName1, " vs ", OtherNames[i], " (1-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(peakName1, " vs ", OtherNames[i], " (1-kb shift)"))
  plot(lz_10kb, main = paste0(peakName1, " vs ", OtherNames[i], " (10-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(peakName1, " vs ", OtherNames[i], " (10-kb shift)"))
  plot(lz_custom, main = paste0(peakName1, " vs ", OtherNames[i], " (~", win, "-kb shift)"))
  mtext(side = 3, at = 2, text = paste0(peakName1, " vs ", OtherNames[i], " (~", win, "-kb shift)"))
  dev.off()
}
