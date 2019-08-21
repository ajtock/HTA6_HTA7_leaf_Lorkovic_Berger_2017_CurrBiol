#!/applications/R/R-3.3.2/bin/Rscript

# Separate arm and pericentromeric peaks
# Convert to GRanges objects
# (with overlapping peaks unmerged or merged)
# sort by decreasing -log10(qval) for use with weeder2

# Usage:
# ./5_rangerPeaks_arm_peri_TreadsNormCreads_sort_by_minuslog10Q_peakSummits.R HTA7_ChIP_SRR5298546

args <- commandArgs(trailingOnly = TRUE)
ChIP_prefix <- args[1]

library(GenomicRanges)
library(dplyr)

# Definitions
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

rangerPeaks <- read.table(paste0("./", ChIP_prefix,
                                 "_rangerPeaks_TreadsNormCreads.narrowPeak"))
rangerPeaks <- cbind(rangerPeaks[,1:3],
                     rangerPeaks[,7:10])
colnames(rangerPeaks) <- c("chr", "start0based", "end",
                           "sigval", "pval", "qval", "summit0based")
rangerPeaks <- as.data.frame(cbind(rangerPeaks$chr,
                                   rangerPeaks$start0based+1,
                                   rangerPeaks$end,
                                   rangerPeaks$sigval,
                                   rangerPeaks$pval,
                                   rangerPeaks$qval,
                                   rangerPeaks$summit0based))
colnames(rangerPeaks) <- c("chr", "start", "end",
                           "sigval", "pval", "qval", "summit0based")

rangerPeaks_arm <- NULL
rangerPeaks_peri <- NULL
for(i in 1:5) {
  chr.rangerPeaks <- rangerPeaks[rangerPeaks$chr == i,]
  chr.rangerPeaks_arm <- chr.rangerPeaks %>%
    filter(end < pericenStart[i] | start > pericenEnd[i])
  chr.rangerPeaks_peri <- chr.rangerPeaks %>%
    filter(end > pericenStart[i] & start < pericenEnd[i])
  rangerPeaks_arm <- rbind(rangerPeaks_arm,
                           chr.rangerPeaks_arm)
  rangerPeaks_peri <- rbind(rangerPeaks_peri,
                           chr.rangerPeaks_peri)
}
print("Unfiltered arm peaks:")
print(dim(rangerPeaks_arm)[[1]])
print("Unfiltered pericentromeric and centromeric peaks:")
print(dim(rangerPeaks_peri)[[1]])
print("Unfiltered peaks:")
print(dim(rangerPeaks_arm)[[1]] + dim(rangerPeaks_peri)[[1]])

write.table(rangerPeaks_arm,
            file = paste0("./", ChIP_prefix,
                          "_rangerPeaks_arm.txt"))
write.table(rangerPeaks_peri,
            file = paste0("./", ChIP_prefix,
                          "_rangerPeaks_peri.txt"))

# arm
### Create GRanges objects sorted by decreasing -log10(qval)
rangerPeaksGR_arm <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks_arm$chr),
                                  ranges = IRanges(start = rangerPeaks_arm$start,
                                                   end = rangerPeaks_arm$end),
                                  strand = "*",
                                  sigval = rangerPeaks_arm$sigval,
                                  pval = rangerPeaks_arm$pval,
                                  qval = rangerPeaks_arm$qval,
                                  summit0based = rangerPeaks_arm$summit0based),
                          by = ~ qval, decreasing = T)
save(rangerPeaksGR_arm,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_minuslog10Qsorted_noMinWidth.RData"))
# Merge overlapping peaks
rangerPeaksGR_arm_mergedOverlaps <- reduce(rangerPeaksGR_arm)
save(rangerPeaksGR_arm_mergedOverlaps,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
# Remove peaks < 50 bp
rangerPeaksGR_arm_50bpMinWidth <- rangerPeaksGR_arm[width(rangerPeaksGR_arm) >= 50]
save(rangerPeaksGR_arm_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_minuslog10Qsorted_50bpMinWidth.RData"))
rangerPeaksGR_arm_mergedOverlaps_50bpMinWidth <- rangerPeaksGR_arm_mergedOverlaps[width(rangerPeaksGR_arm_mergedOverlaps) >= 50]
save(rangerPeaksGR_arm_mergedOverlaps_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_arm_mergedOverlaps_50bpMinWidth.RData"))

# peri
### Create GRanges objects sorted by decreasing -log10(qval)
rangerPeaksGR_peri <- sort(GRanges(seqnames = paste0("Chr", rangerPeaks_peri$chr),
                                   ranges = IRanges(start = rangerPeaks_peri$start,
                                                    end = rangerPeaks_peri$end),
                                   strand = "*",
                                   sigval = rangerPeaks_peri$sigval,
                                   pval = rangerPeaks_peri$pval,
                                   qval = rangerPeaks_peri$qval,
                                   summit0based = rangerPeaks_peri$summit0based),
                           by = ~ qval, decreasing = T)
save(rangerPeaksGR_peri,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_minuslog10Qsorted_noMinWidth.RData"))
# Merge overlapping peaks
rangerPeaksGR_peri_mergedOverlaps <- reduce(rangerPeaksGR_peri)
save(rangerPeaksGR_peri_mergedOverlaps,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"))
# Remove peaks < 50 bp
rangerPeaksGR_peri_50bpMinWidth <- rangerPeaksGR_peri[width(rangerPeaksGR_peri) >= 50]
save(rangerPeaksGR_peri_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_minuslog10Qsorted_50bpMinWidth.RData"))
rangerPeaksGR_peri_mergedOverlaps_50bpMinWidth <- rangerPeaksGR_peri_mergedOverlaps[width(rangerPeaksGR_peri_mergedOverlaps) >= 50]
save(rangerPeaksGR_peri_mergedOverlaps_50bpMinWidth,
     file = paste0("./", ChIP_prefix,
                   "_rangerPeaksGR_peri_mergedOverlaps_50bpMinWidth.RData"))

