#!/applications/R/R-3.5.0/bin/Rscript

# Plot bar graph of log2(observed:expected) peaks overlapping other features

# Usage:
# ./TE_families_wt_only_HTA7_overlapping_HTA6_arm_peaks.R "Arm HTA7 peaks overlapping HTA6 peaks" "HTA7_HTA6" HTA7_overlapping_HTA6 HTA7_overlapping_HTA6 arm 10000

library(ggplot2)
library(ggthemes)

#plotTitle <- "Arm HTA7 peaks overlapping HTA6 peaks"
#pt1Name <- ""HTA7_HTA6""
#pt1LibName <- "HTA7_overlapping_HTA6" 
#ptOrderLibName <- "HTA7_overlapping_HTA6"
#region <- "arm" 
## Number of permutations (randomisations) performed
#perms <- "10000"

args <- commandArgs(trailingOnly = T)
plotTitle <- args[1]
pt1Name <- args[2]
pt1LibName <- args[3]
ptOrderLibName <- args[4]
region <- args[5]
# Number of permutations (randomisations) performed
perms <- as.character(args[6])

ptOrderDir <- paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546/",
                     "peri/regioneR/noMinWidth_mergedOverlaps/")
pt1Dir <- paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546/",
                 region, "/regioneR/noMinWidth_mergedOverlaps/")

plotDir <- "./bar_graphs/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

famNames <- c("dna", "heli", "ptmari", "mudr", "enspm", "hat", "harbinger",
              "rna", "gypsy", "copia", "linel1", "sine")
famNamesPlot <- c("DNA", "Helitron", "Pogo/Tc1/Mariner", "MuDR", "EnSpm", "hAT", "Harbinger",
                  "RNA", "Gypsy LTR", "Copia LTR", "LINE-1", "SINE")

# Load permutation test results for peak set to be used for ordering
# of other features in bar graph
load(paste0(ptOrderDir, "permTest_", perms, "perms_",
            ptOrderLibName, "_peri_peaks_vs_TEsDNA.RData"))
ptOrder_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL
load(paste0(ptOrderDir, "permTest_", perms, "perms_",
            ptOrderLibName, "_peri_peaks_vs_TEsRNA.RData"))
ptOrder_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL

# Load permutation test results to be used for plotting
# DNA
load(paste0(pt1Dir, "permTest_", perms, "perms_",
            pt1LibName, "_", region, "_peaks_vs_TEsDNA.RData"))
pt1_DNA <- ptPeaksTEsDNAPerChrom
ptPeaksTEsDNAPerChrom <- NULL

# RNA
load(paste0(pt1Dir, "permTest_", perms, "perms_",
            pt1LibName, "_", region, "_peaks_vs_TEsRNA.RData"))
pt1_RNA <- ptPeaksTEsRNAPerChrom
ptPeaksTEsRNAPerChrom <- NULL

# Combined
ptOrder <- c(ptOrder_DNA, ptOrder_RNA)
pt1 <- c(pt1_DNA, pt1_RNA)

# ptOrder
ptOrder_Pval <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$pval
})
ptOrder_Obs <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$observed
})
ptOrder_Perm <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$permuted
})
ptOrder_Exp <- lapply(seq_along(ptOrder), function(x) {
  mean(ptOrder[[x]]$numOverlaps$permuted)
})
ptOrder_log2ObsExp <- lapply(seq_along(ptOrder_Obs), function(x) {
  log2((ptOrder_Obs[[x]]+1)/(ptOrder_Exp[[x]]+1))
})
ptOrder_Zscore <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$zscore
})
ptOrder_AltHyp <- lapply(seq_along(ptOrder), function(x) {
  ptOrder[[x]]$numOverlaps$alternative
})
ptOrder_alpha0.05 <- lapply(seq_along(ptOrder_Perm), function(x) {
  if(ptOrder_AltHyp[[x]] == "greater") {
    quantile(ptOrder_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(ptOrder_Perm[[x]], 0.05)[[1]]
  }
})
ptOrder_log2alpha0.05 <- lapply(seq_along(ptOrder_alpha0.05), function(x) {
  log2((ptOrder_alpha0.05[[x]]+1)/(ptOrder_Exp[[x]]+1))
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
ptOrder_log2ObsExp_sorted <- sort.int(unlist(ptOrder_log2ObsExp),
                                      decreasing = T)
ptOrder_log2alpha0.05_sorted <- unlist(ptOrder_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                                      decreasing = T,
                                                                      index.return = T)$ix])
ptOrder_famNames_sorted <- famNames[sort.int(unlist(ptOrder_log2ObsExp),
                                             decreasing = T,
                                             index.return = T)$ix]
ptOrder_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                     decreasing = T,
                                                     index.return = T)$ix]

# pt1
pt1_Pval <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$pval
})
pt1_Obs <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$observed
})
pt1_Perm <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$permuted
})
pt1_Exp <- lapply(seq_along(pt1), function(x) {
  mean(pt1[[x]]$numOverlaps$permuted)
})
pt1_log2ObsExp <- lapply(seq_along(pt1_Obs), function(x) {
  log2((pt1_Obs[[x]]+1)/(pt1_Exp[[x]]+1))
})
pt1_Zscore <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$zscore
})
pt1_AltHyp <- lapply(seq_along(pt1), function(x) {
  pt1[[x]]$numOverlaps$alternative
})
pt1_alpha0.05 <- lapply(seq_along(pt1_Perm), function(x) {
  if(pt1_AltHyp[[x]] == "greater") {
    quantile(pt1_Perm[[x]], 0.95)[[1]]
  } else {
    quantile(pt1_Perm[[x]], 0.05)[[1]]
  }
})
pt1_log2alpha0.05 <- lapply(seq_along(pt1_alpha0.05), function(x) {
  log2((pt1_alpha0.05[[x]]+1)/(pt1_Exp[[x]]+1))
})

# Order permutation test statistics and feature names by
# decreasing log2(observed/expected) overlaps of wt treatment peaks
pt1_log2ObsExp_sorted <- unlist(pt1_log2ObsExp[sort.int(unlist(ptOrder_log2ObsExp),
                                                        decreasing = T,
                                                        index.return = T)$ix])
pt1_log2alpha0.05_sorted <- unlist(pt1_log2alpha0.05[sort.int(unlist(ptOrder_log2ObsExp),
                                                              decreasing = T,
                                                              index.return = T)$ix])
pt1_famNames_sorted <- famNames[sort.int(unlist(ptOrder_log2ObsExp),
                                         decreasing = T,
                                         index.return = T)$ix]
pt1_famNamesPlot_sorted <- famNamesPlot[sort.int(unlist(ptOrder_log2ObsExp),
                                                 decreasing = T,
                                                 index.return = T)$ix]

# Combine in data.frame
df <- data.frame(Sample = rep(c(pt1Name),
                              each = length(ptOrder_log2ObsExp_sorted)),
                 Transposon_family = rep(pt1_famNamesPlot_sorted, 1),
                 log2ObsExp = c(pt1_log2ObsExp_sorted),
                 log2alpha0.05 = c(pt1_log2alpha0.05_sorted))

df$Transposon_family <- factor(df$Transposon_family,
                               levels = pt1_famNamesPlot_sorted)
df$Sample <- factor(df$Sample,
                    levels = c(pt1Name))

bp <- ggplot(data = df,
             mapping = aes(x = Transposon_family,
                           y = log2ObsExp,
                           fill = Sample)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  scale_fill_manual(name = "",
                    values = c("dodgerblue2"),
                    labels = c(pt1Name)) +
  geom_point(mapping = aes(Transposon_family, log2alpha0.05),
             position = position_dodge(0.9),
             shape = "-", colour  = "grey70", size = 10) +
  labs(y = expression("Log"[2]*"(observed/expected) peak overlap")) +
  scale_y_continuous(limits = c(-4.0, 4.0)) +
  scale_x_discrete(position = "top") +
  guides(fill = guide_legend(direction = "horizontal",
                             label.position = "top",
                             label.theme = element_text(size = 20, hjust = 0, vjust = 0.5, angle = 90),
                             nrow = 1,
                             byrow = TRUE)) +
  theme_bw() +
  theme(axis.line.y = element_line(size = 1, colour = "black"),
        axis.ticks.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", hjust = 0.5, vjust = 0.5, angle = 90),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        #legend.position = c(0.05, 0.30),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 10.5, 5.5), "pt"),
        plot.title = element_text(size = 17, colour = "black", hjust = 0.5)) +
  ggtitle(paste0(plotTitle, " (", prettyNum(as.character(perms),
                                            big.mark = ",", trim = "T"),
                 " permutations)"))
ggsave(paste0(plotDir, "barplot_TE_families_permTestResults_",
              as.character(perms), "perms_",
              "log2_Observed_Expected_",
              pt1LibName, "_",
              region, "_peaks.pdf"),
       plot = bp,
       height = 8, width = 12)
save(bp,
     file = paste0(plotDir, "barplot_TE_families_permTestResults_",
                   as.character(perms), "perms_",
                   "log2_Observed_Expected_",
                   pt1LibName, "_",
                   region, "_peaks.RData"))


