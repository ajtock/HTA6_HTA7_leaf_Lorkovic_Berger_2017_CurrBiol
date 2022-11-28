#!/applications/R/R-4.0.0/bin/Rscript


# Get genes that overlap HTA7 ChIP-seq peaks and define a matched
# set of randomly positioned loci

# Usage:
# ./get_genes_overlapping_HTA7_peaks.R

args <- commandArgs(trailingOnly = T)
 
library(GenomicRanges)

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Load HTA7 peaks
load(paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_",
            region, "_mergedOverlaps_noMinWidth.RData"))
if(region == "arm") {
  HTA7 <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  HTA7 <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
}

# Load HTA6 peaks
load(paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA6_ChIP_SRR5298545_rangerPeaksGR_",
            region, "_mergedOverlaps_noMinWidth.RData"))
if(region == "arm") {
  HTA6 <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  HTA6 <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
}

# Load H3K9me2 peaks
load(paste0("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep1_input_p0.05_q0.05/WT_H3K9me2_ChIP_rangerPeaksGR_",
            region, "_mergedOverlaps_noMinWidth.RData"))
if(region == "arm") {
  H3K9me2 <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  H3K9me2 <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
}

# Get overlapping and nonoverlapping peaks
HTA7_HTA6_overlaps <- findOverlaps(query = HTA6,
                                   subject = HTA7,
                                   ignore.strand = TRUE,
                                   select = "all")
HTA7_H3K9me2_overlaps <- findOverlaps(query = H3K9me2,
                                      subject = HTA7,
                                      ignore.strand = TRUE,
                                      select = "all")
HTA6_H3K9me2_overlaps <- findOverlaps(query = H3K9me2,
                                      subject = HTA6,
                                      ignore.strand = TRUE,
                                      select = "all")
# overlapping
HTA7_HTA6 <- HTA7[unique(subjectHits(HTA7_HTA6_overlaps))]
HTA6_HTA7 <- HTA6[unique(queryHits(HTA7_HTA6_overlaps))]
HTA7_H3K9me2 <- HTA7[unique(subjectHits(HTA7_H3K9me2_overlaps))]
H3K9me2_HTA7 <- H3K9me2[unique(queryHits(HTA7_H3K9me2_overlaps))]
HTA6_H3K9me2 <- HTA6[unique(subjectHits(HTA6_H3K9me2_overlaps))]
H3K9me2_HTA6 <- H3K9me2[unique(queryHits(HTA6_H3K9me2_overlaps))]
# nonoverlapping
HTA7_not_HTA6 <- HTA7[-subjectHits(HTA7_HTA6_overlaps)]
HTA6_not_HTA7 <- HTA6[-queryHits(HTA7_HTA6_overlaps)]
HTA7_not_H3K9me2 <- HTA7[-subjectHits(HTA7_H3K9me2_overlaps)]
H3K9me2_not_HTA7 <- H3K9me2[-queryHits(HTA7_H3K9me2_overlaps)]
HTA6_not_H3K9me2 <- HTA6[-subjectHits(HTA6_H3K9me2_overlaps)]
H3K9me2_not_HTA6 <- H3K9me2[-queryHits(HTA6_H3K9me2_overlaps)]

stopifnot((length(HTA7_HTA6) + length(HTA7_not_HTA6)) == length(HTA7))
stopifnot((length(HTA6_HTA7) + length(HTA6_not_HTA7)) == length(HTA6))
stopifnot((length(HTA7_H3K9me2) + length(HTA7_not_H3K9me2)) == length(HTA7))
stopifnot((length(HTA6_H3K9me2) + length(HTA6_not_H3K9me2)) == length(HTA6))
stopifnot((length(H3K9me2_HTA7) + length(H3K9me2_not_HTA7)) == length(H3K9me2))
stopifnot((length(H3K9me2_HTA6) + length(H3K9me2_not_HTA6)) == length(H3K9me2))

# Load genes
genes <- read.table("/projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt",
                    header = T)
genesGR <- GRanges(seqnames = genes$chr,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = "*",
                   geneID = sub(pattern = "(AT\\dG\\d+).\\d+",
                                replacement = "\\1",
                                x = genes$gene_model))
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print(length(genesGR))

# For each group of peaks, find overlapping genes
# overlapping peaks
HTA7_HTA6_genes_overlap <- findOverlaps(query = HTA7_HTA6,
                                        subject = genesGR,
                                        ignore.strand = TRUE,
                                        select = "all")
HTA7_HTA6_genesGR <- genesGR[unique(subjectHits(HTA7_HTA6_genes_overlap))]

HTA6_HTA7_genes_overlap <- findOverlaps(query = HTA6_HTA7,
                                        subject = genesGR,
                                        ignore.strand = TRUE,
                                        select = "all")
HTA6_HTA7_genesGR <- genesGR[unique(subjectHits(HTA6_HTA7_genes_overlap))]

HTA7_H3K9me2_genes_overlap <- findOverlaps(query = HTA7_H3K9me2,
                                           subject = genesGR,
                                           ignore.strand = TRUE,
                                           select = "all")
HTA7_H3K9me2_genesGR <- genesGR[unique(subjectHits(HTA7_H3K9me2_genes_overlap))]

H3K9me2_HTA7_genes_overlap <- findOverlaps(query = H3K9me2_HTA7,
                                           subject = genesGR,
                                           ignore.strand = TRUE,
                                           select = "all")
H3K9me2_HTA7_genesGR <- genesGR[unique(subjectHits(H3K9me2_HTA7_genes_overlap))]

HTA6_H3K9me2_genes_overlap <- findOverlaps(query = HTA6_H3K9me2,
                                           subject = genesGR,
                                           ignore.strand = TRUE,
                                           select = "all")
HTA6_H3K9me2_genesGR <- genesGR[unique(subjectHits(HTA6_H3K9me2_genes_overlap))]

H3K9me2_HTA6_genes_overlap <- findOverlaps(query = H3K9me2_HTA6,
                                           subject = genesGR,
                                           ignore.strand = TRUE,
                                           select = "all")
H3K9me2_HTA6_genesGR <- genesGR[unique(subjectHits(H3K9me2_HTA6_genes_overlap))]

# nonoverlapping peaks
HTA7_not_HTA6_genes_overlap <- findOverlaps(query = HTA7_not_HTA6,
                                            subject = genesGR,
                                            ignore.strand = TRUE,
                                            select = "all")
HTA7_not_HTA6_genesGR <- genesGR[unique(subjectHits(HTA7_not_HTA6_genes_overlap))]

HTA6_not_HTA7_genes_overlap <- findOverlaps(query = HTA6_not_HTA7,
                                            subject = genesGR,
                                            ignore.strand = TRUE,
                                            select = "all")
HTA6_not_HTA7_genesGR <- genesGR[unique(subjectHits(HTA6_not_HTA7_genes_overlap))]

HTA7_not_H3K9me2_genes_overlap <- findOverlaps(query = HTA7_not_H3K9me2,
                                               subject = genesGR,
                                               ignore.strand = TRUE,
                                               select = "all")
HTA7_not_H3K9me2_genesGR <- genesGR[unique(subjectHits(HTA7_not_H3K9me2_genes_overlap))]

H3K9me2_not_HTA7_genes_overlap <- findOverlaps(query = H3K9me2_not_HTA7,
                                               subject = genesGR,
                                               ignore.strand = TRUE,
                                               select = "all")
H3K9me2_not_HTA7_genesGR <- genesGR[unique(subjectHits(H3K9me2_not_HTA7_genes_overlap))]

HTA6_not_H3K9me2_genes_overlap <- findOverlaps(query = HTA6_not_H3K9me2,
                                               subject = genesGR,
                                               ignore.strand = TRUE,
                                               select = "all")
HTA6_not_H3K9me2_genesGR <- genesGR[unique(subjectHits(HTA6_not_H3K9me2_genes_overlap))]

H3K9me2_not_HTA6_genes_overlap <- findOverlaps(query = H3K9me2_not_HTA6,
                                               subject = genesGR,
                                               ignore.strand = TRUE,
                                               select = "all")
H3K9me2_not_HTA6_genesGR <- genesGR[unique(subjectHits(H3K9me2_not_HTA6_genes_overlap))]

# Write gene IDs to file
write.table(HTA7_HTA6_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_HTA6_peaks_geneIDs.txt"))
write.table(HTA6_HTA7_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA6_HTA7_peaks_geneIDs.txt"))
write.table(HTA7_H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_H3K9me2_peaks_geneIDs.txt"))
write.table(H3K9me2_HTA7_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_H3K9me2_HTA7_peaks_geneIDs.txt"))
write.table(HTA6_H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA6_H3K9me2_peaks_geneIDs.txt"))
write.table(H3K9me2_HTA6_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_H3K9me2_HTA6_peaks_geneIDs.txt"))

write.table(HTA7_not_HTA6_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_not_HTA6_peaks_geneIDs.txt"))
write.table(HTA6_not_HTA7_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA6_not_HTA7_peaks_geneIDs.txt"))
write.table(HTA7_not_H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_not_H3K9me2_peaks_geneIDs.txt"))
write.table(H3K9me2_not_HTA7_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_H3K9me2_not_HTA7_peaks_geneIDs.txt"))
write.table(HTA6_not_H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA6_not_H3K9me2_peaks_geneIDs.txt"))
write.table(H3K9me2_not_HTA6_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_H3K9me2_not_HTA6_peaks_geneIDs.txt"))
