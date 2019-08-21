#!/applications/R/R-3.5.0/bin/Rscript

# Separately for arms and pericentromeres, find genes overlapping HTA7 and H3K9me2 peaks
# and report gene IDs for analysis of GO term enrichment and
# meiocyte-specifically expressed gene enrichment

# Usage:
# ./find_overlapping_genes.R arm

args <- commandArgs(trailingOnly = T)
region <- args[1]
 
library(GenomicRanges)

# Genomic definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- makeGRangesFromDataFrame(data.frame(seqnames = chrs,
                                              start = chrStart,
                                              end = chrLens))
# Conditional definition of regions to be masked
if(region == "arm") {
  maskGR <- makeGRangesFromDataFrame(data.frame(seqnames = chrs,
                                                start = pericenStart,
                                                end = pericenEnd))
} else if(region == "peri") {
  maskGR <- makeGRangesFromDataFrame(data.frame(seqnames = rep(chrs, 2),
                                                start = c(chrStart, pericenEnd+1),
                                                end = c(pericenStart-1, chrLens)))
} else {
  maskGR <- GRanges()
}


# Load HTA7 peaks
load(paste0("../HTA7_ChIP_SRR5298546_rangerPeaksGR_", region,
            "_mergedOverlaps_noMinWidth.RData"))
if(region == "arm") {
  HTA7peaksGR <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  HTA7peaksGR <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
}

load(paste0("/projects/ajt200/BAM_masters/H3K9me2/WT/peaks/PeakRanger1.18/ranger/REC8_MYC_Rep1_input_p0.05_q0.05/WT_H3K9me2_ChIP_rangerPeaksGR_",
            region, "_mergedOverlaps_noMinWidth.RData"))
if(region == "arm") {
  H3K9me2peaksGR <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  H3K9me2peaksGR <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
}

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

# Find genes overlapping HTA7 and H3K9me2 peaks
HTA7_genes_overlap <- findOverlaps(query = HTA7peaksGR,
                                   subject = genesGR,
                                   ignore.strand = TRUE,
                                   select = "all")
HTA7_genesGR <- genesGR[unique(subjectHits(HTA7_genes_overlap))]

H3K9me2_genes_overlap <- findOverlaps(query = H3K9me2peaksGR,
                                      subject = genesGR,
                                      ignore.strand = TRUE,
                                      select = "all")
H3K9me2_genesGR <- genesGR[unique(subjectHits(H3K9me2_genes_overlap))]

HTA7_H3K9me2_genesGR <- HTA7_genesGR[HTA7_genesGR$geneID %in% H3K9me2_genesGR$geneID]

write.table(HTA7_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_peaks_geneIDs.txt"))
write.table(H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_H3K9me2_peaks_geneIDs.txt"))
write.table(HTA7_H3K9me2_genesGR$geneID,
            file = paste0("./genes_overlapping_", region, "_HTA7_and_H3K9me2_peaks_geneIDs.txt"))
