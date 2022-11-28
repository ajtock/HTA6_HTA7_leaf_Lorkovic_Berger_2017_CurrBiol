#!/applications/R/R-4.0.0/bin/Rscript

# Get genes overlapping HTA7 peaks and matched random loci

# Usage:
# ./extract_gene_coords_genes_overlapping_HTA7_peaks.R genomewide

#region <- "genomewide"

args <- commandArgs(trailingOnly = T)
region <- args[1]

options(stringsAsFactors = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(data.table)

# Load HTA7 peaks
load(paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
load(paste0("/home/ajt200/analysis/HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/peaks/PeakRanger1.18/ranger/H2A_input_p0.05_q0.05/HTA7_ChIP_SRR5298546_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"))

if(region == "arm") {
  HTA7 <- rangerPeaksGR_arm_mergedOverlaps
  rangerPeaksGR_arm_mergedOverlaps <- NULL
} else if(region == "peri") {
  HTA7 <- rangerPeaksGR_peri_mergedOverlaps
  rangerPeaksGR_peri_mergedOverlaps <- NULL
} else {
  HTA7 <- c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps)
  HTA7 <- sortSeqlevels(HTA7)
  HTA7 <- sort(HTA7)
}
seqlevels(HTA7) <- sub("^Chr", "", seqlevels(HTA7))

# Load gene IDs for genes up-regulated in hta7 mutant
hta7_upreg <- read.table("/home/ajt200/analysis/Pallas_RNAseq_tetoolkit/snakemake_RNAseq_tetoolkit/DESeq2/FDR0.05_L2FC0.0/hta7_v_wt/genes/res_hta7_v_wt_FDR0.05_L2FC0.0_chr_upRegSorted_genes_featureIDs.txt", header = F)[,1]

# Load table of gene coordinates (GFF3)
genes <- readGFF("/data/public_data/arabidopsis/Araport11/Araport11_GFF3_genes_transposons.201606.gff")
genes <- genes[genes$seqid %in% paste0("Chr", 1:5),]
genes <- genes[genes$type == "mRNA",]
print(dim(genes))
#[1] 52060    25

# Obtain frequency of occurrence of each gene parent ID
n_occur <- data.frame(table(unlist(genes$Parent)))

# Obtain gene records for which the gene parent ID occurs only once
genes_unique <-  as.data.frame(
  genes[ unlist(genes$Parent)
    %in% n_occur$Var1[n_occur$Freq == 1],
  ]
)
# Obtain gene records for which the gene parent ID occurs more than once
genes_multi <- as.data.frame(
  genes[ unlist(genes$Parent)
    %in% n_occur$Var1[n_occur$Freq > 1],
  ]
)

# For each gene parent ID in genes_multi, obtain the gene record with the
# longest transcript
# If multiple genes records have the longest transcript,
# keep the first reported one only
genes_multi_list <- mclapply(seq_along(genes_multi[,1]), function(h) {
  genes_multi_ID_all <- genes_multi[ unlist(genes_multi$Parent)
                          == unlist(genes_multi[h,]$Parent),
                        ]
  genes_multi_ID_all[ genes_multi_ID_all$end-genes_multi_ID_all$start
    == max(genes_multi_ID_all$end-genes_multi_ID_all$start),
  ][1,]
}, mc.cores = detectCores())

# Collapse genes_multi_list into single data.frame and remove duplicates
genes_multi_dup <- rbindlist(genes_multi_list)
genes_multi_rep <- unique(as.data.frame(genes_multi_dup))

# Combine into one representative set of gene entries, order,
# and output in GFF3 and BED formats
genes_rep <- rbind(genes_unique, genes_multi_rep)
genes_rep <- genes_rep[ order(genes_rep$seqid,
                              genes_rep$start,
                              genes_rep$end), ]
genes <- genes_rep
genes$seqid <- gsub("^Chr", "", genes$seqid)
nrow(genes)
#[1] 31346

# Convert into GRanges
genesGR <- GRanges(seqnames = genes$seqid,
                   ranges = IRanges(start = genes$start,
                                    end = genes$end),
                   strand = genes$strand,
                   geneID = genes$Parent)

# For each group of peaks, find overlapping genes
# overlapping peaks
HTA7_genes_overlap <- findOverlaps(query = HTA7,
                                   subject = genesGR,
                                   ignore.strand = TRUE,
                                   select = "all")
HTA7_genesGR <- genesGR[unique(subjectHits(HTA7_genes_overlap))]

# Get genes up-regulated in hta7
hta7_upreg_genesGR <- genesGR[genesGR$geneID %in% hta7_upreg]


# Write to BED format
genes_bed <- data.frame(chr = as.character(seqnames(genesGR)),
                        start = as.integer(start(genesGR)-1),
                        end = as.integer(end(genesGR)),
                        name = as.character(genesGR$geneID),
                        score = "NA",
                        strand = as.character(strand(genesGR))) 
HTA7_genes_bed <- data.frame(chr = as.character(seqnames(HTA7_genesGR)),
                             start = as.integer(start(HTA7_genesGR)-1),
                             end = as.integer(end(HTA7_genesGR)),
                             name = as.character(HTA7_genesGR$geneID),
                             score = "NA",
                             strand = as.character(strand(HTA7_genesGR))) 
hta7_upreg_genes_bed <- data.frame(chr = as.character(seqnames(hta7_upreg_genesGR)),
                                   start = as.integer(start(hta7_upreg_genesGR)-1),
                                   end = as.integer(end(hta7_upreg_genesGR)),
                                   name = as.character(hta7_upreg_genesGR$geneID),
                                   score = "NA",
                                   strand = as.character(strand(hta7_upreg_genesGR))) 
write.table(genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(HTA7_genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_overlapping_HTA7_peaks_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(hta7_upreg_genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_upregulated_in_hta7_", region, ".bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)


# Genomic definitions
fai <- read.table("/projects/meiosis/ajt200/TAIR10/TAIR10_chr_all.fa.fai", header = F)
chrs <- fai$V1[1:5]
chrLens <- fai$V2[1:5]

regionGR = GRanges(seqnames = chrs,
                   ranges = IRanges(start = 1, end = chrLens),
                   strand = "*")

# Define function to select randomly positioned loci of the same
# width distribution as features of interest
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLoc_genesGR contains the same number of loci per chromosome as genesGR
ranLoc_genesGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_regionGR <- regionGR[seqnames(regionGR) == chrs[j]]

  chr_genesGR <- genesGR[seqnames(genesGR) == chrs[j]]

  if(length(chr_genesGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_ranLoc_genes_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_regionGR), function(x) {
                                                                      ( start(chr_regionGR[x]) + max(width(chr_genesGR)) + 2000 ) :
                                                                      ( end(chr_regionGR[x]) - max(width(chr_genesGR)) - 2000 )
                                                                    })),
                                             n = length(chr_genesGR))
    chr_ranLoc_genesGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_ranLoc_genes_Start,
                                                width = width(chr_genesGR)),
                               strand = strand(chr_genesGR))
    ranLoc_genesGR <- append(ranLoc_genesGR, chr_ranLoc_genesGR)
  }
}
stopifnot(identical(width(ranLoc_genesGR), width(genesGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLoc_HTA7_genesGR contains the same number of loci per chromosome as HTA7_genesGR
ranLoc_HTA7_genesGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_regionGR <- regionGR[seqnames(regionGR) == chrs[j]]

  chr_HTA7_genesGR <- HTA7_genesGR[seqnames(HTA7_genesGR) == chrs[j]]

  if(length(chr_HTA7_genesGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_ranLoc_HTA7_genes_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_regionGR), function(x) {
                                                                      ( start(chr_regionGR[x]) + max(width(chr_HTA7_genesGR)) + 2000 ) :
                                                                      ( end(chr_regionGR[x]) - max(width(chr_HTA7_genesGR)) - 2000 )
                                                                    })),
                                             n = length(chr_HTA7_genesGR))
    chr_ranLoc_HTA7_genesGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_ranLoc_HTA7_genes_Start,
                                                width = width(chr_HTA7_genesGR)),
                               strand = strand(chr_HTA7_genesGR))
    ranLoc_HTA7_genesGR <- append(ranLoc_HTA7_genesGR, chr_ranLoc_HTA7_genesGR)
  }
}
stopifnot(identical(width(ranLoc_HTA7_genesGR), width(HTA7_genesGR)))

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLoc_hta7_upreg_genesGR contains the same number of loci per chromosome as hta7_upreg_genesGR
ranLoc_hta7_upreg_genesGR <- GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_regionGR <- regionGR[seqnames(regionGR) == chrs[j]]

  chr_hta7_upreg_genesGR <- hta7_upreg_genesGR[seqnames(hta7_upreg_genesGR) == chrs[j]]

  if(length(chr_hta7_upreg_genesGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_ranLoc_hta7_upreg_genes_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_regionGR), function(x) {
                                                                      ( start(chr_regionGR[x]) + max(width(chr_hta7_upreg_genesGR)) + 2000 ) :
                                                                      ( end(chr_regionGR[x]) - max(width(chr_hta7_upreg_genesGR)) - 2000 )
                                                                    })),
                                             n = length(chr_hta7_upreg_genesGR))
    chr_ranLoc_hta7_upreg_genesGR <- GRanges(seqnames = chrs[j],
                               ranges = IRanges(start = chr_ranLoc_hta7_upreg_genes_Start,
                                                width = width(chr_hta7_upreg_genesGR)),
                               strand = strand(chr_hta7_upreg_genesGR))
    ranLoc_hta7_upreg_genesGR <- append(ranLoc_hta7_upreg_genesGR, chr_ranLoc_hta7_upreg_genesGR)
  }
}
stopifnot(identical(width(ranLoc_hta7_upreg_genesGR), width(hta7_upreg_genesGR)))



# Write to BED format
ranLoc_genes_bed <- data.frame(chr = as.character(seqnames(ranLoc_genesGR)),
                               start = as.integer(start(ranLoc_genesGR)-1),
                               end = as.integer(end(ranLoc_genesGR)),
                               name = as.integer(1:length(ranLoc_genesGR)),
                               score = "NA",
                               strand = as.character(strand(ranLoc_genesGR))) 
ranLoc_HTA7_genes_bed <- data.frame(chr = as.character(seqnames(ranLoc_HTA7_genesGR)),
                                    start = as.integer(start(ranLoc_HTA7_genesGR)-1),
                                    end = as.integer(end(ranLoc_HTA7_genesGR)),
                                    name = as.integer(1:length(ranLoc_HTA7_genesGR)),
                                    score = "NA",
                                    strand = as.character(strand(ranLoc_HTA7_genesGR))) 
ranLoc_hta7_upreg_genes_bed <- data.frame(chr = as.character(seqnames(ranLoc_hta7_upreg_genesGR)),
                                          start = as.integer(start(ranLoc_hta7_upreg_genesGR)-1),
                                          end = as.integer(end(ranLoc_hta7_upreg_genesGR)),
                                          name = as.integer(1:length(ranLoc_hta7_upreg_genesGR)),
                                          score = "NA",
                                          strand = as.character(strand(ranLoc_hta7_upreg_genesGR))) 
write.table(ranLoc_genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_", region, "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(ranLoc_HTA7_genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_overlapping_HTA7_peaks_", region, "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(ranLoc_hta7_upreg_genes_bed,
            file = paste0("TAIR10_chr_all_representative_genes_upregulated_in_hta7_", region, "_randomLoci.bed"),
            quote = F, sep = "\t", row.names = F, col.names = F)
