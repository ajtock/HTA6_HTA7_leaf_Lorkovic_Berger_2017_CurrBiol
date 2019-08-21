#!/applications/R/R-3.4.0/bin/Rscript

### Perform hypergeometric tests to determine whether
# genes up-regulated in meiocytes are over-represented or
# under-represented among genes overlapping HTA7 and/or H3K9me2 peaks

# Over-representation:
# P-value is the probability of drawing >= length(geneIDintersection) [x] genes
# in a sample size of length(geneIDs) [k] from a total gene set consisting of 
# length(upRegMeio_geneIDs) [m] + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs)) [n]

# Under-representation:
# P-value is the probability of drawing <= length(geneIDintersection) [q] genes
# in a sample size of length(geneIDs) [k] from a total gene set consisting of
# length(upRegMeio_geneIDs) [m] + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs)) [n]

# Usage: ./proportion_upReg_meiocyte_hypergeometricTest.R genes_overlapping_arm_HTA7_peaks_geneIDs.txt arm upReg_meiocyte FDR0.01_L2FC0.0 100000

library(methods)
library(GenomicRanges)

#geneIDsFile <- "genes_overlapping_arm_HTA7_peaks_geneIDs.txt"
#region <- "arm"
#outDir <- "upReg_meiocyte"
#meioFDR <- "FDR0.01_L2FC0.0"
#samples <- 100000

args <- commandArgs(trailingOnly = TRUE)
geneIDsFile <- args[1]
region <- args[2]
outDir <- args[3]
meioFDR <- args[4]
# Number of randomisations to perform
samples <- as.numeric(args[5])

outDir <- paste0(dirname(geneIDsFile), "/", outDir, "/")
plotDir <- paste0(dirname(geneIDsFile), "/", outDir, "/plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

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
# Find and remove genes overlapping masked regions
mask_genes_overlap <- findOverlaps(query = maskGR,
                                   subject = genesGR,
                                   ignore.strand = TRUE,
                                   select = "all")
genesGR <- genesGR[-subjectHits(mask_genes_overlap)]

TAIR10_geneIDsChr <- genesGR$geneID

# Genes up-regulated in meiocytes (on chromosomes only)
meioDir <- paste0("/projects/ajt200/BAM_masters/RNAseq_meiocyte_Walker_Feng_2018_NatGenet/WT/snakemake_RNAseq_tetoolkit/DESeq2/",
                  meioFDR, "/wt_meiocyte_v_wt_leaf/genes/")
upRegMeio_geneIDs <- as.character(read.table(paste0(meioDir,
                                                    "res_wt_meiocyte_v_wt_leaf_", as.character(meioFDR),
                                                    "_chr_upRegSorted_genes_featureIDs.txt"))$V1)
# Subset to exclude those in masked regions
upRegMeio_geneIDs <- upRegMeio_geneIDs[upRegMeio_geneIDs %in% TAIR10_geneIDsChr]
print(length(upRegMeio_geneIDs))

# Query gene IDs
geneIDs <- as.character(read.table(geneIDsFile)$x)

# Obtain intersection of gene IDs in query gene set and
# gene IDs of genes up-regulated in meiocytes 
geneIDintersection <- intersect(geneIDs, upRegMeio_geneIDs)
print(length(geneIDintersection))
#[1] 262
# WARNING: ASSUMES INPUT FILE HAS 3-LETTER EXTENSION
baseName <- basename(geneIDsFile)
baseName <- substr(baseName, 1, nchar(baseName)-4)
write.table(geneIDintersection,
            file = paste0(outDir,
                          baseName,
                          "_res_wt_meiocyte_v_wt_leaf_", as.character(meioFDR),
                          "_chr_upRegSorted_genes_featureIDs.txt"))

# Calculate proportion of query genes up-regulated in meiocytes
proportion_upRegMeio <- length(geneIDintersection)/length(geneIDs)

## Use the hypergeometric distribution to calculate P-values for
# over-representation or under-representation of genes up-regulated
# in meiocytes among genes overlapping peaks

# Over-representation:
# P-value is the probability of drawing >= length(geneIDintersection) [x] genes
# in a sample size of length(geneIDs) [k] from a total gene set consisting of 
# length(upRegMeio_geneIDs) [m] + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs)) [n]
m <- length(upRegMeio_geneIDs); print(m)
n <- length(TAIR10_geneIDsChr)-m; print(n)
k <- length(geneIDs); print(k)

# From Karl Broman's answer at
# https://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c:
# dhyper(x, m, n, k) gives the probability of drawing exactly x.
# So P-value is given by sum of the probabilities of drawing
# length(geneIDintersection) to length(geneIDs)
Pval_overrep <- sum(dhyper(x = length(geneIDintersection):k,
                           m = m,
                           n = n,
                           k = k))
print(Pval_overrep)
#[1] 0.9851184
# Or by 1 minus the sum of probabilities of drawing 0:(length(geneIDintersection)-1)
print(1 - sum(dhyper(x = 0:(length(geneIDintersection)-1),
                     m = m,
                     n = n,
                     k = k)))
#[1] 0.9851184
# phyper(q, m, n, k) gives the probability of drawing <= q,
# so phyper(q, m, n, k) is the same as sum(dhyper(0:x, m, n, k)),
# where q and x are equal.
# phyper(q, m, n, k, lower.tail = FALSE) is the same as 
# 1 - phyper(q, m, n, k), and so is the probablity of drawing >= q+1:
print(phyper(q = length(geneIDintersection)-1,
             m = m,
             n = n,
             k = k,
             lower.tail = FALSE))
#[1] 0.9851184
print(1 - phyper(q = length(geneIDintersection)-1,
                 m = m,
                 n = n,
                 k = k))
#[1] 0.9851184

# Under-representation:
# P-value is the probability of drawing <= length(geneIDintersection) [q] genes
# in a sample size of length(geneIDs) [k] from a total gene set consisting of
# length(upRegMeio_geneIDs) [m] + (length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs)) [n]
Pval_underrep <- phyper(q = length(geneIDintersection),
                        m = m,
                        n = n,
                        k = k,
                        lower.tail = TRUE)
print(Pval_underrep)

# Sample without replacement
# (set seed for reproducible sampling)
set.seed(834753)
hgDist <- rhyper(nn = samples,
                 m = length(upRegMeio_geneIDs),
                 n = length(TAIR10_geneIDsChr)-length(upRegMeio_geneIDs),
                 k = length(geneIDs))
random_proportions_upRegMeio <- hgDist/length(geneIDs)

setClass("hypergeomTest",
         representation(Pval_overrep = "numeric",
                        Pval_underrep = "numeric",
                        proportion_upRegMeio = "numeric",
                        random_proportions_upRegMeio = "numeric",
                        geneIDintersectionLength = "numeric",
                        hypergeometricDistribution = "numeric"))
hgTestResults <- new("hypergeomTest",
                     Pval_overrep = Pval_overrep,
                     Pval_underrep = Pval_underrep,
                     proportion_upRegMeio = proportion_upRegMeio,
                     random_proportions_upRegMeio = random_proportions_upRegMeio, 
                     geneIDintersectionLength = length(geneIDintersection),
                     hypergeometricDistribution = hgDist)

basenameFile <- basename(geneIDsFile)
len <- nchar(basenameFile)
# WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
basenameFile <- substr(basenameFile, 1, len-4)
save(hgTestResults,
     file = paste0(outDir, basenameFile,
                   "_hypergeometricTestResults_upRegMeioSig", as.character(meioFDR), ".RData"))

library(plotrix)
# Generate histogram
pdf(paste0(plotDir, "hist_", basenameFile,
           "_hypergeometricTestResults_upRegMeioSig", as.character(meioFDR), ".pdf"),
           height = 4, width = 5)
# Calculate max density
maxDensityPlus <- max(density(hgTestResults@random_proportions_upRegMeio)$y)*1.2
# Conditionally define P-value and alpha0.05 tail to be plotted
if(hgTestResults@proportion_upRegMeio > mean(hgTestResults@random_proportions_upRegMeio)) {
  alpha0.05 <- quantile(hgTestResults@random_proportions_upRegMeio, 0.95)[[1]]
  xlim <- c(pmax(0, min(hgTestResults@random_proportions_upRegMeio)-.1),
            pmax(hgTestResults@proportion_upRegMeio+.1, alpha0.05+.1))
  Pval <- as.character(round(Pval_overrep, 8))
} else {
  alpha0.05 <- quantile(hgTestResults@random_proportions_upRegMeio, 0.05)[[1]]
  xlim <- c(pmax(0, pmin(hgTestResults@proportion_upRegMeio-.1, alpha0.05-.1)),
            max(hgTestResults@random_proportions_upRegMeio)+.1)
  Pval <- as.character(round(Pval_underrep, 8))
}

par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(hgTestResults@random_proportions_upRegMeio,
     freq = FALSE,
     col = "grey70",
     border = NA,
     lwd = 2,
     xlim = xlim,
     ylim = c(0,
              maxDensityPlus),
     xlab = paste0("Proportion of genes \nup-regulated in meiocytes (",
                   as.character(meioFDR), ")"),
     ylab = "Density",
     main = "",
     cex.lab = 1, cex.axis = 1)
## Disable scientific notation (e.g., 100000 samples rather than 1e+05 samples)
options(scipen = 100)
titleText <- list(bquote(.(basenameFile)),
                  bquote(italic("P")~"="~.(Pval)),
                  bquote("Samples (hypergeometric distribution) ="~
                         .(prettyNum(samples,
                                     big.mark = ",", trim = T))))
mtext(do.call(expression, titleText), side = 3, line = 2:0, cex = 1)
lines(density(hgTestResults@random_proportions_upRegMeio), lwd = 1.5)
ablineclip(v = mean(hgTestResults@random_proportions_upRegMeio),
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2)
ablineclip(v = hgTestResults@proportion_upRegMeio,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, col = "forestgreen")
ablineclip(v = alpha0.05,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2, lty = 5, col = "red")
text(x = c(pmax(0.05, min(hgTestResults@random_proportions_upRegMeio)-.05),
           mean(hgTestResults@random_proportions_upRegMeio),
           hgTestResults@proportion_upRegMeio,
           alpha0.05),
     y = c(maxDensityPlus*.95,
           maxDensityPlus,
           maxDensityPlus,
           maxDensityPlus*.95),
     labels = c("Simulated",
                "Expected",
                "Observed",
                expression(alpha~"= 0.05")),
     col = c("grey70",
             "black",
             "forestgreen",
             "red"),
     cex = 0.7)
box(lwd = 2)
dev.off()

