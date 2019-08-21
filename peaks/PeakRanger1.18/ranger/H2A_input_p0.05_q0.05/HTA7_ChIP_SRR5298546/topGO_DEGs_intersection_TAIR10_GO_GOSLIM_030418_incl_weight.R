#!/applications/R/R-3.3.2/bin/Rscript

########################################
# Analyse genes for GO term enrichment #
# relative to all TAIR10 genes         #
########################################

# This script is based on a useful post by Avril Coghlan:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html

# Example usage:
# ./topGO_DEGs_intersection_TAIR10_GO_GOSLIM_030418_incl_weight.R BP genes_overlapping_arm_HTA7_peaks_geneIDs.txt 0.05 arm

#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
ont <- args[1]
target <- args[2]
sigLevel <- args[3]
region <- args[4]

suppressMessages(library(GenomicRanges))
suppressMessages(library(topGO))
suppressMessages(library(Rgraphviz))

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


# Read in GO annotations for TAIR10 genes to define "gene universe"
geneID2GO <- readMappings(file = paste0("/projects/ajt200/TAIR10/TAIR10_ATH_GO_GOSLIM_030418_geneID_GOann_reshaped_noC_noM.txt"))
geneID2GO <- geneID2GO[names(geneID2GO) %in% genesGR$geneID]
geneUniverse <- names(geneID2GO)

# Define list of genes of interest; file should contain a single column of gene identifiers
genesOfInterest <- as.character(read.table(target)$x)

# Specify where genes of interest appear in the gene universe vector
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# Build GOdata object in topGO
capture.output(GOdata <- new("topGOdata", description = "Differentially expressed genes", ontology = ont,
                             allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO),
               file="/dev/null")

# Access list of genes of interest
#sg <- sigGenes(GOdata)
#print(str(sg))
#print(numSigGenes(GOdata))

# Run Fisher's exact tests to determine GO term enrichment
capture.output(resultClassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher"),
               file="/dev/null")
capture.output(resultElim <- runTest(GOdata, algorithm = "elim", statistic = "fisher"),
               file="/dev/null")
capture.output(resultWeight <- runTest(GOdata, algorithm = "weight", statistic = "fisher"),
               file="/dev/null")
capture.output(resultTopGO <- runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
               file="/dev/null")

# Count number of results where weight01 gives a P-value <= sigLevel (args[3])
mySummary <- summary(attributes(resultTopGO)$score <= as.numeric(sigLevel))
numSignif <- as.integer(mySummary[[3]])

# List significant results and write to file
capture.output(enrichRes <- GenTable(GOdata, classicFisher = resultClassic,
                                     elimFisher = resultElim,
                                     weightFisher = resultWeight,
                                     topGOFisher = resultTopGO,
                                     orderBy = "topGOFisher",
                                     ranksOf = "elimFisher", topNodes = numSignif),
               file="/dev/null")

# WARNING: ASSUMES INPUT FILE HAS A 3-LETTER EXTENSION
basename <- basename(target)
len <- nchar(basename)
basename <- substr(basename, 1, len-4)

out_name <- paste(basename, "GO", ont, "enrichment.tsv", sep="_")
folder <- paste0(dirname(target), "/GO")
system(paste0("[ -d ", folder, " ] || mkdir ", folder))
folder2 <- paste0(folder, "/", basename, "_GO_", ont)
system(paste0("[ -d ", folder2, " ] || mkdir ", folder2))

capture.output(write.table(enrichRes, file = file.path(folder, out_name), sep = "\t",
                           row.names = FALSE, col.names = TRUE, quote = FALSE),
               file="/dev/null")

# Visualise the positions of the top 5 statistically significant GO terms in the GO hierarchy
out_name2 <- paste(basename, "GO", ont, "enrichment", sep="_")
printGraph(GOdata, resultTopGO, firstSigNodes = 5,
           fn.prefix = file.path(folder, out_name2), useInfo = "all", pdfSW = TRUE)

# Extract gene IDs annotated with significantly enriched GO terms
myTerms <- enrichRes$GO.ID
myGenes <- genesInTerm(GOdata, myTerms)
for(i in 1:length(myTerms)) {
  myTerm <- myTerms[i]
  myGenesForTerm <- myGenes[myTerm][[1]]
  myFactor <- myGenesForTerm %in% genesOfInterest
  myGenesForTermT <- myGenesForTerm[myFactor == TRUE]
  myGenesForTermT <- paste(myGenesForTermT, collapse = ",")
  myGenesForTermT <- paste(myTerm, myGenesForTermT, sep = "\t")
  out_name3 <- paste0(basename, "_GO_", ont, "_enrichment_", myTerm, ".txt")  
  write(myGenesForTermT, file = file.path(folder2, out_name3))
}
