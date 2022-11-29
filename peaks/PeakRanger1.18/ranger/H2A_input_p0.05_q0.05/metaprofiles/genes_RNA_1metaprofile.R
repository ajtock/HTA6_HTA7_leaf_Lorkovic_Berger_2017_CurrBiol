#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 15.09.2022

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)

# Usage:
# conda activate R-4.0.0
# ./genes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_SPO11oligos_Rep1' 'WT_gDNA_Rep1_R1' '160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos' '150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos' 'SPO11-1-oligos' 'gDNA' 'dodgerblue2' 'dodgerblue2' 'TAIR10_chr_all' genomewide
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#binName <- "10bp"
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("WT_SPO11oligos_Rep1",
#                             split = ","))
#controlNames <- unlist(strsplit("WT_gDNA_Rep1_R1",
#                                split = ","))
#ChIPNamesDir <- unlist(strsplit("160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos",
#                                split = ","))
#controlNamesDir <- unlist(strsplit("150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos",
#                                   split = ","))
#ChIPNamesPlot <- unlist(strsplit("SPO11-1-oligos",
#                                 split = ","))
#controlNamesPlot <- unlist(strsplit("gDNA",
#                                    split = ","))
#ChIPColours <- unlist(strsplit("dodgerblue1",
#                               split = ","))
#controlColours <- unlist(strsplit("dodgerblue1",
#                                  split = ","))
#refbase <- unlist(strsplit("TAIR10_chr_all",
#                           split = ","))
#regionName <- "genomewide"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
upstream <- as.numeric(args[4])
downstream <- as.numeric(args[4])
flankName <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
legendPos <- as.numeric(unlist(strsplit(args[8],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[9],
                             split = ","))
controlNames <- unlist(strsplit(args[10],
                                split = ","))
ChIPNamesDir <- unlist(strsplit(args[11],
                                split = ","))
controlNamesDir <- unlist(strsplit(args[12],
                                   split = ","))
ChIPNamesPlot <- unlist(strsplit(args[13],
                                 split = ","))
controlNamesPlot <- unlist(strsplit(args[14],
                                    split = ","))
ChIPColours <- unlist(strsplit(args[15],
                               split = ","))
controlColours <- unlist(strsplit(args[16],
                                  split = ","))
refbase <- unlist(strsplit(args[17],
                           split = ","))
regionName <- args[18]

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

#if(regionName == "genomewide") {
#  regionNamePlot <- "Genome-wide" 
#} else if(regionName == "arm") {
#  regionNamePlot <- "Arm" 
#} else if(regionName == "peri") {
#  regionNamePlot <- "Peri" 
#} else {
#  stop("regionName is not genomewide, arm, or peri")
#}

HTA7_genesNamePlot <- "Overlapping HTA7 peaks"
HTA7_genes_ranLocNamePlot <- "Random loci"
hta7_upreg_genesNamePlot <- "Upregulated in hta7"
hta7_upreg_genes_ranLocNamePlot <- "Random loci"
genesNamePlot <- "All genes"
genes_ranLocNamePlot <- "Random loci"

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/")
})

## ChIP
# HTA7_genes
ChIP_HTA7_genesMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_HTA7_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))
# HTA7_genes_ranLoc
ChIP_HTA7_genes_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_HTA7_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# hta7_upreg_genes
ChIP_hta7_upreg_genesMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_hta7_upreg_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))
# hta7_upreg_genes_ranLoc
ChIP_hta7_upreg_genes_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_hta7_upreg_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))

# genes
ChIP_genesMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))
# genes_ranLoc
ChIP_genes_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "geneProfiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(ChIPNames))


## control
# HTA7_genes
control_HTA7_genesMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_HTA7_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))
# HTA7_genes_ranLoc
control_HTA7_genes_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_HTA7_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))

# hta7_upreg_genes
control_hta7_upreg_genesMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_hta7_upreg_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))
# hta7_upreg_genes_ranLoc
control_hta7_upreg_genes_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_hta7_upreg_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))

# genes
control_genesMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_genes_", regionName,
                                "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))
# genes_ranLoc
control_genes_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
    as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_genes_", regionName,
                                "_ranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
}, mc.cores = length(controlNames))


# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# HTA7_genes
log2ChIP_HTA7_genesMats <- mclapply(seq_along(ChIP_HTA7_genesMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_HTA7_genesMats[[x]]+1)/(control_HTA7_genesMats[[x]]+1))
}, mc.cores = length(ChIP_HTA7_genesMats))
# HTA7_genes_ranLoc
log2ChIP_HTA7_genes_ranLocMats <- mclapply(seq_along(ChIP_HTA7_genes_ranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_HTA7_genes_ranLocMats[[x]]+1)/(control_HTA7_genes_ranLocMats[[x]]+1))
}, mc.cores = length(ChIP_HTA7_genes_ranLocMats))

# hta7_upreg_genes
log2ChIP_hta7_upreg_genesMats <- mclapply(seq_along(ChIP_hta7_upreg_genesMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_hta7_upreg_genesMats[[x]]+1)/(control_hta7_upreg_genesMats[[x]]+1))
}, mc.cores = length(ChIP_hta7_upreg_genesMats))
# hta7_upreg_genes_ranLoc
log2ChIP_hta7_upreg_genes_ranLocMats <- mclapply(seq_along(ChIP_hta7_upreg_genes_ranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_hta7_upreg_genes_ranLocMats[[x]]+1)/(control_hta7_upreg_genes_ranLocMats[[x]]+1))
}, mc.cores = length(ChIP_hta7_upreg_genes_ranLocMats))

# genes
log2ChIP_genesMats <- mclapply(seq_along(ChIP_genesMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_genesMats[[x]]+1)/(control_genesMats[[x]]+1))
}, mc.cores = length(ChIP_genesMats))
# genes_ranLoc
log2ChIP_genes_ranLocMats <- mclapply(seq_along(ChIP_genes_ranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_genes_ranLocMats[[x]]+1)/(control_genes_ranLocMats[[x]]+1))
}, mc.cores = length(ChIP_genes_ranLocMats))



# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_HTA7_genesMats)) {
  colnames(log2ChIP_HTA7_genesMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                              paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                              paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_HTA7_genes_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                     paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                     paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_hta7_upreg_genesMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                    paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                    paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_hta7_upreg_genes_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                           paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                           paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))

  colnames(log2ChIP_genesMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                         paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                         paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
  colnames(log2ChIP_genes_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                                paste0("t", ((upstream/binSize)+1):((upstream+bodyLength)/binSize)),
                                                paste0("d", (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIP_HTA7_genesMats), function(x) {
  list(
       # HTA7_genes
       log2ChIP_HTA7_genesMats[[x]],
       # HTA7_genes_ranLocs
       log2ChIP_HTA7_genes_ranLocMats[[x]],
       # hta7_upreg_genes
       log2ChIP_hta7_upreg_genesMats[[x]],
       # hta7_upreg_genes_ranLocs
       log2ChIP_hta7_upreg_genes_ranLocMats[[x]],
       # genes
       log2ChIP_genesMats[[x]],
       # genes_ranLocs
       log2ChIP_genes_ranLocMats[[x]]
      )
}, mc.cores = length(log2ChIP_HTA7_genesMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                              levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting
HTA7_genesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
HTA7_genes_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
hta7_upreg_genesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
hta7_upreg_genes_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
genesTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[5]]
})
genes_ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[6]]
})
names(HTA7_genesTmp) <- log2ChIPNamesPlot
names(HTA7_genes_ranLocTmp) <- log2ChIPNamesPlot
names(hta7_upreg_genesTmp) <- log2ChIPNamesPlot
names(hta7_upreg_genes_ranLocTmp) <- log2ChIPNamesPlot
names(genesTmp) <- log2ChIPNamesPlot
names(genes_ranLocTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(HTA7_genesTmp, .id = "libName"),
  bind_rows(HTA7_genes_ranLocTmp, .id = "libName"),
  bind_rows(hta7_upreg_genesTmp, .id = "libName"),
  bind_rows(hta7_upreg_genes_ranLocTmp, .id = "libName"),
  bind_rows(genesTmp, .id = "libName"),
  bind_rows(genes_ranLocTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
                       summaryDFfeature_log2ChIP[[4]]$CI_lower,
                       summaryDFfeature_log2ChIP[[5]]$CI_lower,
                       summaryDFfeature_log2ChIP[[6]]$CI_lower),
                     na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
                       summaryDFfeature_log2ChIP[[4]]$CI_upper,
                       summaryDFfeature_log2ChIP[[5]]$CI_upper,
                       summaryDFfeature_log2ChIP[[6]]$CI_upper),
                     na.rm = T)


# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 22)))
})


# Plot average profiles with 95% CI ribbon

## HTA7_genes
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(HTA7_genesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## HTA7_genes_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(HTA7_genes_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## hta7_upreg_genes
summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(hta7_upreg_genesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## hta7_upreg_genes_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(hta7_upreg_genes_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## genes
summaryDFfeature <- summaryDFfeature_log2ChIP[[5]]
ggObj5_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            geneStartLab,
                            geneEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(genesNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## genes_ranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[6]]
ggObj6_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%4.2f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[6]])[1]/length(log2ChIPNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(genes_ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP,
                                              ggObj5_combined_log2ChIP,
                                              ggObj6_combined_log2ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5,
                                                       6
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_",
              paste0(log2ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_HTA7_genes__hta7_upreg_genes__all_genes_", regionName, refbase,
              "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*6, limitsize = FALSE)
