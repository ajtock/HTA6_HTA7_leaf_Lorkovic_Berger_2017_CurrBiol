#!/applications/R/R-4.0.0/bin/Rscript

# Plot heatmaps of features sorted by coverage levels between start and end sites

# Usage:
# /applications/R/R-4.0.0/bin/Rscript features_avgProfileRibbon_noCorr_incl_ranLoc.R 2000 2kb '2 kb' 20 20bp 'SPO11_1_oligos_RPI1' 'SPO11-1' 'dodgerblue2' '0.02,0.94'

#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#libNames <- unlist(strsplit("SPO11_1_oligos_RPI1",
#                            split = ","))
#libNamesPlot <- unlist(strsplit("SPO11-1",
#                                split = ","))
#colours <- unlist(strsplit("dodgerblue2",
#                           split = ","))
### Custom annotation legends using textGrob()
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.94",
#                                        split = ",")))
## middle left upper
#legendPos <- as.numeric(unlist(strsplit("0.02,0.70",
#                                        split = ",")))
## middle left 
#legendPos <- as.numeric(unlist(strsplit("0.02,0.52",
#                                        split = ",")))
## middle left lower
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
## middle right
#legendPos <- as.numeric(unlist(strsplit("0.70,0.54",
#                                        split = ",")))

### ggplot2 legends
## top left
#legendPos <- as.numeric(unlist(strsplit("0.165,0.92",
#                                        split = ",")))
## OR bottom right
#legendPos <- as.numeric(unlist(strsplit("0.84,0.08",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
upstream <- as.numeric(args[1])
downstream <- as.numeric(args[1])
flankName <- args[2]
flankNamePlot <- args[3]
binSize <- as.numeric(args[4])
binName <- args[5]
libNames <- unlist(strsplit(args[6],
                            split = ","))
libNamesPlot <- unlist(strsplit(args[7],
                                split = ","))
colours <- unlist(strsplit(args[8],
                           split = ","))
legendPos <- as.numeric(unlist(strsplit(args[9],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

superfamNames <- system(paste0("ls ", matDir, "SPO11_1_oligos_RPI1*TEs_overlapping_HTA7_HTA6_genomewide_peaks_smoothed*dataframe.txt"),
                        intern = T)
superfamNames <- sub(pattern = "./matrices/SPO11_1_oligos_RPI1_norm_cov_(\\w+)_TEs_overlapping_HTA7_HTA6_genomewide_peaks_smoothed_\\w+_dataframe.txt",
                     replacement = "\\1",
                     x = superfamNames)

for(y in seq_along(superfamNames)) {
  # Load feature1 coverage matrix for each ChIP-seq dataset
  feature1Mats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],                            
                                "_norm_cov_",
                                superfamNames[y],
                                "_TEs_overlapping_HTA7_HTA6_genomewide_peaks_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  
  # Load ranLoc1 coverage matrix for each ChIP-seq dataset
  ranLoc1Mats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],                            
                                "_norm_cov_",
                                superfamNames[y],
                                "_TEs_overlapping_HTA7_HTA6_genomewide_peaks_ranLoc_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  
  ## feature1
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFfeature1_list <- mclapply(seq_along(feature1Mats), function(x) {
    data.frame(window = colnames(feature1Mats[[x]]),
               t(feature1Mats[[x]]))
  }, mc.cores = length(feature1Mats))
  
  # Convert into tidy data.frame (long format)
  tidyDFfeature1_list  <- mclapply(seq_along(wideDFfeature1_list), function(x) {
    gather(data  = wideDFfeature1_list[[x]],
           key   = feature1,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFfeature1_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFfeature1_list)) {
    tidyDFfeature1_list[[x]]$window <- factor(tidyDFfeature1_list[[x]]$window,
                                             levels = as.character(wideDFfeature1_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (feature1s) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFfeature1_list  <- mclapply(seq_along(tidyDFfeature1_list), function(x) {
    data.frame(window = as.character(wideDFfeature1_list[[x]]$window),
               n      = tapply(X     = tidyDFfeature1_list[[x]]$coverage,
                               INDEX = tidyDFfeature1_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature1_list[[x]]$coverage,
                               INDEX = tidyDFfeature1_list[[x]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFfeature1_list[[x]]$coverage,
                               INDEX = tidyDFfeature1_list[[x]]$window,
                               FUN   = sd))
  }, mc.cores = length(tidyDFfeature1_list))
  
  for(x in seq_along(summaryDFfeature1_list)) {
    summaryDFfeature1_list[[x]]$window <- factor(summaryDFfeature1_list[[x]]$window,
                                                levels = as.character(wideDFfeature1_list[[x]]$window))
    summaryDFfeature1_list[[x]]$winNo <- factor(1:dim(summaryDFfeature1_list[[x]])[1])
    summaryDFfeature1_list[[x]]$sem <- summaryDFfeature1_list[[x]]$sd/sqrt(summaryDFfeature1_list[[x]]$n-1)
    summaryDFfeature1_list[[x]]$CI_lower <- summaryDFfeature1_list[[x]]$mean -
      qt(0.975, df = summaryDFfeature1_list[[x]]$n-1)*summaryDFfeature1_list[[x]]$sem
    summaryDFfeature1_list[[x]]$CI_upper <- summaryDFfeature1_list[[x]]$mean +
      qt(0.975, df = summaryDFfeature1_list[[x]]$n-1)*summaryDFfeature1_list[[x]]$sem
  }
  
  names(summaryDFfeature1_list) <- libNamesPlot
  
  # Convert list summaryDFfeature1_list into a single data.frame for plotting
  summaryDFfeature1 <- bind_rows(summaryDFfeature1_list, .id = "libName")
  summaryDFfeature1$libName <- factor(summaryDFfeature1$libName,
                                     levels = names(summaryDFfeature1_list))
  
  
  ## ranLoc1
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFranLoc1_list <- mclapply(seq_along(ranLoc1Mats), function(x) {
    data.frame(window = colnames(ranLoc1Mats[[x]]),
               t(ranLoc1Mats[[x]]))
  }, mc.cores = length(ranLoc1Mats))
  
  # Convert into tidy data.frame (long format)
  tidyDFranLoc1_list  <- mclapply(seq_along(wideDFranLoc1_list), function(x) {
    gather(data  = wideDFranLoc1_list[[x]],
           key   = ranLoc1,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFranLoc1_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFranLoc1_list)) {
    tidyDFranLoc1_list[[x]]$window <- factor(tidyDFranLoc1_list[[x]]$window,
                                            levels = as.character(wideDFranLoc1_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (ranLoc1s) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFranLoc1_list  <- mclapply(seq_along(tidyDFranLoc1_list), function(x) {
    data.frame(window = as.character(wideDFranLoc1_list[[x]]$window),
               n      = tapply(X     = tidyDFranLoc1_list[[x]]$coverage,
                               INDEX = tidyDFranLoc1_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFranLoc1_list[[x]]$coverage,
                               INDEX = tidyDFranLoc1_list[[x]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFranLoc1_list[[x]]$coverage,
                               INDEX = tidyDFranLoc1_list[[x]]$window,
                               FUN   = sd))
  }, mc.cores = length(tidyDFranLoc1_list))
  
  for(x in seq_along(summaryDFranLoc1_list)) {
    summaryDFranLoc1_list[[x]]$window <- factor(summaryDFranLoc1_list[[x]]$window,
                                               levels = as.character(wideDFranLoc1_list[[x]]$window))
    summaryDFranLoc1_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc1_list[[x]])[1])
    summaryDFranLoc1_list[[x]]$sem <- summaryDFranLoc1_list[[x]]$sd/sqrt(summaryDFranLoc1_list[[x]]$n-1)
    summaryDFranLoc1_list[[x]]$CI_lower <- summaryDFranLoc1_list[[x]]$mean -
      qt(0.975, df = summaryDFranLoc1_list[[x]]$n-1)*summaryDFranLoc1_list[[x]]$sem
    summaryDFranLoc1_list[[x]]$CI_upper <- summaryDFranLoc1_list[[x]]$mean +
      qt(0.975, df = summaryDFranLoc1_list[[x]]$n-1)*summaryDFranLoc1_list[[x]]$sem
  }
  
  names(summaryDFranLoc1_list) <- libNamesPlot
  
  # Convert list summaryDFranLoc1_list into a single data.frame for plotting
  summaryDFranLoc1 <- bind_rows(summaryDFranLoc1_list, .id = "libName")
  summaryDFranLoc1$libName <- factor(summaryDFranLoc1$libName,
                                    levels = names(summaryDFranLoc1_list))


  # Load feature2 coverage matrix for each ChIP-seq dataset
  feature2Mats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],                            
                                "_norm_cov_",
                                superfamNames[y],
                                "_TEs_overlapping_HTA7_not_HTA6_genomewide_peaks_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  
  # Load ranLoc2 coverage matrix for each ChIP-seq dataset
  ranLoc2Mats <- mclapply(seq_along(libNames), function(x) {
    as.matrix(read.table(paste0(matDir,
                                libNames[x],                            
                                "_norm_cov_",
                                superfamNames[y],
                                "_TEs_overlapping_HTA7_not_HTA6_genomewide_peaks_ranLoc_smoothed_target_and_",
                                flankName, "_flank_dataframe.txt"),
                         header = T))
  }, mc.cores = length(libNames))
  
  ## feature2
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFfeature2_list <- mclapply(seq_along(feature2Mats), function(x) {
    data.frame(window = colnames(feature2Mats[[x]]),
               t(feature2Mats[[x]]))
  }, mc.cores = length(feature2Mats))
  
  # Convert into tidy data.frame (long format)
  tidyDFfeature2_list  <- mclapply(seq_along(wideDFfeature2_list), function(x) {
    gather(data  = wideDFfeature2_list[[x]],
           key   = feature2,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFfeature2_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFfeature2_list)) {
    tidyDFfeature2_list[[x]]$window <- factor(tidyDFfeature2_list[[x]]$window,
                                             levels = as.character(wideDFfeature2_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (feature2s) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFfeature2_list  <- mclapply(seq_along(tidyDFfeature2_list), function(x) {
    data.frame(window = as.character(wideDFfeature2_list[[x]]$window),
               n      = tapply(X     = tidyDFfeature2_list[[x]]$coverage,
                               INDEX = tidyDFfeature2_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature2_list[[x]]$coverage,
                               INDEX = tidyDFfeature2_list[[x]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFfeature2_list[[x]]$coverage,
                               INDEX = tidyDFfeature2_list[[x]]$window,
                               FUN   = sd))
  }, mc.cores = length(tidyDFfeature2_list))
  
  for(x in seq_along(summaryDFfeature2_list)) {
    summaryDFfeature2_list[[x]]$window <- factor(summaryDFfeature2_list[[x]]$window,
                                                levels = as.character(wideDFfeature2_list[[x]]$window))
    summaryDFfeature2_list[[x]]$winNo <- factor(1:dim(summaryDFfeature2_list[[x]])[1])
    summaryDFfeature2_list[[x]]$sem <- summaryDFfeature2_list[[x]]$sd/sqrt(summaryDFfeature2_list[[x]]$n-1)
    summaryDFfeature2_list[[x]]$CI_lower <- summaryDFfeature2_list[[x]]$mean -
      qt(0.975, df = summaryDFfeature2_list[[x]]$n-1)*summaryDFfeature2_list[[x]]$sem
    summaryDFfeature2_list[[x]]$CI_upper <- summaryDFfeature2_list[[x]]$mean +
      qt(0.975, df = summaryDFfeature2_list[[x]]$n-1)*summaryDFfeature2_list[[x]]$sem
  }
  
  names(summaryDFfeature2_list) <- libNamesPlot
  
  # Convert list summaryDFfeature2_list into a single data.frame for plotting
  summaryDFfeature2 <- bind_rows(summaryDFfeature2_list, .id = "libName")
  summaryDFfeature2$libName <- factor(summaryDFfeature2$libName,
                                     levels = names(summaryDFfeature2_list))
  
  
  ## ranLoc2
  # Transpose matrix and convert to dataframe
  # in which first column is window name
  wideDFranLoc2_list <- mclapply(seq_along(ranLoc2Mats), function(x) {
    data.frame(window = colnames(ranLoc2Mats[[x]]),
               t(ranLoc2Mats[[x]]))
  }, mc.cores = length(ranLoc2Mats))
  
  # Convert into tidy data.frame (long format)
  tidyDFranLoc2_list  <- mclapply(seq_along(wideDFranLoc2_list), function(x) {
    gather(data  = wideDFranLoc2_list[[x]],
           key   = ranLoc2,
           value = coverage,
           -window)
  }, mc.cores = length(wideDFranLoc2_list))
  
  # Order levels of factor "window" so that sequential levels
  # correspond to sequential windows
  for(x in seq_along(tidyDFranLoc2_list)) {
    tidyDFranLoc2_list[[x]]$window <- factor(tidyDFranLoc2_list[[x]]$window,
                                            levels = as.character(wideDFranLoc2_list[[x]]$window))
  }
  
  # Create summary data.frame in which each row corresponds to a window (Column 1),
  # Column2 is the number of coverage values (ranLoc2s) per window,
  # Column3 is the mean of coverage values per window,
  # Column4 is the standard deviation of coverage values per window,
  # Column5 is the standard error of the mean of coverage values per window,
  # Column6 is the lower bound of the 95% confidence interval, and
  # Column7 is the upper bound of the 95% confidence interval
  summaryDFranLoc2_list  <- mclapply(seq_along(tidyDFranLoc2_list), function(x) {
    data.frame(window = as.character(wideDFranLoc2_list[[x]]$window),
               n      = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = mean),
               sd     = tapply(X     = tidyDFranLoc2_list[[x]]$coverage,
                               INDEX = tidyDFranLoc2_list[[x]]$window,
                               FUN   = sd))
  }, mc.cores = length(tidyDFranLoc2_list))
  
  for(x in seq_along(summaryDFranLoc2_list)) {
    summaryDFranLoc2_list[[x]]$window <- factor(summaryDFranLoc2_list[[x]]$window,
                                               levels = as.character(wideDFranLoc2_list[[x]]$window))
    summaryDFranLoc2_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc2_list[[x]])[1])
    summaryDFranLoc2_list[[x]]$sem <- summaryDFranLoc2_list[[x]]$sd/sqrt(summaryDFranLoc2_list[[x]]$n-1)
    summaryDFranLoc2_list[[x]]$CI_lower <- summaryDFranLoc2_list[[x]]$mean -
      qt(0.975, df = summaryDFranLoc2_list[[x]]$n-1)*summaryDFranLoc2_list[[x]]$sem
    summaryDFranLoc2_list[[x]]$CI_upper <- summaryDFranLoc2_list[[x]]$mean +
      qt(0.975, df = summaryDFranLoc2_list[[x]]$n-1)*summaryDFranLoc2_list[[x]]$sem
  }
  
  names(summaryDFranLoc2_list) <- libNamesPlot
  
  # Convert list summaryDFranLoc2_list into a single data.frame for plotting
  summaryDFranLoc2 <- bind_rows(summaryDFranLoc2_list, .id = "libName")
  summaryDFranLoc2$libName <- factor(summaryDFranLoc2$libName,
                                    levels = names(summaryDFranLoc2_list))
  
  
  featureStartLab <- "Start"
  featureEndLab <- "End"
  
  # Define y-axis limits
  ymin <- min(c(summaryDFfeature1$CI_lower,
                summaryDFranLoc1$CI_lower,
                summaryDFfeature2$CI_lower,
                summaryDFranLoc2$CI_lower))
  ymax <- max(c(summaryDFfeature1$CI_upper,
                summaryDFranLoc1$CI_upper,
                summaryDFfeature2$CI_upper,
                summaryDFranLoc2$CI_upper))
  
  # Function for formatting y-axis labels
  # with a given number of decimals
  fmt_decimals <- function(decimals) {
    function(x) format(x, nsmall = decimals, scientific = FALSE)
  }
  
#  # Define legend labels
#  legendLab1 <- grobTree(textGrob(bquote(.(libNamesPlot[1])),
#                                  x = legendPos[1], y = legendPos[2], just = "left",
#                                  gp = gpar(col = colours[1], fontsize = 18))) 
  # Define legend labels
  legendLabs <- lapply(seq_along(libNamesPlot), function(x) {
    grobTree(textGrob(bquote(.(libNamesPlot[x])),
                      x = legendPos[1], y = legendPos[2]-((x-1)*0.08), just = "left",
                      gp = gpar(col = colours[x], fontsize = 22)))
  })

  # Plot average coverage profiles with 95% CI ribbon
  ## feature1
  ggObjGA <- NULL
  ggObj1 <- NULL
  ggObj1 <- ggplot(data = summaryDFfeature1,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_line(data = summaryDFfeature1,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFfeature1,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFfeature1_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFfeature1_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                featureStartLab,
                                featureEndLab,
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFfeature1_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = expression("Log"[2]*"(ChIP/input)")) +
    annotation_custom(legendLabs[[1]]) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("HTA7_HTA6" ~ .(superfamNames[y]) ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFfeature1$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ## ranLoc1
  ggObj2 <- NULL
  ggObj2 <- ggplot(data = summaryDFranLoc1,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_line(data = summaryDFranLoc1,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFranLoc1,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFranLoc1_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFranLoc1_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                "Start",
                                "End",
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFranLoc1_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "") +
    annotation_custom(legendLabs[[1]]) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("Random loci" ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFranLoc1$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ## feature2
  ggObj3 <- NULL
  ggObj3 <- ggplot(data = summaryDFfeature2,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_line(data = summaryDFfeature2,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFfeature2,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFfeature2_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFfeature2_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                featureStartLab,
                                featureEndLab,
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFfeature2_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "") +
    annotation_custom(legendLabs[[1]]) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("HTA7_not_HTA6" ~ .(superfamNames[y]) ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFfeature2$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ## ranLoc2
  ggObj4 <- NULL
  ggObj4 <- ggplot(data = summaryDFranLoc2,
                   mapping = aes(x = winNo,
                                 y = mean,
                                 group = libName),
                  ) +
    geom_line(data = summaryDFranLoc2,
              mapping = aes(colour = libName),
              size = 1) +
    scale_colour_manual(values = colours) +
    geom_ribbon(data = summaryDFranLoc2,
                #mapping = aes(ymin = mean-sem,
                #              ymax = mean+sem,
                mapping = aes(ymin = CI_lower,
                              ymax = CI_upper,
                              fill = libName),
                alpha = 0.4) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(ymin, ymax),
                       labels = function(x) sprintf("%5.2f", x)) +
    scale_x_discrete(breaks = c(1,
                                (upstream/binSize)+1,
                                (dim(summaryDFranLoc2_list[[1]])[1])-(downstream/binSize),
                                dim(summaryDFranLoc2_list[[1]])[1]),
                     labels = c(paste0("-", flankNamePlot),
                                "Start",
                                "End",
                                paste0("+", flankNamePlot))) +
    geom_vline(xintercept = c((upstream/binSize)+1,
                              (dim(summaryDFranLoc2_list[[1]])[1])-(downstream/binSize)),
               linetype = "dashed",
               size = 1) +
    labs(x = "",
         y = "") +
    annotation_custom(legendLabs[[1]]) +
    theme_bw() +
    theme(
          axis.ticks = element_line(size = 1.0, colour = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size = 22, colour = "black"),
          axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
          axis.title = element_text(size = 30, colour = "black"),
          legend.position = "none",
          #legend.text = element_text(size = 10),
          #legend.background = element_rect(fill = "transparent"),
          #legend.key = element_rect(colour = "transparent",
          #                          fill = "transparent"),
          #legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 3.5, colour = "black"),
          panel.background = element_blank(),
          plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
          plot.title = element_text(hjust = 0.5, size = 30)) +
    ggtitle(bquote("Random loci" ~
                   "(" * italic("n") ~ "=" ~
                   .(prettyNum(summaryDFranLoc2$n[1],
                               big.mark = ",", trim = T)) *
                   ")"))
  ggObjGA <- grid.arrange(ggObj1, ggObj2, ggObj3, ggObj4, nrow = 1, ncol = 4)
  ggsave(paste0(plotDir,
                paste(libNames, collapse = "_"),
                "_around_genomewide_HTA7_HTA6_and_HTA7_not_HTA6_peak_overlapping_", superfamNames[y], "_TEs.pdf"),
         plot = ggObjGA,
         height = 6.5, width = 7*4)
}
