#!/usr/bin/env Rscript

# Convert GRanges objects to gff files

# Usage:
# source activate R-4.0.0
# ./GRangesToBEDandGFF3.R HTA6_ChIP_SRR5298545 arm 'Chr1,Chr2,Chr3,Chr4,Chr5'
# ./GRangesToBEDandGFF3.R HTA7_ChIP_SRR5298546 arm 'Chr1,Chr2,Chr3,Chr4,Chr5'
# conda deactivate

#libName = "HTA6_ChIP_SRR5298545"
#region = "arm"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))

args = commandArgs(trailingOnly=T)
libName = args[1]
region = args[2]
chrName <- unlist(strsplit(args[3],
                           split = ","))

gff_dir = "GFF3/"
bed_dir = "BED/"
system(paste0("[ -d ", gff_dir, " ] || mkdir -p ", gff_dir))
system(paste0("[ -d ", bed_dir, " ] || mkdir -p ", bed_dir))

options(stringsAsFactors=F)
library(GenomicRanges)

# Genomic definitions
fai = read.table("/projects/meiosis/ajt200/TAIR10/bowtie2_index/TAIR10_chr_all.fa.fai", header=F)[1:5,]
fai$V1 = paste0("Chr", fai$V1) 
chrs = fai$V1
chrLens = fai$V2
cen_pericen  = read.csv("/projects/meiosis/ajt200/TAIR10/Ziolkowski_2017_cen_pericen_coords.csv", header=T)[,2:6]
colnames(cen_pericen)[1:5] = c("chr", "pericen_start", "cen_start", "cen_end", "pericen_end")
cen_pericen$chr = paste0("Chr", cen_pericen$chr)
pericenStart = cen_pericen$pericen_start
pericenEnd = cen_pericen$pericen_end


# Define region to be analysed
if(region == "arm") {
  regionGR = GRanges(seqnames=rep(chrs, 2),
                     ranges=IRanges(start=c(rep(1, length(chrs)),
                                            pericenEnd+1),
                                    end=c(pericenStart-1,
                                          chrLens)),
                     strand="*")
  maskGR = GRanges(seqnames=chrs,
                   ranges=IRanges(start=pericenStart,
                                  end=pericenEnd),
                   strand="*")
} else if(region == "peri") {
  regionGR = GRanges(seqnames=chrs,
                     ranges=IRanges(start=pericenStart,
                                    end=pericenEnd),
                     strand="*")
  maskGR = GRanges(seqnames=rep(chrs, 2),
                   ranges=IRanges(start=c(rep(1, length(chrs)),
                                          pericenEnd+1),
                                  end=c(pericenStart-1,
                                        chrLens)),
                   strand="*")
} else if(region == "genomewide") {
  regionGR = GRanges(seqnames=chrs,
                     ranges=IRanges(start=1,
                                    end=chrLens),
                     strand="*")
  maskGR = GRanges()
} else {
  stop("region is not arm, peri, or genomewide")
}


# Import peaks GRanges object
if(region == "arm") {
  load(paste0(libName, "_rangerPeaksGR_", region, "_mergedOverlaps_noMinWidth.RData"))
  peaksGR = rangerPeaksGR_arm_mergedOverlaps
  rm(rangerPeaksGR_arm_mergedOverlaps)
} else if(region == "peri") {
  load(paste0(libName, "_rangerPeaksGR_", region, "_mergedOverlaps_noMinWidth.RData"))
  peaksGR = rangerPeaksGR_peri_mergedOverlaps
  rm(rangerPeaksGR_peri_mergedOverlaps)
} else {
  load(paste0(libName, "_rangerPeaksGR_arm_mergedOverlaps_noMinWidth.RData"))
  load(paste0(libName, "_rangerPeaksGR_peri_mergedOverlaps_noMinWidth.RData"))
  peaksGR = sort(c(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps))
  rm(rangerPeaksGR_arm_mergedOverlaps, rangerPeaksGR_peri_mergedOverlaps)
}
strand(peaksGR) = "*"


# Define function to select randomly positioned loci of the same
# width distribution as peaksGR
ranLocStartSelect = function(coordinates, n) {
  sample(x=coordinates,
         size=n,
         replace=FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen=100)

# Apply ranLocStartSelect() on a per-chromosome basis so that
# ranLocGR contains the same number of loci per chromosome as peaksGR
chrs = chrs[which(chrs %in% chrName)]
peaksGR = peaksGR[seqnames(peaksGR) %in% chrName]
ranLocGR = GRanges()
for(j in 1:length(chrs)) {
  print(chrs[j])

  chr_peaksGR = peaksGR[seqnames(peaksGR) == chrs[j]]

  chr_regionGR = regionGR[seqnames(regionGR) == chrs[j]]

  if(length(chr_peaksGR) > 0) {
    # Define seed so that random selections are reproducible
    set.seed(93750174)
    chr_ranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_regionGR), function(x) {
                                                                 ( start(chr_regionGR[x]) + max(width(chr_peaksGR)) + 2000 ) :
                                                                 ( end(chr_regionGR[x]) - max(width(chr_peaksGR)) - 2000 )
                                                               })),
                                          n = length(chr_peaksGR))
    chr_ranLocGR <- GRanges(seqnames = chrs[j],
                            ranges = IRanges(start = chr_ranLoc_Start,
                                             width = width(chr_peaksGR)),
                            strand = strand(chr_peaksGR))
    ranLocGR <- append(ranLocGR, chr_ranLocGR)
  }
}
stopifnot(identical(width(ranLocGR), width(peaksGR)))


# peaks
# Convert into GFF3 and BED formats
peaks = data.frame(peaksGR)
# Remove "Chr" for compatibility with TAIR10 files
peaks$seqnames = gsub("Chr", "", peaks$seqnames)

# GFF3
peaks_gff = data.frame(seqid=as.character(peaks$seqnames),
                       source=as.character("."),
                       type=as.character(paste0(libName, "_peak")),
                       start=as.integer(peaks$start),
                       end=as.integer(peaks$end),
                       score=as.character("."),
                       strand=as.character("."),
                       phase=as.character("."),
                       attribute=as.character("."))
write.table(peaks_gff,
            file=paste0(gff_dir, libName, "_peaks_", region, "_",  
                        paste0(chrName, collapse="_"), ".gff3"),
            quote=F, sep="\t", row.names=F, col.names=F)

# BED
peaks_bed = data.frame(chr=as.character(peaks$seqnames),
                       start=as.integer(peaks$start-1),
                       end=as.integer(peaks$end),
                       name=as.integer(1:nrow(peaks)),
                       score=rep("NA", nrow(peaks)),
                       strand=as.character(peaks$strand))
write.table(peaks_bed,
            file=paste0(bed_dir, libName, "_peaks_", region, "_",  
                        paste0(chrName, collapse="_"), ".bed"),
            quote=F, sep="\t", row.names=F, col.names=F)


# ranLoc
# Convert into GFF3 and BED formats
ranLoc = data.frame(ranLocGR)
# Remove "Chr" for compatibility with TAIR10 files
ranLoc$seqnames = gsub("Chr", "", ranLoc$seqnames)

# GFF3
ranLoc_gff = data.frame(seqid=as.character(ranLoc$seqnames),
                        source=as.character("."),
                        type=as.character(paste0(libName, "_ranLoc")),
                        start=as.integer(ranLoc$start),
                        end=as.integer(ranLoc$end),
                        score=as.character("."),
                        strand=as.character("."),
                        phase=as.character("."),
                        attribute=as.character("."))
write.table(ranLoc_gff,
            file=paste0(gff_dir, libName, "_peaks_", region, "_",  
                        paste0(chrName, collapse="_"), "_randomLoci.gff3"),
            quote=F, sep="\t", row.names=F, col.names=F)

# BED
ranLoc_bed = data.frame(chr=as.character(ranLoc$seqnames),
                        start=as.integer(ranLoc$start-1),
                        end=as.integer(ranLoc$end),
                        name=as.integer(1:nrow(ranLoc)),
                        score=rep("NA", nrow(ranLoc)),
                        strand=as.character(ranLoc$strand))
write.table(ranLoc_bed,
            file=paste0(bed_dir, libName, "_peaks_", region, "_",  
                        paste0(chrName, collapse="_"), "_randomLoci.bed"),
            quote=F, sep="\t", row.names=F, col.names=F)
