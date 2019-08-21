#!/applications/R/R-3.5.0/bin/Rscript

# Adjust "signalValue" (region treatment reads; column 7) by region control reads
# signalValue = (region treatment reads+1)/(region control reads+1)
# Use this signalValue with caution - unknown whether these read counts are
# normalised by library size by PeakRanger
# Preferrable to use qValue (FDR) or pValue in downstream analyses,
# including IDR calculation

# -log10 transform p-values and q-values

# Usage:
# ./3_narrowPeak_sigVal_TreadsNormCreads_minuslog10PQ.R HTA7_ChIP_SRR5298546

args <- commandArgs(trailingOnly = TRUE)
ChIP_prefix <- args[1]

peaks <- read.table(paste0("./", ChIP_prefix, "_rangerPeaks_",
                           "Treads_Creads.narrowPeak.UntransformedPQ"))
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "signalVal", "pValUntrans", "qValUntrans",
                     "summit0based", "Treads", "Creads")
peaks <- cbind(peaks[,1:6], (peaks[,11]+1)/(peaks[,12]+1),
               -log10(peaks[,8]), -log10(peaks[,9]), peaks[,10])
colnames(peaks) <- c("chr", "start0based", "end",
                     "name", "score", "strand",
                     "TreadsNormCreads", "pVal", "qVal",
                     "summit0based")
peaks$pVal[which(!is.finite(peaks$pVal))] <- 323
peaks$qVal[which(!is.finite(peaks$qVal))] <- 323
head(peaks)
write.table(peaks,
            file = paste0("./", ChIP_prefix, "_rangerPeaks_",
                          "TreadsNormCreads.narrowPeak"),
            col.names = F, row.names = F, sep = "\t", quote = F)
