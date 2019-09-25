#
# Code to group AHR ChIP peaks by "mode"
#
# Sudin Bhattacharya 5/22/18


# System command to identify AHR and write out peaks not overlapping with DHSs
system("bedtools intersect -a AHRpeaks_FDR1_noHeader.bed -b ENCFF001WPX_DHS_broadPeaks.bed -v > AHRpeaks_not-in-DHSs.bed")

# Read in AHR peaks info
peaks_info <- read.csv("AHRpeaks_features.csv", header = TRUE)

# Read in AHR peaks bed file
peaks_bed <- read.table("AHRpeaks_FDR1_noHeader.bed", sep="\t")

# Filter bed file by # of AHR peaks
peaks_noAHRE_bed <- peaks_bed[peaks_info$num_dres == 0, ]
peaks_oneAHRE_bed <- peaks_bed[peaks_info$num_dres == 1, ]
peaks_withAHRE_bed <- peaks_bed[peaks_info$num_dres > 0, ]

# Write different peak groups to beds file
write.table(peaks_noAHRE_bed, file="AHRpeaks_noAHREs.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

write.table(peaks_oneAHRE_bed, file="AHRpeaks_oneAHRE.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

write.table(peaks_withAHRE_bed, file="AHRpeaks_withAHREs.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

