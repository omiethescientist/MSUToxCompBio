# 
# Code to read in bound and unbound AHREs in MCF-7 cells, 
# and get the ChIP reads associated with these AHREs
# 
#   Sudin Bhattacharya 4/26/18

# Read in bound AHREs from CSV file
AHREs_bound <- read.csv(file = "AHRpeak_AHREs_19bp.csv", 
                        header=TRUE, sep=",")

# Read in unbound AHREs from CSV file
AHREs_unbound <- read.csv(file = "DHS-minus-AHR_peaks_AHREs_19bp.csv", 
                          header=TRUE, sep=",")

# Read in bound AHRES BED file
AHREs_bound_bed <- read.csv(file = "AHRpeak_AHREs_19bp.bed", 
                        header=FALSE, sep="\t", skip = 1)

# Drop scores column
cols_to_drop <- c("V4")
AHREs_bound_bed2 <- AHREs_bound_bed[, !(names(AHREs_bound_bed) %in% cols_to_drop)]

# Add unique ID column
num_ids <- dim(AHREs_bound_bed)[1]
ids <- vector()

for (i in 1: num_ids) {
  ids <- c(ids, paste0("AHRE", i))
}

AHREs_bound_bed2$ID <- ids

# Write to new bed file
write.table(AHREs_bound_bed2, file="AHR_peaks_AHREs_19bp_2.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)


