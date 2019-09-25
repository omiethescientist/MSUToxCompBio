#
# Code to identify singleton AHREs under DHSs and use them for classification
#
# Sudin Bhattacharya 5/19/18

library(seqLogo)
library(seqinr)

# Read in AHREs that are under AHR ChIP peaks with only one AHRE
singletonAHREs <- read.csv("singletonAHREs_info.csv", header = TRUE)
head(singletonAHREs)

# ---- Write these AHREs to bed file ------
# Create new data frame from longSites with only necessary fields
singletonAHREs2 <- data.frame(chr = singletonAHREs$chrom, 
                              start = singletonAHREs$start_coord,
                              end = singletonAHREs$end_coord,
                              score = singletonAHREs$siteMS_score)
singletonAHREs2$score <- format(singletonAHREs2$score, digits = 3)
head(singletonAHREs2)

# Write to bed file
write.table(singletonAHREs2, file="singletonAHREs_info.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)


#--- System calls to Unix
# Check bed file
system("head singletonAHREs_info.bed")

# Compute intersect of singleton AHREs bed file and DHS bed file
system("bedtools intersect -a singletonAHREs_info.bed -b ENCFF001WPX_DHS_broadPeaks.bed > singletonAHREs_DHS_overlap.bed")

# Check bed file
system("head singletonAHREs_DHS_overlap.bed")
system("wc singletonAHREs_DHS_overlap.bed")

# Read in singleton AHRE-DHS overlap regions
singletonAHREs_DHS_overlap <- read.table("singletonAHREs_DHS_overlap.bed", sep="\t")
# Rename columns
colnames(singletonAHREs_DHS_overlap) <- c("chrom", "start_coord", "end_coord", "siteMS_score")
# View data frame
head(singletonAHREs_DHS_overlap)

# Filter singleton AHREs to get sequences of those in open chromatin: either start or end overlaps
singletonAHREs_open <- singletonAHREs[((singletonAHREs$start_coord %in% singletonAHREs_DHS_overlap$start_coord)
                                       | (singletonAHREs$end_coord %in% singletonAHREs_DHS_overlap$end_coord)), ]
dim(singletonAHREs_open)

# Save open singleton AHREs to CSV file
write.csv(singletonAHREs_open, file = "singletonAHREs_open.csv", row.names = FALSE)

# Read in pos and neg AHREs
positive_AHREs <- singletonAHREs_open
negative_AHREs_all <- read.csv("DHS-minus-AHR_peaks_AHREs_19bp.csv", header = TRUE)

num_pos_AHREs <- nrow(positive_AHREs)
num_neg_AHREs <- 800L

# Shuffle and draw num_neg_AHREs AHREs from neg AHRE list
negative_AHREs_shuff <- negative_AHREs_all[sample(nrow(negative_AHREs_all)), ]
negative_AHREs <- negative_AHREs_shuff[1:num_neg_AHREs, ]

# Save negative AHREs to CSV file
write.csv(negative_AHREs, file = "negative_AHREs.csv", row.names = FALSE)

## ----------------
#- Create sequence logo for *pos* AHREs
# Initialize vectors of A, C, G, T counts in pos AHREs
AHRES_pos_A_count <- rep(0, 19)
AHRES_pos_C_count <- rep(0, 19)
AHRES_pos_G_count <- rep(0, 19)
AHRES_pos_T_count <- rep(0, 19)

# Loop through pos AHREs
for(AHRE in positive_AHREs$seq) {
  AHRE_string <- s2c(AHRE)
  for(j in 1:length(AHRE_string)) {
    if(AHRE_string[j] == 'A') {
      AHRES_pos_A_count[j] = AHRES_pos_A_count[j] + 1
    } else if (AHRE_string[j] == 'C') {
      AHRES_pos_C_count[j] = AHRES_pos_C_count[j] + 1
    } else if (AHRE_string[j] == 'G') {
      AHRES_pos_G_count[j] = AHRES_pos_G_count[j] + 1
    } else if (AHRE_string[j] == 'T') {
      AHRES_pos_T_count[j] = AHRES_pos_T_count[j] + 1
    } 
  }
}

# Create data frame using the four vectors
pos_df <- data.frame(AHRES_pos_A_count, AHRES_pos_C_count, 
                     AHRES_pos_G_count, AHRES_pos_T_count)

# Define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

# Create position weight matrix and sequence logo figure
pos_pwm <- apply(pos_df, 1, proportion)
pos_pwm <- makePWM(pos_pwm)

# Create sequence logo
seqLogo(pos_pwm, ic.scale = FALSE)


## ----------------
#- Create sequence logo for *neg* AHREs
# Initialize vectors of A, C, G, T counts in neg AHREs
AHRES_neg_A_count <- rep(0, 19)
AHRES_neg_C_count <- rep(0, 19)
AHRES_neg_G_count <- rep(0, 19)
AHRES_neg_T_count <- rep(0, 19)

# Loop through neg AHREs
for(AHRE in negative_AHREs$seq) {
  AHRE_string <- s2c(AHRE)
  for(j in 1:length(AHRE_string)) {
    if(AHRE_string[j] == 'A') {
      AHRES_neg_A_count[j] = AHRES_neg_A_count[j] + 1
    } else if (AHRE_string[j] == 'C') {
      AHRES_neg_C_count[j] = AHRES_neg_C_count[j] + 1
    } else if (AHRE_string[j] == 'G') {
      AHRES_neg_G_count[j] = AHRES_neg_G_count[j] + 1
    } else if (AHRE_string[j] == 'T') {
      AHRES_neg_T_count[j] = AHRES_neg_T_count[j] + 1
    } 
  }
}

# Create data frame using the four vectors
neg_df <- data.frame(AHRES_neg_A_count, AHRES_neg_C_count, 
                     AHRES_neg_G_count, AHRES_neg_T_count)

# # Define function that divides the frequency by the row sum i.e. proportions
# proportion <- function(x){
#   rs <- sum(x);
#   return(x / rs);
# }

# Create position weight matrix and sequence logo figure
neg_pwm <- apply(neg_df, 1, proportion)
neg_pwm <- makePWM(neg_pwm)

# Create sequence logo
seqLogo(neg_pwm, ic.scale = FALSE)


## ----------------
# Code pos AHREs numerically (one hot encoding)
for (AHRE in positive_AHREs$seq) {
  AHRE_string <- s2c(AHRE)
  AHRE_coded <- character()  # create an empty char vector
  for (n in AHRE_string) {
    if (n == "A") {
      AHRE_coded <- c(AHRE_coded, c(1,0,0,0)) 
    } else if (n == "T") {
      AHRE_coded <- c(AHRE_coded, c(0,1,0,0))
    } else if (n == "G") {
      AHRE_coded <- c(AHRE_coded, c(0,0,1,0))
    } else if (n == "C") {
      AHRE_coded <- c(AHRE_coded, c(0,0,0,1))
    } else {
      print("Warning! Found aberrant character in AHRE sequence!!")
    }
  }
  # append coded site to file
  write(AHRE_coded, file="pos_AHREs-open_coded.csv", ncolumns = length(AHRE_coded), 
        append = TRUE, sep=",")
}

# code neg AHREs numerically (one hot encoding)
for (AHRE in negative_AHREs$seq) {
  AHRE_string <- s2c(AHRE)
  AHRE_coded <- character()  # create an empty char vector
  for (n in AHRE_string) {
    if (n == "A") {
      AHRE_coded <- c(AHRE_coded, c(1,0,0,0)) 
    } else if (n == "T") {
      AHRE_coded <- c(AHRE_coded, c(0,1,0,0))
    } else if (n == "G") {
      AHRE_coded <- c(AHRE_coded, c(0,0,1,0))
    } else if (n == "C") {
      AHRE_coded <- c(AHRE_coded, c(0,0,0,1))
    } else {
      print("Warning! Found aberrant character in AHRE sequence!!")
    }
  }
  # append coded site to file
  write(AHRE_coded, file="neg_AHREs-open_coded.csv", ncolumns = length(AHRE_coded), 
        append = TRUE, sep=",")
}



