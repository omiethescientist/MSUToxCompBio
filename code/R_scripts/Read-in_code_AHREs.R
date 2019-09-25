#
# Code to read in positive (under AHR peak) and negative (not under AHR peak) AHREs
# Also, create seqLogos
#
# Sudin Bhattacharya 5/15/18

library(plotly)
library(seqinr)
library(ggbiplot)
library(seqLogo)

# Read in pos and neg AHREs
positive_AHREs <- read.csv("singletonAHREs_info.csv", header = TRUE)
negative_AHREs_all <- read.csv("DHS-minus-AHR_peaks_AHREs_19bp.csv", header = TRUE)

num_pos_AHREs <- nrow(positive_AHREs)
num_neg_AHREs <- 1000

# Shuffle and draw 900 AHREs from neg AHRE list
negative_AHREs_shuff <- negative_AHREs_all[sample(nrow(negative_AHREs_all)), ]
negative_AHREs <- negative_AHREs_shuff[1:num_neg_AHREs, ]

# code pos AHREs numerically (one hot encoding)
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
  write(AHRE_coded, file="pos_AHREs_coded.csv", ncolumns = length(AHRE_coded), 
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
  write(AHRE_coded, file="neg_AHREs_coded.csv", ncolumns = length(AHRE_coded), 
        append = TRUE, sep=",")
}

# Read in coded AHREs
pos_AHREs_coded <- read.csv("pos_AHREs_coded.csv", header = FALSE)
neg_AHREs_coded <- read.csv("neg_AHREs_coded.csv", header = FALSE)

# Add Bound_status column to coded AHREs (1 if bound; 0 if unbound)
pos_AHREs_coded$bound_status <- rep(1L, nrow(pos_AHREs_coded))
neg_AHREs_coded$bound_status <- rep(0L, nrow(neg_AHREs_coded))

# Combine pos and neg AHREs into one data frame
all_AHREs_coded <- rbind(pos_AHREs_coded, neg_AHREs_coded)

# Carry out PCA on coded AHRE sites (removing constant middle GCGTG core)
all_AHREs_coded.pca <- prcomp(all_AHREs_coded[ , c(1:28, 49:76)],
                              center = TRUE,
                              scale. = TRUE) 
summary(all_AHREs_coded.pca)

# Scree plot
plot(all_AHREs_coded.pca,type="lines")

# Generate biplot of PCA results
ggp1 <- ggbiplot(all_AHREs_coded.pca, obs.scale = 1, var.scale = 1, 
              groups = as.factor(all_AHREs_coded$bound_status), alpha = 0.5, var.axes = FALSE,
              ellipse = FALSE, circle = FALSE) + 
  scale_color_hue(name = "AHREs", labels = c("unbound", "bound")) 

ggp1


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


#--------
# Prepare AHR  peak area data for regression model

# Read in relevant files
# peak_sums_df <- read.table("AHRpeaks_FDR1_readSums.tab", sep="\t")
# peaks_info <- read.csv("AHRpeaks_features.csv", header = TRUE)
# AHREs_all <- read.csv("AHRpeak_AHREs_19bp_withPeakNum.csv", header = TRUE)
# 
# # Add peak read sums and averages to peaks_info
# colnames(peak_sums_df)[4] <- "read_sums"
# peaks_info$peak_sum <- peak_sums_df$read_sums
# 
# # Get peaks with only one AHRE
# peaks_withSingleAHRE <- peaks_info[peaks_info$num_dres == 1, ]
# 
# # Next, get the AHREs under those peaks
# singletonAHREs <- read.csv("singletonAHREs_info.csv", header = TRUE)


