#
# Code to compute shape features of AHREs and write to file for classification
#
# Sudin Bhattacharya 5/20/18

library(DNAshapeR)
library(Biostrings)

# Read in AHREs and generate fasta file for input to DNAshapeR
pos_AHREs <- read.csv("singletonAHREs_open.csv", header = TRUE)
neg_AHREs <- read.csv("negative_AHREs.csv", header = TRUE)

# Write pos AHREs to file in fasta format
counter = 1
for (AHRE in pos_AHREs$seq) {
  str1 <- paste(">posAHRE_", counter, sep = "")
  write(str1, file = "pos_AHREs_fasta.fa", append = TRUE)
  write(AHRE, file = "pos_AHREs_fasta.fa", append = TRUE)
  counter = counter + 1
}

# Write neg AHREs to file in fasta format
counter = 1
for (AHRE in neg_AHREs$seq) {
  str1 <- paste(">negAHRE_", counter, sep = "")
  write(str1, file = "neg_AHREs_fasta.fa", append = TRUE)
  write(AHRE, file = "neg_AHREs_fasta.fa", append = TRUE)
  counter = counter + 1
}

# Read in DRE sequences in fasta format for input to DNAshapeR for next line
# ** note that for the input commands below to work, fasta files have to be in 
# /home/sudin/R/x86_64-pc-linux-gnu-library/3.3/DNAshapeR/extdata/  (Linux)
# or /Library/Frameworks/R.framework/Resources/library/DNAshapeR/extdata/ (MacOS)
# or equivalent
pos_AHREs_sh0 <- system.file("extdata", "pos_AHREs_fasta.fa", package = "DNAshapeR")
neg_AHREs_sh0 <- system.file("extdata", "neg_AHREs_fasta.fa", package = "DNAshapeR")

pos_AHREs_sh1 <- getShape(pos_AHREs_sh0)
neg_AHREs_sh1 <- getShape(neg_AHREs_sh0)

# Visualize DNA shape predictions for pos and neg AHREs
plotShape(pos_AHREs_sh1$MGW)
plotShape(neg_AHREs_sh1$MGW)

plotShape(pos_AHREs_sh1$ProT)
plotShape(neg_AHREs_sh1$ProT)

plotShape(pos_AHREs_sh1$Roll)
plotShape(neg_AHREs_sh1$Roll)

plotShape(pos_AHREs_sh1$HelT)
plotShape(neg_AHREs_sh1$HelT)

# Generate shape feature vectors for positive and negative DREs
featureType <- c("1-MGW", "1-ProT", "1-Roll", "1-HelT")
featureVector_pos <- encodeSeqShape(pos_AHREs_sh0, pos_AHREs_sh1, featureType)
featureVector_neg <- encodeSeqShape(neg_AHREs_sh0, neg_AHREs_sh1, featureType)

# Write shape feature vectors to file
write.csv(featureVector_pos, file = "AHREs_pos_shape_features.csv", row.names = FALSE)
write.csv(featureVector_neg, file = "AHREs_neg_shape_features.csv", row.names = FALSE)



