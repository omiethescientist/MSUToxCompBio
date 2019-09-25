# 
# Code to scan AHR ChIP peak fasta seqs, search for *all* potential DREs, and list them
# 
#   Sudin Bhattacharya 4/25/18

library(seqinr)
library(Biostrings)
library(TFBSTools)
library(seqLogo)
library(biomaRt)
library(JASPAR2016)
library(reshape2)
library(ggplot2)
library(nucleR)

# Get the AhR PFM from JASPAR 2016 and convert to PWM
Ahr_pfm <- getMatrixByName(JASPAR2016, "Ahr::Arnt")
Ahr_pwm <- toPWM(Ahr_pfm)

# Modify the AhR PFM matrix to include just DRE core "GCGTG"
Ahr_pfm_mod <- Ahr_pfm
Matrix(Ahr_pfm_mod) <- matrix(c( 0L,  0L,  0L,  0L,  0L, 
                                 0L, 24L,  0L,  0L,  0L,
                                 24L,  0L, 24L,  0L, 24L,
                                 0L,  0L,  0L, 24L,  0L), 
                              byrow=TRUE, nrow=4, 
                              dimnames=list(c("A", "C", "G", "T")))
Ahr_pwm_mod <- toPWM(Ahr_pfm_mod)

# Consensus index vector Ci based on 13 bona fide DREs in Dere et al, CRT 2011
Ci <- c(33,37,21,24,37,57,42,100,100,100,100,100,57,37,42,28,47,34,24)

# PWM of relative base frequencies based on the 13 bona fide DREs
A_freq <- c(8,15,39,23,0,0,15,0,0,0,0,0,54,15,69,54,0,31,8)
C_freq <- c(62,46,23,23,54,8,8,0,100,0,0,0,46,39,8,8,31,46,23)
G_freq <- c(15,39,31,46,23,15,8,100,0,100,0,100,0,46,8,23,62,23,23)
T_freq <- c(15,0,8,8,23,77,69,0,0,0,100,0,0,0,15,15,8,0,46)

length_cons_mat <- length(Ci)  # length of consensus matrix

# Create max_score vector
max_score <- vector()
for (i in 1:length_cons_mat) {
  max_score[i] <- max(A_freq[i], C_freq[i], G_freq[i], T_freq[i])
}

# Set threshold MS score for DREs (value from Dere et al, CRT 2011)
#MS_score_threshold = 0.8473


##-----------
# Create an empty data frame to store long (19-bp) site, location, MS score 
# and associated peak no. for all putative DREs in sequence
longSites2 <- data.frame(seq = character(), 
                        chrom = character(),
                        start_coord = integer(),
                        end_coord = integer(),
                        siteMS_score = numeric(),
                        peak_assoc_with = integer(),
                        stringsAsFactors = FALSE)

peaks_info <- data.frame(id = integer(),
                         num_dres = integer(),
                         density_dres = numeric(),
                         stringsAsFactors = FALSE)


##----------
# Read in AHR peak seqs fasta file
peak_seqs <- read.fasta("AHRpeaks_seqs.fa")
cat("--- No. of AHR peaks found = ", length(peak_seqs), "\n")

##----------
# Loop over all peak sequences read in
for (i in 1:length(peak_seqs)) {
  # Read in whole sequence
  mySeq <- peak_seqs[[i]]  
  
  peak_annot <- attr(mySeq, "Annot")  # Get the annotation of the peak
  peak_chrom1 <- strsplit(peak_annot, ":")[[1]][1]  # Get chromosome # of peak
  peak_chrom2 <- strsplit(peak_chrom1, ">")[[1]][2]
  
  peak_range <- strsplit(peak_annot, ":")[[1]][2]  # Get chromosomal range of peak
  peak_start <- as.numeric(strsplit(peak_range, "-")[[1]][1])  # Get start coord of peak
  
  
  # Convert sequence to DNA string
  mySeqDNAString <- DNAString(c2s(mySeq))
  
  # Search sequence for exct match to *modified* AhR PWM     
  siteSet <- searchSeq(Ahr_pwm_mod, mySeqDNAString, 
                       seqname = "seq1", min.score = "100%", strand = "*")
  sites_in_myseq <- writeGFF3(siteSet)
  
  # If binding sites found, add all instances to dataframe of sites
  numSites_in_myseq <- dim(sites_in_myseq)[1]
  cat("* Found ", numSites_in_myseq, " putative DREs in peak # ", i, " \n")
  
  if (numSites_in_myseq > 0) {
    # initialize list of DRE start and end locations
    starts = numeric()
    ends = numeric()
    
    # Loop over DREs
    for (j in 1:numSites_in_myseq) {
      # Find start and end coords of site
      start <- sites_in_myseq$start[j]
      end <- sites_in_myseq$end[j]
      
      # Ignore sites too close to end of sequence string
      if ((start <= 7) || (end >= length(mySeqDNAString) - 7)) {
        next
      }
      
      if (sites_in_myseq$strand[j] == "+") {
        long_site <- mySeqDNAString[(start - 7):(end + 7)]
      } else {
        long_site <- reverseComplement(mySeqDNAString[(start - 7):(end + 7)])
      }
      
      # Compute MS score of long_site
      numer <- 0  # Initialize numerator to zero
      denom <- 0  # Initialize denominator to zero
      for (k in 1:length(long_site)) {
        if (long_site[k] == DNAString("A")) {
          base_score <- A_freq[k]
        } else if (long_site[k] == DNAString("C")) {
          base_score <- C_freq[k]
        } else if (long_site[k] == DNAString("G")) {
          base_score <- G_freq[k]
        } else if (long_site[k] == DNAString("T")) {
          base_score <- T_freq[k]
        }
        
        numer <- numer + Ci[k]*base_score
        denom <- denom + Ci[k]*max_score[k]
      }
      
      MS_score <- numer/denom
      #cat("MS_score = ", MS_score, "\n")
      
      # Add long_site info and MS score to longSites dataframe
      #cat("DRE no. = ", j, "\n")
      long_site_string <- toString(long_site)
      start_long_site = (start - 7)
      end_long_site = (end + 7) 
      
      # Add to list of starts and ends for current peak
      starts <- c(starts, start_long_site)
      ends <- c(ends, end_long_site)
      
      # Create temporary data frame of long sites
      temp_df <- data.frame(seq = long_site_string, 
                            chrom = peak_chrom2,
                            start_coord = start_long_site + peak_start, 
                            end_coord = end_long_site + peak_start, 
                            siteMS_score = MS_score,
                            peak_assoc_with = i)
      
      # Append to list of long sites
      longSites2 <- rbind(longSites2, temp_df) 
      
    }  # end of loop over DREs
    
    # Calculate spread of DREs under peak
    spread_dres = max(ends) - min(starts)
    
    # Create temporary data frame of AHR peaks
    if (spread_dres != 0) {
      temp_df_peaks <- data.frame(id = i,
                                  num_dres = numSites_in_myseq,
                                  density_dres = numSites_in_myseq / spread_dres)
    } else {
      cat("**  Spread of DREs = 0!! \n")
    }
   
    
  } else {
    # No. of DREs in peak = 0
    temp_df_peaks <- data.frame(id = i,
                                num_dres = 0,
                                density_dres = 0)
  }
  
  # Append to list of peaks
  peaks_info <- rbind(peaks_info, temp_df_peaks) 
  
  
}  # end loop over peak seqs


##----------
# Save long-site seqs, scores and peak nos. to CSV file
write.csv(longSites2, file = "AHRpeak_AHREs_19bp_withPeakNum.csv", row.names = FALSE)

# Savepeak info. to CSV file
write.csv(peaks_info, file = "AHRpeaks_features.csv", row.names = FALSE)

# Create histogram of no. of DREs among AHR peaks
hist(peaks_info$num_dres, 
     main = "Distribution of DREs among AHR peaks",
     xlab = "No. of DREs", ylab = "No. of AHR peaks",
     col = "cyan",
     border = "blue")

##----------
# Create histogram of no. of DREs among AHR peaks with ggplot2
ggp <- ggplot(data = peaks_info, aes(x = num_dres)) + 
  geom_histogram(color="darkblue", fill="lightblue", binwidth = 1) + 
  ggtitle("Distribution of AHREs among AHR ChIP peaks") + 
  labs(x = "Number of AHREs per peak", y = "Number of AHR peaks") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 14), axis.text=element_text(size = 12))
ggp

##----------
# System calls to Unix
# Remove header lines from AHR peaks bed file
system("tail -n +3 AHRpeaks_FDR1.bed > AHRpeaks_FDR1_noHeader.bed")

# Get peak read counts from AHR peak regions
system("bigWigAverageOverBed TCDD_AHR_24h_treat.bigwig AHRpeaks_FDR1_noHeader.bed AHRpeaks_FDR1_readSums.tab")

peak_sums_df <- read.table("AHRpeaks_FDR1_readSums.tab", sep="\t")
colnames(peak_sums_df)[4] <- "read_sums"
colnames(peak_sums_df)[5] <- "read_avg"

# Add peak read sums and averages to peaks_info
peaks_info$peak_sum <- peak_sums_df$read_sums
peaks_info$peak_avg <- peak_sums_df$read_avg


##----------
# Plot scatter plot of peak sums

# Plot scatter plot of peak sums vs no. DREs
ggp2 <- ggplot(data = peaks_info, aes(x = num_dres, y = peak_sum)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Number of AHREs per AHR peak", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp2

# Plot scatter plot of peak sums vs density of DREs
ggp3 <- ggplot(data = peaks_info, aes(x = density_dres, y = peak_sum)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Density of AHREs under AHR peak", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp3

# --Filter out peaks with no AHREs before plotting histograms etc.
peaks_info_filt <- peaks_info[peaks_info$num_dres > 0, ]

# Plot scatter plot of peak sums vs no. DREs *after filtering out peaks with no DREs*
ggp4 <- ggplot(data = peaks_info_filt, aes(x = num_dres, y = peak_sum)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "Number of AHREs under AHR peak", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp4

# Plot scatter plot of peak sums vs density of DREs *after filtering out peaks with no DREs*
ggp5 <- ggplot(data = peaks_info_filt, aes(x = density_dres, y = peak_sum)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10(limits = c(10, 3.5e5)) +
  labs(x = "Density of AHREs under AHR peak", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp5

#--
# Plot *boxplot* of AHR peak areas grouped by no. DREs per peak
ggp4a <- ggplot(data = peaks_info_filt, aes(x = as.factor(num_dres), y = peak_sum)) + 
  geom_boxplot(fill='lightblue', color="black") + scale_y_log10() +
  labs(x = "Number of AHREs under AHR peak", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp4a


##----------
# Plot scatter plot of peak averages

# Plot scatter plot of peak averages vs no. DREs
ggp6 <- ggplot(data = peaks_info, aes(x = num_dres, y = peak_avg)) + 
  geom_point(alpha = 1/3) + 
  labs(x = "Number of AHREs per AHR peak", y = "Averaged area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp6

# Plot scatter plot of peak averages vs no. DREs *after filtering out peaks with no DREs*
ggp7 <- ggplot(data = peaks_info_filt, aes(x = num_dres, y = peak_avg)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "Number of AHREs under AHR peak", y = "Averaged area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp7

# Plot scatter plot of peak averages vs density of DREs *after filtering out peaks with no DREs*
ggp8 <- ggplot(data = peaks_info_filt, aes(x = density_dres, y = peak_avg)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "Density of AHREs under AHR peak", y = "Averaged area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp8

#--- Try this:
# Plot scatter plot of peak averages vs peak areas *after filtering out peaks with no DREs*
ggp9 <- ggplot(data = peaks_info_filt, aes(x = peak_sum, y = peak_avg)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_x_log10() +  scale_y_log10() +
  labs(x = "Area under AHR peak", y = "Averaged area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp9


#------
# Get singleton AHREs (AHREs under peaks with only one AHRE)

# First, get peaks with only one AHRE
peaks_withSingleAHRE <- peaks_info[peaks_info$num_dres == 1, ]

# Next, get the AHREs under those peaks
singletonAHREs <- longSites2[longSites2$peak_assoc_with %in% peaks_withSingleAHRE$id, ]

# -- NOTE: 837 singleton_AHREs, but 844 peaks_withSingleAHRE => some AHREs shared by peaks?

#--- See if there is any correlation between AHRE MS_score and AHR peak area?

# Filter peaks_withSingleAHRE such that one peak assocaited with each singleton_AHRE
peaks_withSingleAHRE_filt <- peaks_withSingleAHRE[peaks_withSingleAHRE$id %in% 
                                                    singletonAHREs$peak_assoc_with, ]

# Add peak_sum column to singleton_AHRE
singletonAHREs$peak_sum <- peaks_withSingleAHRE_filt$peak_sum

# Draw scatter-plot of AHR peak area vs. AHRE MS score
# Plot scatter plot of peak sums vs no. DREs *after filtering out peaks with no DREs*
ggp10 <- ggplot(data = singletonAHREs, aes(x = siteMS_score, y = peak_sum)) + 
  geom_point(colour = "darkblue", alpha = 1/3) + scale_y_log10() +
  labs(x = "AHRE MS score", y = "Area under AHR peak") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0, size = 14)) +
  theme(axis.title = element_text(size = 12), axis.text=element_text(size = 12))
ggp10

# --> Little correlation between singleton_AHRE MS_score and AHR peak area

# Save singletonAHREs dataframe to csv file
write.csv(singletonAHREs, file = "singletonAHREs_info.csv", row.names = FALSE)



