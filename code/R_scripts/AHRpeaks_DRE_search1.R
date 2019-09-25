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
# Create an empty data frame to store long (19-bp) site, location and MS score
# for all putative DREs in sequence
longSites <- data.frame(seq = character(), 
                        chrom = character(),
                        start_coord = integer(),
                        end_coord = integer(),
                        siteMS_score = numeric(),
                        stringsAsFactors=FALSE)


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
  
  if (numSites_in_myseq > 0) {
    cat("* Found ", numSites_in_myseq, " putative DREs in peak # ", i, " \n")
    
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
      temp_df <- data.frame(seq = long_site_string, 
                            chrom = peak_chrom2,
                            start_coord = start_long_site + peak_start, 
                            end_coord = end_long_site + peak_start, 
                            siteMS_score = MS_score)
      longSites <- rbind(longSites, temp_df) 
      
    }  # end of loop over DREs
    
  }
  
}  # end loop over peak seqs


##----------
# Save long-site seqs and scores to CSV file
write.csv(longSites, file = "AHRpeak_AHREs_19bp.csv", row.names = FALSE)


##----------
# Create bed file from long-site coords

# Create new data frame from longSites with only necessary fields
longSites2 <- data.frame(chr = longSites$chrom,
                         start = longSites$start_coord,
                         end = longSites$end_coord,
                         score = longSites$siteMS_score)
longSites2$score <- format(longSites2$score, digits = 3)
head(longSites2)

# Write to bed file
f <- file("AHRpeak_AHREs_19bp.bed", "w") 
line = 'track name="AHREs"  description="AHREs under AHR peaks"  useScore=1'
writeLines(line, con = f)
close(f)
write.table(longSites2, file="AHRpeak_AHREs_19bp.bed",
            quote=F, sep="\t", row.names=F, col.names=F, append = T)

