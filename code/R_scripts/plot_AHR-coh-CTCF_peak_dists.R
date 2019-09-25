#
# Code to read in AHR-cohesin-CTCF peak separation distances calculated by David
# and create histograms / box-violin plots
#
# Sudin Bhattacharya 5/29/18

library(ggplot2)

# Read in AHR-coh co-bound peak separation dists
AHR_coh_peak_dists <- read.csv("DavidF_motif_anls/AHRpeaks_RAD21-coh-peaks_overlap_distances.csv", 
                               header = TRUE)

# Read in AHR-coh-CTCF co-bound peak separation dists
AHR_coh_CTCF_peak_dists <- read.csv("DavidF_motif_anls/AHRpeaks_CTCFpeaks_RAD21-coh-peaks_overlap_distances.csv", 
                               header = TRUE)

# --- Make individual data frames (168 hr.) ---
AHR_coh_dist_df <- data.frame(group = "AHR-cohesin", 
                              dist = AHR_coh_peak_dists$distance)
AHR_coh_CTCF_dist_df <- data.frame(group = "AHR-cohesin-CTCF", 
                                   dist = AHR_coh_CTCF_peak_dists$distance)

# Combine into single long data frame
peaks_dist_df <- rbind(AHR_coh_dist_df, AHR_coh_CTCF_dist_df)

# Boxplot using ggplot2
gp1 <- ggplot(peaks_dist_df, aes(x=group, y=dist, fill=group)) + 
  geom_boxplot(width=0.3) + 
  scale_y_log10() +
  xlab("Co-bound ChIP peaks") + ylab("Inter-peak distance") + 
  ggtitle("") +
  scale_fill_brewer(palette = "Pastel1") + theme(legend.position="none") + 
  theme(axis.text.x = element_text(size=11, color = "blue"), 
        axis.text.y = element_text(size=11))

gp1

# Boxp and violin lot using ggplot2
gp2 <- ggplot(peaks_dist_df, aes(x=group, y=dist, fill=group)) + 
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=.2, aes(fill = group)) + 
  scale_y_log10() +
  xlab("Co-bound ChIP peaks") + ylab("Inter-peak distance") + 
  ggtitle("") +
  scale_fill_brewer(palette = "Pastel1") + theme(legend.position="none")

gp2
