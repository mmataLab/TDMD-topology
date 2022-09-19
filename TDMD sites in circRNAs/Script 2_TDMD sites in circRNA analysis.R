### Script 2: Analyze the sites obtained through scanMiR in circRNAs
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse", "ggpubr", "viridis")

# Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# ------------------------------------------------------------------------
# Set the working directory and load the necessary tables. 
# ------------------------------------------------------------------------

# Set working directory
setwd("~/Influence of target RNA topology on microRNA stability/Bioinformatics/")

# Load tables (you need to have exported this table with script 1 from the miRNA sponging coefficient folder)
  # Sponging coefficient table
most_sponged_miRNAs <- read_csv("./miRNA sponging coefficient/Results/miRNA sponging coefficient.csv")
  # circRNA expression tables
circRNA_expression_H9 <- read.csv("./miRNA sponging coefficient/circRNA_expression_H9_GSE73325.csv", sep = ";", header = TRUE)
circRNA_expression_FB <- read.csv("./miRNA sponging coefficient/circRNA_expression_FB_GSE73325.csv", sep = ";", header = TRUE)
  # Sites obtained with scanMiR using Script 2
circ_miR_matches_TDMD <- read_csv("./TDMD sites in circRNAs/Results/circ_miR_matches_TDMD_all.csv")

# ------------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------------

# Correct miRNA names at the most sponged miRNA table
most_sponged_miRNAs$miRNAname <- gsub("hsa-miR-137", "hsa-miR-137-3p",
                                      gsub("hsa-miR-217", "hsa-miR-217-5p", 
                                           gsub("hsa-miR-320a", "hsa-miR-320a-3p",
                                                gsub("hsa-miR-375", "hsa-miR-375-3p",
                                                     most_sponged_miRNAs$miRNAname))))

# Merge the circRNA tables
circRNA_expression <- left_join(circRNA_expression_FB, circRNA_expression_H9[,c(1,2,3,8)], by = c("Chromosome", "Start", "End"))
circRNA_expression <- circRNA_expression[,c(1,2,3,4,5,6,7,9,10,11,12,8)]
colnames(circRNA_expression)[11:12] <- c("H9_JReads", "FB_JReads")

# String split the circ-miR matches to get the coordinates of the circRNAs
coord <- as.data.frame(str_split_fixed(circ_miR_matches_TDMD$transcript, pattern = "\\|", n = 4)) 
coord <- as.data.frame(str_split_fixed(coord$V2, pattern = ":", n = 2)) 
coord2 <- as.data.frame(str_split_fixed(coord$V2, pattern = "-", n = 2)) 
coord2$V2 <- gsub("\\+", "", 
                  gsub("-", "", coord2$V2))
colnames(coord) <- c("Chromosome")
colnames(coord2) <- c("Start", "End")

# Add it to the original table
circ_miR_matches_TDMD <- cbind(circ_miR_matches_TDMD[,c(1:2,5:7)], coord[1], coord2)
circ_miR_matches_TDMD$Start <- as.integer(circ_miR_matches_TDMD$Start)
circ_miR_matches_TDMD$End <- as.integer(circ_miR_matches_TDMD$End)

# Now add the circRNA expression data
circ_miR_full <- left_join(circ_miR_matches_TDMD, circRNA_expression)

# Keep only those circRNAs reported in the circRNA expression data (GSE73325)
circ_miR_filtered <- circ_miR_full %>% drop_na(Strand)

# Add the sponging coefficient data
circ_miR_filtered <- left_join(circ_miR_filtered, most_sponged_miRNAs[,c(1:4,7)], by = c("miRNA" = "miRNAname"))

# See how many sites each miRNA has out of the miRNA list
    # Only consider circRNAs that are expressed before differentiation
sites_per_miR <- subset(circ_miR_full, H9_JReads > 0) %>% group_by(miRNA) %>% summarise(sites = n())
sites_per_miR <- left_join(most_sponged_miRNAs[,c(1:4,7)], sites_per_miR, by = c("miRNAname" = "miRNA"))
sites_per_miR[is.na(sites_per_miR)] <- 0

# Order for figures
sites_per_miR$quartile <- factor(sites_per_miR$quartile, levels = c("- sponged",
                                                                    "+ sponged",
                                                                    "++ sponged",
                                                                    "+++ sponged"))

# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

### Compare the amount of sites each miRNA belonging to each quartile has
gts <- ggviolin(sites_per_miR,
         x = "quartile",
         y = "sites",
         fill = "quartile",
         add = c("jitter"),
         add.params = list(size = 0.4),
         alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, guide = "none", option = "D", direction = -1) +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("- sponged", "+ sponged"),
                                        c("- sponged", "++ sponged"),
                                        c("- sponged", "+++ sponged")),
                     tip.length = 0,
                     step.increase = 0.15,
                     label = "p.signif",
                     size = 4) +
  xlab("") + ylab("Predicted TDMD-like sites") +
  ylim(-1, 10) +
  theme_pubr(legend = "none",
             x.text.angle = 45,
             base_family = "Myriad Pro", base_size = 18) + 
  theme(plot.margin = margin(c(40,20,0,10)))
gts

# ------------------------------------------------------------------------
# Export table with TDMD-like site number per miRNA
# ------------------------------------------------------------------------

write_csv(sites_per_miR, "./TDMD sites in circRNAs/Results/TDMD sites on circs per miR.csv")
