### Script 2b: Analyze the sites obtained through scanMiR in linRNAs
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse", "ggpubr", "viridis", "emmeans", "rstatix", "Homo.sapiens")

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
most_sponged_miRNAs <- read_csv("./miRNA sponging coefficient/Results/miRNA sponging coefficient_2.0.csv")
  # mRNA expression tables
mRNA_expression <- read_csv2("./miRNA sponging coefficient/mRNA_expression_GSE73325.csv", col_names = TRUE)
  # Sites obtained with scanMiR using Script 1b
lin_miR_matches_TDMD <- read_rds("./TDMD sites in circRNAs/Results/hg19_TDMD_spongedMiRs_UTR_ORF.rds")
    # Make it a data frame and keep only 3'UTR sites
lin_miR_matches_TDMD_df <- as.data.frame(lin_miR_matches_TDMD) %>% filter(ORF == FALSE)

# ------------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------------

# Correct miRNA names at the most sponged miRNA table
most_sponged_miRNAs$miRNAname <- gsub("hsa-miR-137", "hsa-miR-137-3p",
                                      gsub("hsa-miR-217", "hsa-miR-217-5p", 
                                           gsub("hsa-miR-320a", "hsa-miR-320a-3p",
                                                gsub("hsa-miR-375", "hsa-miR-375-3p",
                                                     most_sponged_miRNAs$miRNAname))))

# Convert from UCSC gene ID to Symbols to add expression
annots <- select(Homo.sapiens, keys = as.character(lin_miR_matches_TDMD_df$seqnames), "SYMBOL", "TXNAME")
Symbol <- as.vector(annots$SYMBOL)
lin_miR_matches_TDMD_df <- cbind(lin_miR_matches_TDMD_df, Symbol)

# Add the circRNA expression data to the TDMD matches table
lin_miR_full <- left_join(lin_miR_matches_TDMD_df, mRNA_expression)

# Keep only those circRNAs reported in the circRNA expression data
lin_miR_filtered <- lin_miR_full %>% drop_na(H9_RPKM)

# Add the sponging coefficient data
lin_miR_filtered <- left_join(lin_miR_filtered, most_sponged_miRNAs[,c(1:4,7)], by = c("miRNA" = "miRNAname"))

# See how many sites each miRNA has out of the miRNA list
    # Only consider linRNAs that are expressed before differentiation
sites_per_miR <- subset(lin_miR_full, H9_RPKM > 5) %>% group_by(miRNA) %>% summarise(sites = n())
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
    # Test for normality first
shapiro.test(sites_per_miR$sites) # If you run it, you'll see that the distribution is not normal.
    # Do a general linear model, not assuming normality. 
glmodel <- glm(data = sites_per_miR,
               formula = sites ~ quartile)
        # Follow up with estimated marginal means -emmeans- and contrasts.
stat.test <- emmeans(glmodel, pairwise ~ quartile)
stat.c <- cbind(data.frame(group1 = c("- sponged", "- sponged", "- sponged", "+ sponged", "+ sponged", "++ sponged"),
                           group2 = c("+ sponged", "++ sponged", "+++ sponged", "++ sponged", "+++ sponged", "+++ sponged")), as_tibble(stat.test$contrasts)) %>% add_significance()
stat.c <- stat.c[1:3,]

gts <- ggviolin(sites_per_miR,
         x = "quartile",
         y = "sites",
         fill = "quartile",
         add = c("jitter"),
         add.params = list(size = 0.4),
         alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, guide = "none", option = "D", direction = -1) +
  xlab("") + ylab("Predicted TDMD-like sites\non linear RNAs") +
  theme_pubr(legend = "none",
             x.text.angle = 45,
             base_family = "Myriad Pro", base_size = 18) + 
  theme(plot.margin = margin(c(40,20,0,10)))
gts

# ------------------------------------------------------------------------
# Export table with TDMD-like site number per miRNA
# ------------------------------------------------------------------------

write_csv(sites_per_miR, "./TDMD sites in circRNAs/Results/TDMD sites on lins per miR_2.0.csv")
