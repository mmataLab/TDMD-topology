### Script 1: General landscape of miRNA expression and miRNA sponging
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse", "viridis", "fuzzyjoin", "ggpubr", "ggrepel", "emmeans", "rstatix")

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
# Set the working directory, where the tables of circRNA, miRNA and mRNA expression are; and load them. 
# ------------------------------------------------------------------------

# Set working directory
setwd("~/Influence of target RNA topology on microRNA stability/Bioinformatics/")

# Load tables
mRNA_expression <- read_csv2("./miRNA sponging coefficient/mRNA_expression_GSE73325.csv", col_names = TRUE)
miRNA_expression <- read_csv2("./miRNA sponging coefficient/miRNA_expression_GSE56152_GSE63709.csv", col_names = TRUE)
circRNA_expression_H9 <- read.csv("./miRNA sponging coefficient/circRNA_expression_H9_GSE73325.csv", sep = ";", header = TRUE)
circRNA_expression_FB <- read.csv("./miRNA sponging coefficient/circRNA_expression_FB_GSE73325.csv", sep = ";", header = TRUE)
miRNA_sites <- read_tsv("./miRNA sponging coefficient/starBaseV3_hg19_CLIP-seq_all_all.txt")

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Starbase has many sites duplicated, due to the way it's annotated; we need to keep only one to avoid over-estimation of the number of interactions
miRNA_sites <- miRNA_sites %>% group_by(miRNAid, miRNAname, chromosome, start, end, strand) %>% slice(1)

# Merge the circRNA expression tables
circRNA_expression <- left_join(circRNA_expression_FB, circRNA_expression_H9[,c(1,2,3,8)], by = c("Chromosome", "Start", "End"))
circRNA_expression <- circRNA_expression[,c(1,2,3,4,5,6,7,9,10,11,12,8)]
colnames(circRNA_expression)[11:12] <- c("H9_JReads", "FB_JReads")

# Add the miRNA expression data to the sites provided by STARBASE
miRNA_circRNA <- left_join(miRNA_sites, miRNA_expression, by = c("miRNAname" = "ID"))

# Join circRNA expression table with the miRNA site table using Fuzzyjoin, which allows joining using genomic coordinates
miRNA_circRNA_Overlaps <- genome_left_join(circRNA_expression, miRNA_circRNA, by = c("Chromosome" = 'chromosome', 'Start' = 'start', 'End' = 'end'))
miRNA_circRNA_Overlaps <- miRNA_circRNA_Overlaps[!is.na(miRNA_circRNA_Overlaps$miRNAname),]  # Remove NAs
miRNA_circRNA_Overlaps <- miRNA_circRNA_Overlaps[!is.na(miRNA_circRNA_Overlaps$FB_RPM),]  # Remove NAs

# Calculate log2FoldChanges for circRNAs, mRNAs and miRNAs
mRNA_expression$log2FoldChange_mRNA <- log2(mRNA_expression$FB_RPKM/mRNA_expression$H9_RPKM) # Add a column with log2 fold changes of mRNAs
miRNA_circRNA_Overlaps <- left_join(miRNA_circRNA_Overlaps, mRNA_expression, by = c("Gene.symbol" = "Symbol"))
miRNA_circRNA_Overlaps[is.na(miRNA_circRNA_Overlaps)] <- 0  # Swap NAs for 0
miRNA_circRNA_Overlaps$H9_JReads <- miRNA_circRNA_Overlaps$H9_JReads + 1  # Add pseudocounts to be able to calculate log2 fold changes of circRNAs
miRNA_circRNA_Overlaps$FB_JReads <- miRNA_circRNA_Overlaps$FB_JReads + 1  # Add pseudocounts to be able to calculate log2 fold changes of circRNAs
miRNA_circRNA_Overlaps$H9_RPM <- miRNA_circRNA_Overlaps$H9_RPM + 1  # Add pseudocounts to be able to calculate log2 fold changes of miRNAs
miRNA_circRNA_Overlaps$FB_RPM <- miRNA_circRNA_Overlaps$FB_RPM + 1  # Add pseudocounts to be able to calculate log2 fold changes of miRNAs
miRNA_circRNA_Overlaps$log2FoldChange_circ <- log2(miRNA_circRNA_Overlaps$FB_JReads/miRNA_circRNA_Overlaps$H9_JReads) # Add a column with log2 fold changes of circRNAs
miRNA_circRNA_Overlaps$log2FoldChange_miR <- log2(miRNA_circRNA_Overlaps$FB_RPM/miRNA_circRNA_Overlaps$H9_RPM) # Add a column with log2 fold changes of miRNAs

# Count the number of sites each miRNA has on each circRNA
Sites_per_circ <- miRNA_circRNA_Overlaps %>% group_by(Chromosome, Start, End, Strand, geneName, miRNAname) %>% summarise(count=n())
colnames(Sites_per_circ)[7] <- c("Sites_per_circ") 
  # Add it to the previous table
miRNA_circRNA_Overlaps <- left_join(miRNA_circRNA_Overlaps, Sites_per_circ, by = c("Chromosome", "Start", "End", "Strand", "geneName", "miRNAname"))

# Calculate a relative amount of target sites based on how much the circRNA is expressed before diff. and how many sites it has for a given miRNA
miRNA_circRNA_Overlaps$Relative_TargetSites <- miRNA_circRNA_Overlaps$H9_JReads * miRNA_circRNA_Overlaps$Sites_per_circ

# Calculate the total amount of sites each miRNA has overall the circRNAs
Total_Targetsites <- miRNA_circRNA_Overlaps %>% group_by(Chromosome, Start, End, Strand, geneName, miRNAname, Sites_per_circ, Relative_TargetSites) %>% summarise()
Total_Targetsites <- Total_Targetsites %>% group_by(miRNAname) %>% summarise(Total_TargetSites = sum(Relative_TargetSites), Total_Sites_Per_miR = sum(Sites_per_circ))
  # Add it to the previous table
miRNA_circRNA_Overlaps <- left_join(miRNA_circRNA_Overlaps, Total_Targetsites, by = "miRNAname")

# Calculate how much a given circle influences relative to the total sites that miRNA has
miRNA_circRNA_Overlaps$Ratio_TargetSites <- miRNA_circRNA_Overlaps$Relative_TargetSites / miRNA_circRNA_Overlaps$Total_TargetSites

# Summarise and calculate sponging coefficient 
most_sponged_miRNAs <- miRNA_circRNA_Overlaps %>% group_by(miRNAname, H9_RPM, FB_RPM, log2FoldChange_miR, Total_TargetSites) %>% summarise()
  # The sponging coefficient comes from dividing the total amount of available sites by the expression of the miRNA
most_sponged_miRNAs$Sponging_Coefficient <- most_sponged_miRNAs$Total_TargetSites / most_sponged_miRNAs$H9_RPM

# Split all miRNAs in quartiles, separated by their sponging coefficient
most_sponged_miRNAs <- within(most_sponged_miRNAs, quartile <- as.integer(cut(Sponging_Coefficient, quantile(Sponging_Coefficient, type = 8), include.lowest=TRUE)))
most_sponged_miRNAs$quartile <- gsub(pattern = "1", replacement = "- sponged", most_sponged_miRNAs$quartile)
most_sponged_miRNAs$quartile <- gsub(pattern = "2", replacement = "+ sponged", most_sponged_miRNAs$quartile)
most_sponged_miRNAs$quartile <- gsub(pattern = "3", replacement = "++ sponged", most_sponged_miRNAs$quartile)
most_sponged_miRNAs$quartile <- gsub(pattern = "4", replacement = "+++ sponged", most_sponged_miRNAs$quartile)

# ------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------

# General landscape of microRNA expression and the number of sites each miRNA has available

ggl <- ggplot(data = most_sponged_miRNAs, aes(x = Total_TargetSites, y = log2FoldChange_miR)) +
  stat_cor(method = "pearson") +
  geom_point(color = "grey30", alpha = 0.8) + ylim(-20,20) +
  xlab("miRNA-specific 'Effective' sites on all circRNAs") + ylab("miRNA Fold change (log2)") + 
  geom_text_repel(data = subset(most_sponged_miRNAs, Total_TargetSites > 1000),
                  aes(x = log10(H9_RPM), y = Total_TargetSites),
                  label = subset(most_sponged_miRNAs, Total_TargetSites > 1000)$miRNAname,
                  size = 3, box.padding = 1, alpha = 0.8, segment.size = 0.2, segment.alpha = 0.7, color = "black") +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) +
  theme(legend.text = element_text(size = 9),
        legend.title = element_text(size = 14),
        legend.key.width = unit(0.6, 'cm'),
        plot.margin = margin(c(10,30,10,10))) +
  stat_cor(label.x = 3, label.y = -25) 
ggl

## microRNAs and sponging coefficient

gspo <- ggplot(data = most_sponged_miRNAs, aes(x = log10(Sponging_Coefficient), y = log2FoldChange_miR)) +
  geom_point(alpha = 0.8, aes(color = quartile, size = H9_RPM)) +
  stat_cor(method = "pearson") +
  scale_color_viridis(option = "D", direction = -1, discrete = TRUE, guide = "none") +
  scale_size(name = "miRNA expression\nPre-differentiation", limits = c(1,150000)) +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) + 
  geom_text_repel(data = subset(most_sponged_miRNAs, miRNAname == "hsa-miR-7-5p"),
                  aes(x = log10(Sponging_Coefficient), y = log2FoldChange_miR),
                  label = subset(most_sponged_miRNAs, miRNAname == "hsa-miR-7-5p")$miRNAname,
                  size = 3, alpha = 0.9, segment.size = 0.2, segment.alpha = 0.7, color = "black",
                  nudge_x = 1, nudge_y = -10, arrow = arrow(angle = 10, length = unit(0.1, "inches"))) +
  theme(legend.text = element_text(size = 9),
        legend.title = element_text(size = 14),
        legend.key.width = unit(0.02, 'cm')) +
  xlab("'Sponging' suffered coefficient (log10)") + ylab("miRNA Fold change (log2)") +
  ylim(-20,20)
gspo

### Boxplot of the quartiles of microRNAs separated by sponging coefficient
    # Test for normality first
shapiro.test(most_sponged_miRNAs$log2FoldChange_miR) # If you run it, you'll see that the distribution is not normal.
    # Do a general linear model, not assuming normality. 
glmodel <- glm(data = most_sponged_miRNAs,
             formula = log2FoldChange_miR ~ quartile)
        # Follow up with estimated marginal means -emmeans- and contrasts.
stat.test <- emmeans(glmodel, pairwise ~ quartile)
stat.c <- cbind(data.frame(group1 = c("- sponged", "- sponged", "- sponged", "+ sponged", "+ sponged", "++ sponged"),
                           group2 = c("+ sponged", "++ sponged", "+++ sponged", "++ sponged", "+++ sponged", "+++ sponged")), as_tibble(stat.test$contrasts)) %>% add_significance()
stat.c <- stat.c[1:3,]

    #plot
gmsmbmi <- ggplot(data = most_sponged_miRNAs, aes(x = quartile, y = log2FoldChange_miR)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE, aes(fill = quartile)) +
  geom_jitter(size = 0.1, alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, guide = "none", option = "D", direction = -1) +
  xlab("") + ylab("miRNA Fold change (log2)") + 
  theme_pubr(base_family = "Myriad Pro", base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(stat.c,
                     label = "p.value.signif", 
                     step_increase = 5, 
                     y.position = c(10,13,16),
                     tip.length = 0) +
  ylim(-20,20) +
  theme(plot.margin = margin(c(50,70,0,10)))
gmsmbmi

# ------------------------------------------------------------------------
# Export the table of miRNAs arranged by sponging coefficient and the full interactions table
# ------------------------------------------------------------------------

write_csv(most_sponged_miRNAs, "./miRNA sponging coefficient/Results/miRNA sponging coefficient.csv")
write_csv(miRNA_circRNA_Overlaps, "./miRNA sponging coefficient/Results/miR_Cir_interactions_fulltable.csv")



