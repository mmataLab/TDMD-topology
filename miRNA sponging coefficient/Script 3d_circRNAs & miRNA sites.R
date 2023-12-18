### Script 3: circRNAs arranged according to the amount of sites they potentially offer to miRNAs
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
miRNA_sponging_coefficients <- read_csv("./miRNA sponging coefficient/Results/miRNA sponging coefficient_2.0.csv")

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

# Summarise the total amount of sites each circRNA has
most_sponging_circRNAs <- miRNA_circRNA_Overlaps %>% 
  group_by(Chromosome, Start, End, Strand,
           geneName, miRNAname) %>%
  slice(1) %>%
  mutate(Sites_per_circ_Rel = ((H9_RPM+FB_RPM)/2)*Sites_per_circ) %>%
  ungroup() %>%
  group_by(Chromosome, Start, End, Strand, H9_JReads, FB_JReads, geneName, H9_RPKM, FB_RPKM, log2FoldChange_mRNA, log2FoldChange_circ) %>% 
  summarise(Total_Sites_Per_circ = sum(Sites_per_circ_Rel))

# Split the circRNAs into quartiles according to the amount of effective sites
# Filter to keep only expressed circRNAs
most_sponging_circRNAs <- filter(most_sponging_circRNAs, ((H9_JReads + FB_JReads)/2) > 1)
most_sponging_circRNAs <- within(most_sponging_circRNAs, quartile <- as.integer(cut(Total_Sites_Per_circ, quantile(Total_Sites_Per_circ, type = 8), include.lowest=TRUE)))
most_sponging_circRNAs$quartile <- gsub(pattern = "1", replacement = "- sites", most_sponging_circRNAs$quartile)
most_sponging_circRNAs$quartile <- gsub(pattern = "2", replacement = "+ sites", most_sponging_circRNAs$quartile)
most_sponging_circRNAs$quartile <- gsub(pattern = "3", replacement = "++ sites", most_sponging_circRNAs$quartile)
most_sponging_circRNAs$quartile <- gsub(pattern = "4", replacement = "+++ sites", most_sponging_circRNAs$quartile)

# Also see the size of the circRNAs that interact with the different quartiles of miRNAs
miRNA_circRNA_Overlaps <- left_join(miRNA_circRNA_Overlaps, miRNA_sponging_coefficients[,c(1,7)])
for (i in 1:nrow(miRNA_circRNA_Overlaps)) {
  miRNA_circRNA_Overlaps$circ.size[i] <- sum(as.numeric(unlist(strsplit(as.character(miRNA_circRNA_Overlaps$Exon.sizes[i]), ','))))
}
miRNA_circRNA_Overlaps$quartile <- factor(miRNA_circRNA_Overlaps$quartile, level = c("- sponged", "+ sponged", "++ sponged", "+++ sponged"))
miRNA_circRNA_Overlaps <- filter(miRNA_circRNA_Overlaps, H9_JReads > 1)

# ------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------

# See what happens to the different quartiles of circRNAs according to the sites they have for miRNAs
      # Test for normality first, use KS because the number is too big for Shapiro
ks.test(x = most_sponging_circRNAs$log2FoldChange_circ, y = 'pnorm', alternative = 'two.sided') # If you run it, you'll see that the distribution is not normal.
     # Do a general linear model, not assuming normality. 
glmodel <- glm(data = most_sponging_circRNAs,
               formula = log2FoldChange_circ ~ quartile)
         # Follow up with estimated marginal means -emmeans- and contrasts.
stat.test <- emmeans(glmodel, pairwise ~ quartile)
stat.c <- cbind(data.frame(group1 = c("- sites", "- sites", "- sites", "+ sites", "+ sites", "++ sites"),
                           group2 = c("+ sites", "++ sites", "+++ sites", "++ sites", "+++ sites", "+++ sites")), as_tibble(stat.test$contrasts)) %>% add_significance()
stat.c <- stat.c[1:3,]

gmscb <- ggplot(data = most_sponging_circRNAs, aes(x = quartile, y = log2FoldChange_circ)) +
  geom_boxplot(outlier.shape = NA, aes(fill = quartile)) +
  geom_jitter(size = 0.1, alpha = 0.1) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, guide = "none", option = "D", direction = -1) +
  xlab("") + ylab("circRNA Fold change (log2)") + 
  ylim(-3,7) +
  stat_pvalue_manual(stat.c,
                     label = "p.value.signif", 
                     step_increase = 5, 
                     y.position = c(4,4.8,5.6),
                     tip.length = 0) +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(c(50,70,22,10)))
gmscb

# See if the size of the circRNAs varies between those that interact with the most sponged miRNAs and others
stat.test.s <- miRNA_circRNA_Overlaps %>%
  compare_means(formula = circ.size ~ quartile,
                method = "t.test",
                ref.group = "- sponged",
                p.adjust.method = "hochberg") %>% add_significance("p.adj")
stat.test.s
  
    # Plot
gcs <- ggboxplot(miRNA_circRNA_Overlaps,
            x = "quartile",
            y = "circ.size",
            fill = "quartile",
            width = 0.5,
            outlier.shape = "",
            add = c("violin"),
            add.params = list(alpha = 0.1, size = 0.5)) +  
  yscale(.scale = "log2") +
  expand_limits(y = c(0, 1300000)) +
  ylab("circRNA size -bp- (log2)") + xlab("") +
  theme_pubr(base_family = "Myriad Pro", base_size = 18, 
             legend = "none") +
  stat_pvalue_manual(stat.test.s,
                     label = "p.adj.signif", 
                     step_increase = 5, 
                     y.position = c(15,16,17),
                     tip.length = 0) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, guide = "none", option = "D", direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(c(50,70,22,10))) 
gcs 

# ------------------------------------------------------------------------
# Export the results
# ------------------------------------------------------------------------

write_csv(most_sponging_circRNAs, "./miRNA sponging coefficient/Results/circRNAs effective sites_2.0d.csv")


