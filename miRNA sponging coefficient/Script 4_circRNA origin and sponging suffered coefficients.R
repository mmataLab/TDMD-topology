### Script 4: Characteristics of the circRNAs
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse", "viridis", "fuzzyjoin", "ggpubr", "ggrepel")

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

# Load the GTF file with circRNA information from UCSC
GTF <- read_tsv("./miRNA sponging coefficient/hg19.refGene.gtf.gz", col_names=c("chromosome", "class", "geneFeature", "start", "end", "X6", "strand"))
# Load the interactions table produced by Script 1
miRNA_circRNA_Overlaps <- read_csv("./miRNA sponging coefficient/Results/miR_Cir_interactions_fulltable.csv")
# Load the table with the sponging coefficients
most_sponged_miRNAs <- read_csv("./miRNA sponging coefficient/Results/miRNA sponging coefficient.csv")

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Filter the GTF file
GTF <- filter(GTF, geneFeature != "transcript")
GTF <- filter(GTF, geneFeature != "exon")

# Add the sponging coefficients to the interactions table
miRNA_circRNA_Overlaps_sp <- left_join(miRNA_circRNA_Overlaps, most_sponged_miRNAs, by = "miRNAname")

# Join the interactions table with the GTF, with fuzzyjoin
miRNA_circRNA_Overlaps_sp <- genome_inner_join(miRNA_circRNA_Overlaps_sp, GTF, by = c('chromosome', 'start', 'end'))

# Build a table with the desired information
table_geneFeat_sp <- table(miRNA_circRNA_Overlaps_sp$geneFeature, miRNA_circRNA_Overlaps_sp$quartile) %>% 
  as.data.frame()

# Change the format
table_geneFeat_sp_w1 <- pivot_wider(table_geneFeat_sp, names_from = "Var2", values_from = "Freq") 
table_geneFeat_sp_w1 <- table_geneFeat_sp_w1 %>% remove_rownames %>% column_to_rownames(var="Var1")

# Test for enrichment before normalizing
chisq.test(table_geneFeat_sp_w1)

###Normalize to percentage
#Normalization function. 
norm <- function(x) {
  100*x / (sum(x))
}

#normalization
table_geneFeat_sp <- table_geneFeat_sp %>%
  group_by(Var2) %>%
  mutate(Percentage = norm(Freq))

# Test for enrichment after normalizing
table_geneFeat_sp$Freq <- NULL
table_geneFeat_sp_w2 <- pivot_wider(table_geneFeat_sp, names_from = "Var2", values_from = "Percentage") 
table_geneFeat_sp_w2 <- table_geneFeat_sp_w2 %>% remove_rownames %>% column_to_rownames(var="Var1")
chisq.test(table_geneFeat_sp_w2)

# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

ggbarplot(table_geneFeat_sp, x = "Var2", y = "Percentage", fill = "Var1", xlab = "Sponging suffered") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6, option = "D", direction = -1) + 
  theme_pubr(base_family = "Myriad Pro", base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(c(50,70,22,10))) +
  theme(legend.title=element_blank())


