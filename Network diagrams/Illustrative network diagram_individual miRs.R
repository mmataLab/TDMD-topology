### Script network diagram: an illustration of the interaction of selected miRNAs with circRNAs
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse", "fuzzyjoin", "igraph")

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
# Analysis: a lot of this is just a repetition of Script 1 from miRNA sponging coefficient.
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

# Choose a particular miRNA and the circRNAs it interacts with, then build the edge list
df <- subset(miRNA_circRNA_Overlaps, miRNAname == "hsa-miR-7-5p" & H9_JReads > 5) # pick your miRNA here
df <- df[,c("miRNAname", "geneName")]
df$miRNAname <- gsub(df$miRNAname, pattern = "hsa-", replacement = "")

# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

# Create the network data
network <- graph_from_data_frame(df, directed = F)         

# Assign colors
V(network)$color <- adjustcolor("#fde725", alpha.f = .2)
V(network)[1]$color <- "#44015470"

# Plot using igraph
plot(network,
     vertex.size = 70,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     vertex.frame.color = "transparent",
     layout = layout.kamada.kawai(network)
)

