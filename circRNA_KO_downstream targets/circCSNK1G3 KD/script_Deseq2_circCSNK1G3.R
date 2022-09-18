### Script DESeq2 to analyze circCSNK1g3 KD's impact on miR-181b/d-5p targets
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("DESeq2", "tidyverse", "viridis", "ggpubr", "ggrepel")

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

# Load tables
    # Load expression data
countdata <- read.table("./circRNA_KO_downstream targets/circCSNK1G3 KD/GSE113121_rnaseq_counts_circCSNK1G3.txt", sep = "\t", header = TRUE, row.names = 1)
    # Load the annotation table
annotations <- as.data.frame(read.csv("./circRNA_KO_downstream targets/circCSNK1G3 KD/Ensembl_Annotations.txt", sep = ",", header = TRUE))
    # Load miR-181-5p target list from TargetScan
miRTargets <- read.table("./circRNA_KO_downstream targets/circCSNK1G3 KD/TargetScan7.2__miR-181-5p.predicted_targets.txt", sep = "\t", header = TRUE)

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Display the experiment conditions in the appropiate way for DESeq2
conditions <- factor(c("circCSNK1G3KD","circCSNK1G3KD","circCSNK1G3OE","circCSNK1G3OEmock","circCSNK1G3OEmock","circCSNK1G3KDscbld","circCSNK1G3KDscbld"))
coldata <- data.frame(row.names = colnames(countdata), conditions)

# Create DeSeq dataset
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~conditions)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Results *in contrast: first put conditions, then the numerator, then the denominator. 
resT <- results(dds, cooksCutoff = FALSE, contrast = c("conditions","circCSNK1G3KD", "circCSNK1G3KDscbld"))

# Add the annotations to the results table
resT <- cbind(Gene.stable.ID.version = rownames(as.data.frame(resT)),as.data.frame(resT)) 
resT <- merge(resT, annotations, by.x = "Gene.stable.ID.version")

# Add a column indicating which genes are targeted by miR-181b/d-5p
    #First filter miR-181 b & d targets
keep <- c("hsa-miR-181d-5p", "hsa-miR-181b-5p")
miRTargets_Filtered<- subset(miRTargets, Representative.miRNA %in% keep)
    # Indicate which genes are targeted by miR-181b/d-5p
miRTargets_Filtered <- miRTargets_Filtered[,c(1,14,16)]
colnames(miRTargets_Filtered) <- c( "HGNC.symbol","miRNA", "Site_Score")
    # Join the table of targets with the results from DESeq
resT_miRTargets <- left_join(resT, miRTargets_Filtered, by = "HGNC.symbol")
    # Identify the genes that are miRNA targets
data <- mutate(resT_miRTargets, target = ifelse(is.na(miRNA), "notarget", "target")) 
    # Identify the genes that have a p-adj < 0.05
data$threshold = as.factor(data$padj < 0.05)
    # lose NAs
data <- data %>% drop_na(threshold)

# ------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------

#KDvsScrbld Volcano Plot

gKD <- ggplot() +
  geom_point(data = data, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.4, size = 0.8, color = "gray") +
  geom_point(data = subset(data, data$target == 'target'), 
             aes(x = log2FoldChange, y = -log10(padj), color = threshold, shape = threshold), 
             alpha = 0.8, size = 2) +
  scale_color_discrete(name = "p-Adj < 0.05", type = c("red", "blue")) +
  scale_shape_discrete(name = "p-Adj < 0.05") +
  xlab("Log2(Fold change)") + ylab("-Log10(p-adj)") + 
  geom_text_repel(data = subset(data, Gene.stable.ID.version == "ENSG00000100307.13"),
                  aes(x = log2FoldChange, y = -log10(padj)),
                  label = "CBX7",
                  box.padding = 5, size = 3, segment.size = 0.2 ) +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) +
  theme(legend.position = c(0.2,0.8), legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
gKD

# Box/violin plot of targets vs non targets from miR-181b/d-5p

gB <- ggboxplot(data = subset(data, data$threshold == TRUE),
               x = "target",
               y = "log2FoldChange",
               fill = "target",
               palette = "grey",
               add = "jitter",
               add.params = list(size = 0.1, alpha = 0.1),
               alpha = 0.3
) +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("notarget","target")), 
                     label = "p.signif", tip.length = 0, 
                     size = 6, label.y = 9) +
  ylim(-10,11) +
  xlab("") + 
  ylab("Fold change (log2)") + 
  scale_x_discrete(labels = c("Non\nTargets", "Targets")) + 
  theme_pubr(base_family = "Myriad Pro", base_size = 18, legend = "none")
gB
   

