### Script DESeq2 to analyze CDR1as KO's impact on miR-7-5p targets
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
WT_1 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443732_CRTX_WT_rep1_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "WT_1"))
WT_2 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443733_CRTX_WT_rep2_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "WT_2"))
WT_3 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443734_CRTX_WT_rep3_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "WT_3"))
KO_1 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443735_CRTX_KO_rep1_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "KO_1"))
KO_2 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443736_CRTX_KO_rep2_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "KO_2"))
KO_3 <- read_delim("./circRNA_KO_downstream targets/CDR1as KO/GSM2443737_CRTX_KO_rep3_read_counts_ens67.txt", delim = "\t", col_names = c("ENSMBL_ID", "KO_3"))
  # Load the annotation table
annotations <- as.data.frame(read.delim("./circRNA_KO_downstream targets/CDR1as KO/Ensembl_to_Genenames.txt", sep = "\t", header = FALSE))
  # Load miR-7-5p target list from TargetScan
miRTargets <- read.table("./circRNA_KO_downstream targets/CDR1as KO/TargetScan7.1__miR-7-5p.predicted_targets_mouse.txt", sep = "\t", header = TRUE)

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Put together the expression table
countdata <- left_join(WT_1,WT_2) %>% left_join(.,WT_3) %>% left_join(.,KO_1) %>% left_join(.,KO_2) %>% left_join(.,KO_3)
countdata <- countdata %>% remove_rownames %>% column_to_rownames(var = "ENSMBL_ID")

# Display the experiment conditions in the appropriate way for DESeq2
conditions <- factor(c("WT","WT","WT","CDR1as KO","CDR1as KO","CDR1as KO"))
coldata <- data.frame(row.names = colnames(countdata), conditions)

# Create DeSeq dataset
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~conditions)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Results *in contrast: first put conditions, then the numerator, then the denominator. 
resT <- results(dds, cooksCutoff = FALSE, contrast = c("conditions", "CDR1as KO", "WT"))

# Add the annotations to the results table
resT <- cbind(V1 = rownames(as.data.frame(resT)),as.data.frame(resT)) 
resT <- merge(resT, annotations, by.x = "V1")

# Add a column indicating which genes are targeted by miR-7-5p
  # Indicate which genes are targeted by miR-7-5p
miRTargets_Filtered <- miRTargets[,c(1,14,16)]
colnames(miRTargets_Filtered) <- c("V2", "miRNA", "Site_Score")
  # Join the table of targets with the results from DESeq
resT_miRTargets <- left_join(resT, miRTargets_Filtered, by = "V2")
  # Identify the genes that are miRNA targets
data <- mutate(resT_miRTargets, target = ifelse(is.na(miRNA), "Non-Targets", "Targets")) 
  # Identify the genes that have a p-adj < 0.05
data$threshold = as.factor(data$padj < 0.05)
  # lose NAs
data <- data %>% drop_na(threshold)

# ------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------

# KO vs scrambled Volcano Plot

gKO <- ggplot() +
  geom_point(data = data, aes(x = log2FoldChange, y = -log10(padj)), alpha = 0.4, size = 0.8, color = "gray") +
  geom_point(data = subset(data, data$target == 'Targets'), 
             aes(x = log2FoldChange, y = -log10(padj), color = threshold, shape = threshold), 
             alpha = 0.8, size = 2) +
  scale_color_discrete(name = "p-Adj < 0.05", type = c("red", "blue")) +
  scale_shape_discrete(name = "p-Adj < 0.05") +
  xlab("Fold change (log2)") + ylab("-Log10(p-adj)") +
  geom_text_repel(data = subset(data, data$threshold == 'TRUE' & data$target == 'Targets'),
                  aes(x = log2FoldChange, y = -log10(padj)),
                  label = subset(data, data$threshold == 'TRUE' & data$target == 'Targets')$V2,
                  box.padding = 1, size = 3, segment.size = 0.3 ) +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) +
  theme(legend.position = c(0.2,0.8), legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
gKO

# Cumulative proportions
  # Do a Wilcoxon test to compare target vs non targets
WX.p.Val <- wilcox.test(log2FoldChange ~ target, data = data)$p.value
WX.p.Val <- WX.p.Val %>% format(digits=3)
WXP <- paste("WX p-value=", WX.p.Val)
  # Plot
gEDCF_Cdr1asKO <- ggplot(data, 
                         aes(x=log2FoldChange, 
                             colour = target)) + 
  stat_ecdf(size=1) + 
  theme_pubr(base_family = "Myriad Pro", base_size = 18) +
  theme(legend.title=element_blank()) +
  labs(x = "Fold change (log2)", y="Cumulative proportion") + 
  annotate("text", x = -0.5, y = 0.9, label = WXP, size = 4) + 
  scale_x_continuous(limits = c(-0.85,0.85), breaks=seq(-0.8, 0.8, 0.4)) + 
  geom_vline(xintercept = 0, linetype="dashed",size=0.5)

gEDCF_Cdr1asKO
   
