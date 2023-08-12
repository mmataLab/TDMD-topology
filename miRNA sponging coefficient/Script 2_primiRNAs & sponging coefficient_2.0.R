### Script 2: pri-miRNA analysis 
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
setwd("C:\\Users\\mmata\\Dropbox\\Mis manuscritos\\2022_circRNA paper\\R\\Bioinformatics")

# Load tables
primiR_expression <- read_delim("miRNA sponging coefficient/primiRNA_expression_GSE56152_GSE63709.txt", col_names = TRUE, delim = "\t") #This one includes pri-miRs
miRNA_sponging_coefficients <- read_csv("miRNA sponging coefficient/Results/miRNA sponging coefficient_2.0.csv")

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Only thing needed is to join the sponging coefficient data produced by script 1 and the pri-miR expression data
df <- left_join(x = miRNA_sponging_coefficients, y = primiR_expression[,c(1,4,6,7,8)], by = c("miRNAname" = "ID"))

# Multiple regression model with interaction 
regr_res <- lm(log2FoldChange_miR ~ log2FoldChange_primiR + quartile + log2FoldChange_primiR*quartile, data = df)
summary(regr_res)

# Calculate slopes from the linear model
c <- regr_res$coefficients
b1 <- as.numeric(c[2])
b2 <- b1 + as.numeric(c[6])
b3 <- b1 + as.numeric(c[7])
b4 <- b1 + as.numeric(c[8])

# Extract p values from the linear model
p <- summary(regr_res)$coefficients[, "Pr(>|t|)"]
p2 <- as.numeric(p[6])
p3 <- as.numeric(p[7])
p4 <- as.numeric(p[8])

# Put the pvalues together in one data frame for plotting
text <- data.frame(
  label = c(paste(" b =", round(b1, 2)), 
            paste("b =", round(b2, 2), ", p =", round(p2, 2)), 
            paste("b =", round(b3, 2), ", p =", round(p3, 2)), 
            paste("b =", round(b4, 2), ", p =", round(p4, 3))),
  quartile = c("- sponged", "+ sponged", "++ sponged", "+++ sponged")
)

# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

# Scatterplot of the pri-miR fold change vs the mature miRNA fold change, faceted by quartile

gpri <- ggplot(data = df, aes(x = log2FoldChange_primiR, y = log2FoldChange_miR)) +
  geom_point(alpha = 0.8, aes(color = quartile)) +
  stat_cor(method = "pearson", size = 4) +
  geom_smooth(se = TRUE, method = lm, color = "gray", alpha = 0.2) +
  scale_color_viridis(discrete = TRUE, alpha = 0.6, option = "D", guide = "none", direction = -1) +
  facet_wrap(~ quartile, ncol = 2) +
  xlab("pri-miR Fold change (log2)") + 
  ylab("miRNA Fold change (log2)") +
  theme_pubr(base_family = "Myriad Pro", base_size = 18) +
  theme(strip.background = element_rect(color = "black", fill = "white")) + 
  geom_text(data = text, 
            aes(label = label), 
            x = -Inf, 
            y = Inf,
            vjust = 18, 
            hjust = -0.1,
            size = 4, 
            color = "grey30",
            inherit.aes = FALSE) 
gpri  

