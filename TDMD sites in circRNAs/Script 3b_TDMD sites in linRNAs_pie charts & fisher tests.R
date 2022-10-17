### Script 3b: Fisher tests and pie charts to assess enrichment of TDMD-like sites on linRNA through the quartiles of sponging coefficient
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("tidyverse")

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


if (!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("rkabacoff/ggpie")
library(ggpie)

# ------------------------------------------------------------------------
# Set the working directory and load the necessary tables. 
# ------------------------------------------------------------------------

# Set working directory
setwd("~/Influence of target RNA topology on microRNA stability/Bioinformatics/")

# Load tables (you need to have exported this table with script 2b)
sites_per_miR <- read_csv("./TDMD sites in circRNAs/Results/TDMD sites on lins per miR.csv")

# ------------------------------------------------------------------------
# Analysis 
# ------------------------------------------------------------------------

# Presence or absence of sites
sites_per_miR$group <- ifelse(sites_per_miR$sites == 0, "No TDMD-like sites found", "≥1 TDMD-like site")

# Count the number of miRNAs that have or not TDMD-like sites
total_sites <- sites_per_miR %>% group_by(quartile, group) %>% summarise(number = n())

# ------------------------------------------------------------------------
# Statistical analysis
# ------------------------------------------------------------------------

# Fisher test for quartile 1 vs 2: - sponged vs + sponged
fq1vs2 <- fisher.test(as.matrix(cbind(subset(total_sites, quartile == "- sponged")$number, 
                                      subset(total_sites, quartile == "+ sponged")$number)),
                      simulate.p.value = TRUE)
fq1vs2
# Fisher test for quartile 1 vs 3: - sponged vs ++ sponged
fq1vs3 <- fisher.test(as.matrix(cbind(subset(total_sites, quartile == "- sponged")$number, 
                                      subset(total_sites, quartile == "++ sponged")$number)),
                      simulate.p.value = TRUE)
fq1vs3
# Fisher test for quartile 1 vs 4: - sponged vs +++ sponged
fq1vs4 <- fisher.test(as.matrix(cbind(subset(total_sites, quartile == "- sponged")$number, 
                                      subset(total_sites, quartile == "+++ sponged")$number)),
                      simulate.p.value = TRUE)
fq1vs4

# Put the pvalues together in one data frame for plotting
pvalue_text <- data.frame(
                          label = c(paste("Fisher's Test, p = ", round(fq1vs2$p.value, 3)), 
                                    paste("Fisher's Test, p = ", round(fq1vs3$p.value, 3)), 
                                    paste("Fisher's Test, p = ", round(fq1vs4$p.value, 3))),
                          quartile = c("+ sponged", "++ sponged", "+++ sponged")
                          )

# ------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------

# Pie chart for
gpie <- ggpie::ggpie(data = sites_per_miR,
      x = group,
      by = quartile,
      nrow = 1,                    # number of rows
      border.color = "white",      # border color
      border.width = 1.5,          # border width
      label.color = "grey15",       # label color 
      label.size = 3,              # label size
      ) +
  scale_fill_manual(values = c("#fde72595", "#44015495")) + 
  geom_text(data = pvalue_text, 
            aes(label = label), 
            x = Inf, 
            y = -Inf,
            vjust = 20.5,
            hjust = 0.38,
            size = 2.5, 
            color = "grey30",
            inherit.aes = FALSE) +
  theme(plot.title.position = "plot",
        text = element_text(family = "Myriad Pro"),
        strip.text = element_text(vjust = 1, ),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        plot.margin = unit(c(1,5,5,5), 'mm')) 
gpie

