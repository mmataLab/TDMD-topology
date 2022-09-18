### Script 1: use scanMiR to find TDMD-like sites on circRNAs
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("scanMiRData", "scanMiR", "Biostrings", "tidyverse")

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
  # Get the desired circRNA seq
seqs <- readDNAStringSet("./TDMD site plots/circCSNK1G3_hsa.fasta")

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# Filter the list provided by scanMiR to keep only the desired miRNAs
mods <- getKdModels("hsa",NULL)
miRNA <- mods["hsa-miR-181d-5p"]

# Run the findSeedMatches function from scanMiR, which will provide a table of predicted sites on all the circRNA seqs provided by circBase 
## Careful here, this step demands a lot of computer power and free space if you want it done in a fair amount of time... 
m <- findSeedMatches(seqs, 
                     miRNA,
                     ret = "GRanges")

# See the table
m_df <- as.data.frame(m)

# Add the range of the site to identify an individual site when there is more than one possible
range <- "129-136"

# ------------------------------------------------------------------------
# Plot alingment
# ------------------------------------------------------------------------

viewTargetAlignment(
  m[ranges(m) == range],
  miRNA$`hsa-miR-181d-5p`$mirseq,
  seqs = seqs,
  UGsub = TRUE,
  flagBulgeMatches = TRUE,
  outputType = c("print")
)




