### Script 1b: use scanMiR to find TDMD-like sites everywhere
### For inquiries: fuchsf@fbmc.fcen.uba.ar

# ------------------------------------------------------------------------
# You are going to need a few packages to do the analysis
# ------------------------------------------------------------------------

# First specify the packages of interest
packages = c("scanMiRData", "scanMiR", "Biostrings", "tidyverse", "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene", "scanMiRApp")

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
  # Extract transcript sequences
anno <- ScanMiRAnno(genome = BSgenome.Hsapiens.UCSC.hg19, ensdb = TxDb.Hsapiens.UCSC.hg19.knownGene)

seqs <- getTranscriptSequence(annotation = anno,extract = c("withORF"))

  # Get the list of microRNAs (you ought to have run Script 1 from the miRNA sponging coefficient folder)
miRs <- read_csv("./miRNA sponging coefficient.csv")

# ------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------

# First filter out or fix the miRNAs not annotated equally by scanMiR
  # Get all miRNAs from scanMiR
mods <- getKdModels("hsa",NULL)
  # Extract their names
miR_list <- names(mods)
  # See which ones do not match
miRs$in_scanMiR <- miRs$miRNAname %in% miR_list
  # Fix the annotation of the miRs (the seq to cross-check the strand was obtained from the Starbase table)
      # Turn the list into a vector
all <- as.vector(miRs$miRNAname)
      # Filter the vector
all_filtered <- all[all != "hsa-miR-137" & 
                    all != "hsa-miR-217" &
                    all != "hsa-miR-320a" &
                    all != "hsa-miR-375"]
      # Add the correct annotations
all_fixed <- c(all_filtered, "hsa-miR-137-3p", "hsa-miR-217-5p", "hsa-miR-320a-3p", "hsa-miR-375-3p")

# Create concatenated values with the list of the miRNAs
x <- dput(all_fixed)

# Filter the list provided by scanMiR to keep only the desired miRNAs
mods <- mods[x]

# Run the findSeedMatches function from scanMiR, which will provide a table of predicted sites on all the seqs provided
## Careful here, this step demands a lot of computer power and free space if you want it done in a fair amount of time... 
m <- findSeedMatches(seqs = seqs, 
                     seeds = mods,
                     onlyCanonical = FALSE,
                     BP = MulticoreParam(12, progressbar=TRUE),
                     maxLogKd=c(-1,-1.5),
                     ret=c("GRanges"),
                     useTmpFiles=TRUE)




# Filter to keep only TDMD-like sites
m_TDMD <- m[m$note == "TDMD?"]


# ------------------------------------------------------------------------
# Export the results
# ------------------------------------------------------------------------

saveRDS(m_TDMD,"./TDMD sites in circRNAs/Results/hg19_TDMD_spongedMiRs_UTR_ORF.rds")






