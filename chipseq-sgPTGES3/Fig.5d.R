# Load the installed packages
library(ggplot2)
library(ggupset)
library(ggimage)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)

organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

ChIPseekerpath <- "~/HOMER/"

# Create a named list with file paths
files <- list(
  sgNTC = paste0(ChIPseekerpath, "AR_LNCaP_sgNTC.filtered.merged_final.bed"),
  sgPTGES3 = paste0(ChIPseekerpath, "AR_LNCaP_sgPTGES3.filtered.merged_final.bed")
)
print(files)

# write pdf
pdf("AR_LNCaP_sgNTC_sgPTGES3_2files.pdf", width = 8, height = 3)

# Annotate peaks
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 100), verbose=FALSE)

# Plot annotation results
plotAnnoBar(peakAnnoList)

# Close the PDF
dev.off()

