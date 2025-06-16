######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### LOAD LIBRARIES ---
library("stringr")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)


### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_masterdata.rds")

file.consensuspeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_consensus.bed")
file.diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks.bed")
file.non_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_non_diffpeaks.bed")

file.matrix_consensuspeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_consensuspeaks_matrix.tsv.gz")
file.matrix_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks_matrix.tsv.gz")
file.matrix_non_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_non_diffpeaks_matrix.tsv.gz")

file.plot_diffpeaks <- file.path(dir.reproduce_plots, "heatmap_fig_4d_1.pdf")
file.plot_non_diffpeaks <- file.path(dir.reproduce_plots, "heatmap_fig_4d_2.pdf")

###########################################################################################################################################################
### FUNCTION: getComuteMatrix() ---
getComuteMatrix <- function(files.bigwig, file.bed, file.matrix, sampleids, n_cores=10){
    cmd <- paste("computeMatrix scale-regions",
                "-S", files.bigwig,
                "-R", file.bed,
                "-o", file.matrix,
                "-m", 500,
                "-b", 500,
                "-a", 500,
                "-bs", 5,
                "--samplesLabel", sampleids,
                "--verbose",
                "-p", n_cores,  
                sep=" ")
    
    return(cmd)
}

### FUNCTION: getPlotHeatmap() ---
getPlotHeatmap <- function(file.matrix, file.plot, plot_height=5, plot_width=2){
    cmd <- paste("plotHeatmap",
                "-m", file.matrix,
                "-o", file.plot,
                "--sortRegions", "descend",
                "--sortUsing", "mean",
                "--sortUsingSamples", "1 2",
                "--colorMap", "Reds",
                #"--whatToShow", "\"heatmap only\"",
                #"--whatToShow", "\"heatmap and colorbar\"",
                #"--whatToShow", "\"plot and heatmap\"",
                "--whatToShow", "\"plot, heatmap and colorbar\"",
                "--startLabel", "\"Start\"",
                "--endLabel", "\"End\"",
                "--legendLocation", "\"best\"",
                "--regionsLabel", "\"ATAC-seq Peaks\"",
                "--xAxisLabel", "\"Peaks (bp)\"", 
                "--heatmapHeight", plot_height,
                "--heatmapWidth", plot_width,
                #"--kmeans 2", 
                "--dpi", 600,
                "--verbose",
                sep=" ")
    
    return(cmd)
}


###########################################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### EXTRACT METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata
#metadata <- subset(metadata, metadata$CONDITION %in% c("sgGAL4","sgAR","sgPTGES3"))
metadata <- metadata[match(c("sgGAL4-1","sgGAL4-2","sgPTGES3-1","sgPTGES3-2","sgAR-1","sgAR-2"), metadata$SAMPLE_ID),]

### GET BIGWIG FILES ---
files.bigwig <- sapply(metadata$SAMPLE_ID, function(x) { paste(dir.wrk, "/processed/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.shifted.bw", sep="") })


###########################################################################################################################################################
### COMPUTE MATRIX -----
cmd1 <- getComuteMatrix(files.bigwig=paste(files.bigwig, collapse=" "), 
                        file.bed=file.diffpeaks, 
                        file.matrix=file.matrix_diffpeaks, 
                        sampleids=paste(metadata$SAMPLE_ID, collapse=" "), 
                        n_cores=50)
system(cmd1)

### PLOT HEATMAP -----
cmd2 <- getPlotHeatmap(file.matrix=file.matrix_diffpeaks, file.plot=file.plot_diffpeaks, plot_height=1, plot_width=2)
system(cmd2)


#### PLOT HEATMAP -----
#cmd2 <- getPlotHeatmap(file.matrix=file.matrix_diffpeaks, file.plot=file.plot_diffpeaks, plot_height=5, plot_width=2)
#system(cmd2)








###########################################################################################################################################################
### COMPUTE MATRIX -----
cmd3 <- getComuteMatrix(files.bigwig=paste(files.bigwig, collapse=" "), 
                        file.bed=file.non_diffpeaks, 
                        file.matrix=file.matrix_non_diffpeaks, 
                        sampleids=paste(metadata$SAMPLE_ID, collapse=" "), 
                        n_cores=50)
system(cmd3)

### PLOT HEATMAP -----
cmd4 <- getPlotHeatmap(file.matrix=file.matrix_non_diffpeaks, file.plot=file.plot_non_diffpeaks, plot_height=8, plot_width=2)
system(cmd4)
