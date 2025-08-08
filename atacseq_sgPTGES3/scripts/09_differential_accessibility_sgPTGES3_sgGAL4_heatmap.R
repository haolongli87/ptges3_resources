######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### LOAD LIBRARIES ---
library("stringr")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)


### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_masterdata.rds")
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_feature_counts.rds")
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks.rds")



###########################################################################################################################################################
### LOAD ATAC-SEQ DIFF PEAKS DATA ---
list.rds_atacseq_diffpeaks <- readRDS(file=file.rds_atacseq_diffpeaks)

### GET DIFF PEAKS ---
dat_diffPeaks <- list.rds_atacseq_diffpeaks$diffPeaks
dat_diffPeaks_pass <- list.rds_atacseq_diffpeaks$diffPeaks_pass



### LOAD ATAC-SEQ FEATURE COUNTS DATA ---
list.rds_atacseq_feature_counts <- readRDS(file=file.rds_atacseq_feature_counts)


### GET PEAK ANNOTATION DATA ---
dat_feature_annotation <- list.rds_atacseq_feature_counts$feature_annotation
dat_feature_annotation <- dat_feature_annotation[match(dat_diffPeaks$FeatureID, dat_feature_annotation$FeatureID),]


### MERGE DIFF PEAKS AND PEAKS ANNOTATION DATA ---
dat <- merge(dat_diffPeaks, dat_feature_annotation, by="FeatureID")
dat <- subset(dat, dat$FeatureID %in% dat_diffPeaks_pass$FeatureID)


### GET FEATURE NORMALIZED COUNTS DATA ---
mat <- list.rds_atacseq_feature_counts$feature_counts_norm
mat <- subset(mat, select=c("sgGAL4-1","sgGAL4-2","sgPTGES3-1","sgPTGES3-2"))
colnames(mat) <- c("sgCTRL-1","sgCTRL-2","sgPTGES3-1","sgPTGES3-2")
mat <- subset(mat, rownames(mat) %in% dat$FeatureID)


###########################################################################################################################################################
### PLOT ---
# GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_plots, "heatmap_fig_3g.pdf")
p <- pheatmap(as.matrix(mat), 
				color = rev(jColFun(100)), 
				kmeans_k = NA, breaks = NA, border_color = "#000000",
				cellwidth = 12, cellheight = NA, 
            	scale = "row", 
				cluster_rows = TRUE, cluster_cols = FALSE, 
				clustering_distance_rows = "euclidean", 
				clustering_distance_cols = "euclidean", 
				clustering_method = "ward.D2",
				cutree_rows = NA, cutree_cols = NA,
				treeheight_row = 5, treeheight_col = 5, 
				legend = TRUE, legend_breaks = NA, legend_labels = NA,
				annotation_row = NA, annotation_col = NA, 
				annotation_colors = FALSE, annotation_legend = FALSE,
				annotation_names_row = FALSE, annotation_names_col = FALSE,
				drop_levels = TRUE, show_rownames = FALSE, show_colnames = TRUE, main = "",
				fontsize = 3, fontsize_row = 3, fontsize_col = 3,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 2, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 2, height = 2,
				silent = FALSE, na_col = "#DDDDDD")

