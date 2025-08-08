######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("rGREAT")
library("dplyr")
library("purrr")
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


#############################################################################################
### FUNCTION: processEnrichment() ---
processEnrichment <- function(dat){
    dat <- subset(dat, select=c("id","genome_fraction","observed_region_hits","observed_gene_hits","fold_enrichment","mean_tss_dist","gene_set_size","p_value","p_adjust"))
    dat <- dat[order(dat$p_value, decreasing=FALSE),]
    return(dat)
}

#############################################################################################
### LOAD ATAC-SEQ DIFFPEAKS ---
list.rds_atacseq_diffpeaks <- readRDS(file=file.rds_atacseq_diffpeaks)

### LOAD DIFFPEAKS: sgPTGES3 vs sgGAL4 ---
dat_diffPeaks <- list.rds_atacseq_diffpeaks$diffPeaks
df_diffpeaks <- list.rds_atacseq_diffpeaks$diffPeaks_pass

### ADD GEMOME FEATURE COORDINATES ---
df_diffpeaks$chr <- unlist(lapply(stringr::str_split(df_diffpeaks$FeatureID, "_"), function(x) x[1]))
df_diffpeaks$start <- as.numeric(unlist(lapply(stringr::str_split(df_diffpeaks$FeatureID, "_"), function(x) x[2])))
df_diffpeaks$end <- as.numeric(unlist(lapply(stringr::str_split(df_diffpeaks$FeatureID, "_"), function(x) x[3])))

### GET GRANGES OBJECT ---
gr_diffpeaks <- GenomicRanges::makeGRangesFromDataFrame(df_diffpeaks, keep.extra.columns=TRUE)


#############################################################################################
### SET SEED ---
set.seed(12345)

# GREAT ENRICHMENT: HALLMARKS ---
obj <- rGREAT::great(gr=gr_diffpeaks, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                            biomart_dataset = NULL, min_gene_set_size = 5, 
                            mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                            extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)

#dat <- processEnrichment(dat=obj@table)

# GET AR HALLMARKS ---
genelist_ar <- obj@gene_sets$HALLMARK_ANDROGEN_RESPONSE

# GET ANNOTATED REGIONS ---
gr_res <- getRegionGeneAssociations(obj)
names(gr_res) <- NULL
list.annotated_genes <- gr_res$annotated_genes

# GET THE INDEXES OF THE GENES THAT ARE IN THE AR HALLMARKS ---
index_arhallmarks <- which( unlist(lapply(list.annotated_genes, function(x) length( which( names(x) %in% genelist_ar ) ) )) != 0)
#list.annotated_genes[index_arhallmarks]

# EXTRACT PEAKS THAT ARE IN THE AR HALLMARKS ---
gr_res_arhallmarks <- gr_res[index_arhallmarks]


### SAVE OBJECT TO RDATA FILE ---
file.rds_res_arhallmark <- file.path(dir.reproduce_data, "sgPTGES3_diffpeaks_res_arhallmark.rds")
saveRDS(object=gr_res_arhallmarks, file=file.rds_res_arhallmark)


#############################################################################################
### LOAD ATAC-SEQ DIFFPEAKS: AR HALLMARKS ---
file.rds_res_arhallmark <- file.path(dir.reproduce_data, "sgPTGES3_diffpeaks_res_arhallmark.rds")
gr_res_arhallmark <- readRDS(file=file.rds_res_arhallmark)


# PARSE THE DIST. TO TSS ---
gr_res_arhallmark$dist_to_TSS <- purrr::map2(gr_res_arhallmark$annotated_genes, gr_res_arhallmark$dist_to_TSS, function(x, y){
                                                                                    names(y) <- names(x)
                                                                                    return(y)
                                                                                } )

gr_res_arhallmark$annotated_genes <- unlist(lapply(gr_res_arhallmark$annotated_genes, function(x) paste(x[names(x) %in% genelist_ar],collapse=":") ))
gr_res_arhallmark$dist_to_TSS <- unlist(lapply(gr_res_arhallmark$dist_to_TSS, function(x) paste(x[names(x) %in% genelist_ar],collapse=":") ))


###################
gr_res_arhallmark_promoter <- gr_res_arhallmark[c(1,4,5,14)]
gr_res_arhallmark_promoter$FeatureID





#############################################################################################
### LOAD ATAC-SEQ FEATURE COUNTS DATA ---
list.rds_atacseq_feature_counts <- readRDS(file=file.rds_atacseq_feature_counts)

### GET PEAK ANNOTATION DATA ---
dat_feature_annotation <- list.rds_atacseq_feature_counts$feature_annotation
dat_feature_annotation <- dat_feature_annotation[match(dat_diffPeaks$FeatureID, dat_feature_annotation$FeatureID),]

### MERGE DIFF PEAKS AND PEAKS ANNOTATION DATA ---
dat <- merge(dat_diffPeaks, dat_feature_annotation, by="FeatureID")

### GET FEATURE NORMALIZED COUNTS DATA ---
mat <- list.rds_atacseq_feature_counts$feature_counts_norm
mat <- subset(mat, select=c("sgGAL4-1","sgGAL4-2","sgPTGES3-1","sgPTGES3-2"))
colnames(mat) <- c("sgCTRL-1","sgCTRL-2","sgPTGES3-1","sgPTGES3-2")
#mat <- subset(mat, rownames(mat) %in% dat$FeatureID)


#############################################################################################
### PROMOTER REGIONS ---
mat_promoter <- subset(mat, rownames(mat) %in% gr_res_arhallmark_promoter$FeatureID)
mat_promoter <- mat_promoter[match(gr_res_arhallmark_promoter$FeatureID, rownames(mat_promoter)),]
rownames(mat_promoter) <- gr_res_arhallmark_promoter$annotated_genes


### PLOT ---
# GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_plots, "heatmap_atacseq_sgPTGES3_sgCTRL_AR_hallmarks_promoter.pdf")
p <- pheatmap(as.matrix(mat_promoter), 
				color = rev(jColFun(100)), 
				kmeans_k = NA, breaks = NA, border_color = "#000000",
				cellwidth = 18, cellheight = 10, 
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
				drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "",
				fontsize = 3, fontsize_row = 5, fontsize_col = 5,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 2, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 3, height = 2,
				silent = FALSE, na_col = "#DDDDDD")




#############################################################################################
### ALL REGIONS ---
mat_ar <- subset(mat, rownames(mat) %in% gr_res_arhallmark$FeatureID)
mat_ar <- mat_ar[match(gr_res_arhallmark$FeatureID, rownames(mat_ar)),]
rownames(mat_ar) <- gr_res_arhallmark$annotated_genes


### PLOT ---
# GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_plots, "heatmap_atacseq_sgPTGES3_sgCTRL_AR_hallmarks.pdf")
p <- pheatmap(as.matrix(mat_ar), 
				color = rev(jColFun(100)), 
				kmeans_k = NA, breaks = NA, border_color = "#000000",
				cellwidth = 18, cellheight = 5, 
            	scale = "row", 
				cluster_rows = TRUE, cluster_cols = FALSE, 
				clustering_distance_rows = "euclidean", 
				clustering_distance_cols = "euclidean", 
				clustering_method = "ward.D2",
				cutree_rows = NA, cutree_cols = NA,
				treeheight_row = 15, treeheight_col = 5, 
				legend = TRUE, legend_breaks = NA, legend_labels = NA,
				annotation_row = NA, annotation_col = NA, 
				annotation_colors = FALSE, annotation_legend = FALSE,
				annotation_names_row = FALSE, annotation_names_col = FALSE,
				drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "",
				fontsize = 3, fontsize_row = 5, fontsize_col = 5,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 2, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 3, height = 3,
				silent = FALSE, na_col = "#DDDDDD")


