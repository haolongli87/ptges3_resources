######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3_inhibitor") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_masterdata.rds")
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_feature_counts.rds")
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_diffpeaks.rds")

file.msgidb_hallmark_pathway <- file.path(dir.reproduce, "atacseq_sgPTGES3/data/msigdb_h.all.v7.5.1.symbols.gmt")
file.background_genes <- file.path(dir.reproduce, "atacseq_sgPTGES3/data/genelist_protein_coding_genecodev28.txt")

###########################################################################################################################################################
### FUNCTION: gmtPathways() ---
#' Returns a list of pathways from a GMT file.
#' @param gmt.file Path to a GMT file.
#' @return A list of vectors with gene sets.
gmtPathways <- function(file.gmt){
    pathwayLines <- strsplit(readLines(file.gmt), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    return(pathways)
}

#### FUNCTION: get_over_representation_analysis() ---
get_over_representation_analysis <- function(list.genesets, genes.queryset, genes.refset, p.threshold){
    # CREATE RESULT DATAFRAME ---
	dat.result <- data.frame(Category=as.character(names(list.genesets)))

    # COMPUTE ONE-TAIL FISHER EXACT TEST ---
	for(i in 1:length(list.genesets)){
		genes.genesets <- sort(unique(unlist( list.genesets[[i]] )), decreasing=FALSE)
		genes.interest <- intersect(genes.queryset, genes.genesets)
  
		yy <- length(intersect(genes.genesets, genes.interest))
		yn <- length(intersect(genes.genesets, setdiff(genes.refset, genes.interest)))
		ny <- length(intersect(setdiff(genes.refset, genes.genesets), genes.interest))
		nn <- length(intersect(setdiff(genes.refset,genes.interest), setdiff(genes.refset, genes.genesets)))
  
		fisherRes <- fisher.test(rbind(c(yy,yn),c(ny,nn)), alternative="greater")
		dat.result$pvalue[i] <- fisherRes$p.value
		dat.result$fdr[i] <- NA

		dat.result$overlap.percent[i] <- round(length(genes.interest)/length(genes.genesets) * 100, digit=2)
		dat.result$overlap.genes[i] <- paste(genes.interest, collapse=":")
	}

	# MULTIPLE TEST CORRECTION ---
	dat.result$fdr <- p.adjust(dat.result$pvalue, method="BH")
	
    # TRIM RESULTS ---
	dat.result <- dat.result[order(dat.result$pvalue, decreasing=F),]
	dat.result <- subset(dat.result, dat.result$fdr <= p.threshold)
  
    # FORMAT PVALUES ---
	#dat.result$pvalue <- format(dat.result$pvalue, scientific = TRUE, digits = 4)
	#dat.result$fdr <- format(dat.result$fdr, scientific = TRUE, digits = 4)
	
    rownames(dat.result) <- NULL
	return(dat.result)
}  


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


### EXTRACT DATA ---
dat <- subset(dat, dat$GeneType == "protein_coding")
dat <- subset(dat, dat$Feature %in% c("Promoter","Distal Intergenic","Intron"))

### FILTER BY FDR AND FOLDCHANGE ---
dat_diffPeaks_pass_dn <- dat[which( (dat$padj <= 0.05) & (dat$log2FoldChange < -1)),]

### GET GENES LOSS ---
genes_loss <- unique(sort(dat_diffPeaks_pass_dn$Gene))





############################################################################################################################################
############################################ OVER REPRESNTATION ANALYSIS ###################################################################
### DEFINE GENES ---
list.genesets <- gmtPathways(file.gmt=file.msgidb_hallmark_pathway)
genes.queryset <- genes_loss
genes.refset <- data.table::fread(file=file.background_genes, sep="\t", header=FALSE, nThread=1,  data.table=FALSE, verbose=FALSE)$V1  #BACKGROUND GENESET


# GET OVER-REPRESENTATION ANALYSIS (HYPER-GEOMERTIC TEST) ---
dat_ora <- get_over_representation_analysis(list.genesets, genes.queryset, genes.refset, p.threshold=0.0001)






############################################################################################################################################
############################################ HEATMAP: AR-RESPONSE PATHWAY ##################################################################
genes_ar_response <- unlist(stringr::str_split( dat_ora$overlap.genes[which(dat_ora$Category == "HALLMARK_ANDROGEN_RESPONSE")], ":" ))
dat_diffPeaks_pass_dn_ar_response <- subset(dat_diffPeaks_pass_dn, dat_diffPeaks_pass_dn$Gene %in% genes_ar_response)

### LOAD ATAC-SEQ FEATURE COUNTS DATA ---
list.rds_atacseq_feature_counts <- readRDS(file=file.rds_atacseq_feature_counts)
dat_feature_counts <- list.rds_atacseq_feature_counts$feature_counts_norm

# SUBSET DATA BY FEATURES ---
mat <- subset(dat_feature_counts, rownames(dat_feature_counts) %in% dat_diffPeaks_pass_dn_ar_response$FeatureID)
mat <- mat[match(dat_diffPeaks_pass_dn_ar_response$FeatureID, rownames(mat)),]

dat_diffPeaks_pass_dn_ar_response$Feature_title <- stringr::str_replace_all( make.names( apply(dat_diffPeaks_pass_dn_ar_response, 1, function(x) paste( x[12], x[8], sep=" " )), unique=TRUE), "[.]", " ")
rownames(mat) <- dat_diffPeaks_pass_dn_ar_response$Feature_title



###############################################################################################################
### PLOT ---
# GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_plots, "heatmap_figure_5h.pdf")
p <- pheatmap(as.matrix(mat), 
				color = rev(jColFun(100)), 
				kmeans_k = NA, breaks = NA, border_color = NA,
				cellwidth = 10, cellheight = 6, 
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
				drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, main = "HALLMARK_ANDROGEN_RESPONSE",
				fontsize = 4, fontsize_row = 5, fontsize_col = 5,
				display_numbers = FALSE, number_format = "%.2f", number_color = "grey30", 
				fontsize_number = 2, 
            	gaps_row = NULL, gaps_col = NULL, 
				labels_row = NULL, labels_col = NULL, 
            	filename = file.plot, width = 3, height = 2.5,
				silent = FALSE, na_col = "#DDDDDD")

