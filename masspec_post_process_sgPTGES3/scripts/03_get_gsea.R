######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### LOAD LIBRARIES ---
library("stringr")
library("data.table")
library("ggplot2")
library("gridExtra")
library("fgsea")
library("pheatmap")
library("RColorBrewer")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/masspec_post_process_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.dat <- file.path(dir.reproduce_data, "PTGES3ko_masspec_results_ttest.tsv")
file.gmt_hallmarks <- file.path(dir.reproduce_data, "h.all.v7.5.1.symbols.gmt") # DOWNLOADED FROM MSIGDB


#####################################################################################
### LOAD DESIGN TABLE ---
dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### GET A VECTOR OF FOLDCHANGE ---
vec_foldchange <- dat$FoldChange
names(vec_foldchange) <- dat$Protein
vec_foldchange <- sort(vec_foldchange, decreasing=TRUE)


#####################################################################################
### LOAD GENESETS ---
list.genesets_hallmarks <- fgsea::gmtPathways(gmt.file=file.gmt_hallmarks)

### RUN GSEA ---
res_gsea <- fgsea::fgsea(pathways=list.genesets_hallmarks, stats=vec_foldchange, minSize=2, maxSize=500, gseaParam=1)
res_gsea_pass <- subset(res_gsea, res_gsea$pval < 0.05)
res_gsea_pass <- res_gsea_pass[order(res_gsea_pass$ES, decreasing=TRUE),]

### SAVE OBJECT TO RDATA FILE ---
file.rds_masspec_gsea <- file.path(dir.reproduce_data, "sgPTGES3_masspec_gsea.rds")
saveRDS(object=res_gsea_pass, file=file.rds_masspec_gsea)



#####################################################################################
### ADD GROUP ---
res_gsea_pass$Group <- NA
res_gsea_pass$Group[which(res_gsea_pass$ES > 0)] <- "UP"
res_gsea_pass$Group[which(res_gsea_pass$ES < 0)] <- "DN"

### FACTORIZE ---
res_gsea_pass$Group <- factor(res_gsea_pass$Group, levels=c("UP","DN"))
res_gsea_pass$pathway <- factor(res_gsea_pass$pathway, levels=rev(res_gsea_pass$pathway))


# PLOT ---
p <- ggplot(res_gsea_pass, aes(x=ES, y=pathway)) +
            geom_bar(aes(fill=Group), stat="identity", color="#000000", width=0.8, size=0.2) +
            geom_vline(xintercept=0, color="#000000", linetype=1) +
            scale_fill_manual(values=c("#e31a1c","#1f78b4")) +
            coord_cartesian(xlim=c(-0.5,0.3)) +
            scale_x_continuous(breaks=seq(-0.5,0.3, by=0.1), labels=seq(-0.5,0.3, by=0.1)) +
		    theme(
			    axis.text.x = element_text(size = 5, color="#000000"),
			    axis.text.y = element_text(size = 5, color="#000000"),
			    axis.title = element_text(size = 5, color="#000000"),
			    plot.title = element_text(size = 5, color="#000000", hjust=0.5),
			    panel.grid.major.y = element_blank(),
                panel.grid.major.x = element_blank(),
			    panel.grid.minor = element_blank(),
			    axis.ticks = element_line(size = 0.2, color="#000000"),	
			    panel.background = element_rect(fill = "#FFFFFF", color = "#000000"),
			    legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none",
                legend.box = "horizontal") + 
            guides(fill = guide_legend(nrow = 2)) +                   
		    ylab("") +            
		    xlab("Enrichment Score") + 
            ggtitle("") 



# GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_plots, "PTGES3ko_masspec_results_ttest_gsea_barplot.pdf")
pdf(file.plot, height=3, width=4)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()






###############################################################################################################
# RE-LOAD MASSPEC GSEA DATA ---
file.rds_masspec_gsea <- file.path(dir.reproduce_data, "sgPTGES3_masspec_gsea.rds")
res_gsea_pass <- readRDS(file=file.rds_masspec_gsea)

# GET AR RESPONSE GENES ---
genes_ar_response <- res_gsea_pass$leadingEdge[which(res_gsea_pass$pathway == "HALLMARK_ANDROGEN_RESPONSE")][[1]]

### LOAD MASSPEC DATA ---
file.rds_masspec <- file.path(dir.reproduce_data, "sgPTGES3_masspec_expression.rds")
list.rds_masspec <- readRDS(file=file.rds_masspec)

expr <- list.rds_masspec$expr_norm
annot <- list.rds_masspec$annotation


### GET AR RESPONSE GENES EXPR ---
annot_ar_response <- subset(annot, annot$Protein %in% genes_ar_response)
expr_ar_response <- subset(expr, rownames(expr) %in% annot_ar_response$UniprotID)
rownames(expr_ar_response) <- annot_ar_response$Protein
colnames(expr_ar_response) <- stringr::str_replace_all(colnames(expr_ar_response), "sgGAL4","sgCTRL")



###############################################################################################################
### PLOT ---
# GENERATE COLOR PALETTE ---
jColFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))

# PLOT ---
file.plot <- file.path(dir.reproduce_plots, "heatmap_extended_data_fig_5b.pdf")
p <- pheatmap(as.matrix(expr_ar_response), 
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
            	filename = file.plot, width = 2, height = 2.5,
				silent = FALSE, na_col = "#DDDDDD")

