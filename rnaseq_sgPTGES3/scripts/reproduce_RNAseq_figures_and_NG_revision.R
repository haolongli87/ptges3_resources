# 20241219 Meng Zhang
# This script loads gene counts and perform DESeq2,
# and then generates figures for Haolong Li AR_PTGES3 RNAseq experiments
# Also runs gene set enrichment analysis 
set.seed(1)
##################################################################
library(stringr)
library(dplyr)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(gg.gap)
library(ggbreak)
library(ComplexHeatmap)
library(msigdbr)
library(xlsx)

######################################################
###### directories, filenames and variables ##########
######################################################
base_dir <- "C:/Users/mengz/Box\ Sync/manuscript/AR_PTGES3_Haolong_Li/RNAseq_figures_MZ/"
data_dir <- paste0(base_dir, "data_files/")
de_dir <- paste0(data_dir, "sgPTGES3_sgGAL4_LNCaP_DMSO/")
output_dir <- paste0(base_dir, "figures_and_NG_response_2024/")
hsp_output_dir <- paste0(output_dir, "hsp_analysis/")

##################
# Load data files
##################
fn_gene_count <- paste0(data_dir, "raw_gene_counts_matrix_gencode_v30.txt")
### The raw gene counts can be downloaded from our GEO deposit

fn_hallmark_pathways <- paste0(data_dir, "h.all.v6.2.symbols.gmt")
fn_ensemblinfo <- paste0(data_dir, "gencode.v30_gene_annotation_table.txt") # FILE DOWNLOADE FORM: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz

ensemblinfo <- read.delim(fn_ensemblinfo, header = TRUE, stringsAsFactors = FALSE)
rownames(ensemblinfo) <- ensemblinfo$Geneid
ensemblinfo$Geneid <- NULL
ensemblinfo$GeneSymbol <- as.character(lapply(ensemblinfo$GeneSymbol, function(x) gsub(" ", "", x, fixed = T)))
pathways.hallmark <- gmtPathways(fn_hallmark_pathways)

gene_count <- read.delim(fn_gene_count, check.names = F)
rownames(gene_count) <- gene_count$feature_id
gene_count$feature_id <- NULL

# load msigdb pathways
wiki_pathway = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") %>%  dplyr::select(gs_name, gene_symbol)
gomf_pathway = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%  dplyr::select(gs_name, gene_symbol)
reactome_pathway = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%  dplyr::select(gs_name, gene_symbol)
custom_pathway_df <- rbind(wiki_pathway, gomf_pathway, reactome_pathway)
custom_pathway = split(x = custom_pathway_df$gene_symbol, f = custom_pathway_df$gs_name)

##################
# Set threshold
##################
padj_threshold <- 0.05
log2fc_threshold <- 1
gsea_padj_threshold <- 0.01

##################
# DESeq2 analysis
##################
# make the sample table
sample_info <- data.frame(sample_id = colnames(gene_count))
sample_info$condition <- str_split_fixed(sample_info$sample_id, "-", 3)[,2]
rownames(sample_info) <- sample_info$sample_id
all(rownames(sample_info) == colnames(gene_count))

# construct dds for DE analysis from count matrix
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = sample_info, design = ~ condition)

# prefiltering of low count genes and 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- factor(dds$condition, levels = c("sgGAL4_DMSO", "sgPTGES3_DMSO"))

# DESeq2
dds <- DESeq(dds)
res <- results(dds, alpha = padj_threshold)
summary(res)

resOrdered <- res[order(res$padj),]
sum(res$padj < padj_threshold, na.rm = T)

# combine DE results, gene info, and normalized counts into one table
summary <- data.frame()

normalized_counts <- counts(dds, normalized=TRUE)
colnames(normalized_counts) <- lapply(colnames(normalized_counts), function(x) paste(x, "_normalized_counts", sep = ""))

summary <- transform(merge(as.data.frame(resOrdered), ensemblinfo, by = 0), row.names = Row.names, Row.names = NULL)
summary <- transform(merge(summary, normalized_counts, by = 0), row.names = Row.names, Row.names = NULL)

# order the list by adjusted pvalue
summary <- summary[order(summary$padj),]

# remove genes with NA padj values (these are genes with low expression levels)
summary_filtered <- summary[!is.na(summary$padj),]  

# rename columns of normalized counts
rownames(normalized_counts) <- ensemblinfo$GeneSymbol[match(rownames(normalized_counts), rownames(ensemblinfo))]
normalized_counts_log2 <- log2(normalized_counts + 1)
colnames(normalized_counts_log2) <- c("sgGAL4_rep1", "sgGAL4_rep2", "sgPTGES3_rep1", "sgPTGES3_rep2")
normalized_counts_mean <- data.frame(sgGAL4_mean = rowMeans(normalized_counts[, c(1,2)]),
                                     sgPTGES3_mean = rowMeans(normalized_counts[, c(3,4)]))
normalized_counts_mean$sgGAL4_mean_log2 <- log2(normalized_counts_mean$sgGAL4_mean + 1)
normalized_counts_mean$sgPTGES3_mean_log2 <- log2(normalized_counts_mean$sgPTGES3_mean + 1)
normalized_counts_mean$geneName <- rownames(normalized_counts)
################################################################################
# Volcano plot highlighting top 10 upregulated and downregulated genes + AR
################################################################################
fn_volcano <- paste0(output_dir, "volcano_plot_RNAseq_PTGES3_GAL4_LNCaP.pdf")
summary_filtered <- mutate(summary_filtered, sig=ifelse(summary_filtered$padj < padj_threshold &
                                                          abs(summary_filtered$log2FoldChange) > log2fc_threshold, sprintf("padj<%s, abs(FoldChange)>2", padj_threshold), "Not Sig")) 
summary_filtered_ranked <- summary_filtered[with(summary_filtered, order(stat)),]
label_df <- rbind.data.frame(head(summary_filtered_ranked, 10), tail(summary_filtered_ranked, 10), 
                             summary_filtered_ranked[which(summary_filtered_ranked$GeneSymbol == "AR"),])
volc = ggplot(summary_filtered, aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + #add points colored by significance
  scale_color_manual(values=c("gray", "red")) +
  xlim(-4, 4) +
  ylim(0, 120) +
  geom_text_repel(data=label_df,
                  aes(label=GeneSymbol), size = 4) +
  theme_bw() + theme(text = element_text(size = 15), legend.title = element_blank()) 

pdf(fn_volcano, width = 8, height = 5)
print(volc)
dev.off()

################################################################################
# GSEA showing AR signaling pathway is top downregulated upon LNCaP knockdown
################################################################################
fn_gsea_nes <- paste0(output_dir, "GSEA_RNAseq_PTGES3_GAL4_LNCaP_NES.pdf")

genes_rank <- summary_filtered$stat
names(genes_rank) <- summary_filtered$GeneSymbol
genes_rank <- genes_rank[names(genes_rank)[!duplicated(names(genes_rank))]]

fgseaRes <- fgsea(pathways = pathways.hallmark, stats = genes_rank, nperm = 100000, maxSize = 500, nproc = 6)

theDF <- fgseaRes

theDF$theLabel <- NA
theDF$theLabel[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)] <- theDF$pathway[theDF$padj < gsea_padj_threshold & !is.na(theDF$padj)]

# Make barplot
#theDF_sig <- theDF[order(abs(theDF$NES), decreasing = T)[1:nTop],]
# only plot significantly downregulated pathways
theDF_sig <- theDF[theDF$padj < gsea_padj_threshold &
                     theDF$NES < 0,]

p_sig_nes <- ggplot(theDF_sig, aes(reorder(pathway,-NES),NES)) +
  geom_bar(stat="identity", fill="#006094") +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title="") +
  theme(
    axis.text.x = element_text(size = 6, color="#000000"),
    axis.text.y = element_text(size = 6, color="#000000"),
    axis.title = element_text(size = 5, color="#000000"),
    plot.title = element_text(size = 12, color="#000000", hjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size=0.4, color="#000000"), 
    strip.text = element_text(size=10, color="#000000"),
    strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
    panel.background = element_rect(fill="#FFFFFF", color="#000000"),
    legend.text = element_text(size = 10, color="#000000"),
    legend.title = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.position = "none") 
print(p_sig_nes)

pdf(fn_gsea_nes, width = 3, height = 2)
print(p_sig_nes)
dev.off()

################################################################################
# NG responses: 
# GSEA enrichment plot of Hallmark AR response pathway
################################################################################
fn_enrichment <- paste0(output_dir, "GSEA_RNAseq_PTGES_GLA4_LNCaP_enrichment_HALLMARK_AR.pdf")

enrichment_plot <- plotEnrichment(pathway = pathways.hallmark[["HALLMARK_ANDROGEN_RESPONSE"]], 
               gseaParam = 1, ticksSize = 0.3, stats= genes_rank) + 
  labs(title="HALLMARK_ANDROGEN_RESPONSE") + theme(plot.title = element_text(hjust = 0.5, face="bold"))

pdf(fn_enrichment, width = 6, height = 4)
print(enrichment_plot)
dev.off()

################################################################################
# NG responses: 
# Heatmap of genes in the AR response pathways
################################################################################
fn_heatmap_all_genes_ar_pathway <- paste0(output_dir, "Heatmap_all_genes_HALLMARK_AR_sgPTGES3_sgGAL4.pdf")
fn_heatmap_de_genes_ar_pathway <- paste0(output_dir, "Heatmap_DE_genes_HALLMARK_AR_sgPTGES3_sgGAL4.pdf")

ar_genes <- pathways.hallmark[["HALLMARK_ANDROGEN_RESPONSE"]]

# rank the genes by log2fc
ar_genes_ranking <- summary_filtered[which(summary_filtered$GeneSymbol %in% ar_genes),]
ar_genes_ranking <- ar_genes_ranking[with(ar_genes_ranking, order(log2FoldChange)),]
ar_genes_ranking_de <- ar_genes_ranking[which(!ar_genes_ranking$sig == "Not Sig"),]

normalized_counts_log2_select <- normalized_counts_log2[ar_genes_ranking$GeneSymbol,]
normalized_counts_log2_select_de <- normalized_counts_log2[ar_genes_ranking_de$GeneSymbol,]

pdf(fn_heatmap_all_genes_ar_pathway, width = 6, height = 15)
row_ha = rowAnnotation(log2fc = anno_barplot(ar_genes_ranking$log2FoldChange))
Heatmap(normalized_counts_log2_select, name = "log2(normalized gene counts)", 
        cluster_columns = F, cluster_rows = F, show_row_dend = F, show_column_names = T, 
        right_annotation = row_ha, 
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(direction = "horizontal"))

dev.off()

pdf(fn_heatmap_de_genes_ar_pathway, width = 4.5, height = 3)
row_ha = rowAnnotation(log2fc = anno_barplot(ar_genes_ranking_de$log2FoldChange))
Heatmap(normalized_counts_log2_select_de, name = "log2(normalized gene counts)", 
        cluster_columns = F, cluster_rows = F, show_row_dend = F, show_column_names = T, 
        right_annotation = row_ha,
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(direction = "horizontal"))
dev.off()

##############################################
# NG responses: 
# Heatmap of DE genes in enriched pathways
##############################################
fn_heatmap_de_genes_GSEA_enriched_pathways <- paste0(output_dir, "Heatmap_DE_genes_GSEA_enriched_pathways.pdf")

theDF_sig_both_dir <- theDF[theDF$padj < gsea_padj_threshold,]
pathways_sig <- theDF_sig_both_dir$pathway
# create a pathway-gene map
pathway_gene_map_list <- list()
for (pathway in pathways_sig) {
  pathway_gene_map_list[[pathway]] <- data.frame(gene = pathways.hallmark[[pathway]],
                                                 pathway = pathway)
}
pathway_gene_map <- rbindlist(pathway_gene_map_list)
pathway_gene_map$NES <- theDF_sig_both_dir$NES[match(pathway_gene_map$pathway, theDF_sig_both_dir$pathway)]

# Save a table of all DE genes with annotation of enriched pathways from GSEA
fn_de_genes <- paste0(output_dir, "DE_genes_sgPTGES3_sgGAL4.csv")
de_genes <- summary_filtered_ranked[which(!summary_filtered_ranked$sig == "Not Sig"),] #257
de_genes$GESA_enriched_pathway <- pathway_gene_map$pathway[match(de_genes$GeneSymbol, pathway_gene_map$gene)]
de_genes$GSEA_NES <- pathway_gene_map$NES[match(de_genes$GeneSymbol, pathway_gene_map$gene)]
de_genes$GSEA_direction <- de_genes$GSEA_NES
de_genes$GSEA_direction[de_genes$GSEA_NES < 0] <- "Down"
de_genes$GSEA_direction[de_genes$GSEA_NES > 0] <- "Up"
write.csv(de_genes, fn_de_genes)

# rank the genes by log2fc
pathway_genes_ranking <- de_genes[which(!is.na(de_genes$GESA_enriched_pathway)),]
pathway_genes_ranking <- pathway_genes_ranking[with(pathway_genes_ranking, order(GSEA_NES, log2FoldChange)),]

normalized_counts_log2_select_pathway_genes <- normalized_counts_log2[pathway_genes_ranking$GeneSymbol,]

row_ha = rowAnnotation(log2fc = anno_barplot(pathway_genes_ranking$log2FoldChange))
pathway_ha <- rowAnnotation(GSEA_direction = anno_block(gp = gpar(fill = c(rep(3,9),2)),
                               labels = c(rep("Down", 9), "Up"), 
                               labels_gp = gpar(col = "white", fontsize = 10)))
lgd = Legend(labels = c("Down", "Up"), title = "GSEA_direction", legend_gp = gpar(fill = c(3,2)))
ht <- Heatmap(normalized_counts_log2_select_pathway_genes, name = "log2(normalized gene counts)", 
        cluster_columns = F, cluster_rows = F, show_row_dend = F, show_column_names = T, 
        row_split = factor(pathway_genes_ranking$GESA_enriched_pathway, levels = unique(pathway_genes_ranking$GESA_enriched_pathway)), 
        row_title_rot = 0,row_gap = unit(3, "mm"),border = TRUE,
        row_order = rownames(normalized_counts_log2_select_pathway_genes),
        left_annotation = pathway_ha,
        right_annotation = row_ha, 
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(direction = "horizontal"))

pdf(fn_heatmap_de_genes_GSEA_enriched_pathways, width = 10, height = 10)

ComplexHeatmap::draw(ht, annotation_legend_list = lgd)

dev.off()


################################################################################
# NG Response: 
# GSEA showing Heat shock protein pathways
################################################################################
fn_gsea_nes <- paste0(output_dir, "GSEA_RNAseq_PTGES3_GAL4_LNCaP_NES.pdf")

genes_rank <- summary_filtered$stat
names(genes_rank) <- summary_filtered$GeneSymbol
genes_rank <- genes_rank[names(genes_rank)[!duplicated(names(genes_rank))]]

fgseaRes <- fgsea(pathways = custom_pathway, stats = genes_rank, maxSize = 500, nproc = 6)
fgseaRes <- fgseaRes[order(fgseaRes$NES),]

selected_hsp_pathways <- c("REACTOME_HSF1_ACTIVATION", 
                           "REACTOME_HSF1_DEPENDENT_TRANSACTIVATION",
                           "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
                           "WP_FAS_LIGAND_PATHWAY_AND_STRESS_INDUCTION_OF_HEAT_SHOCK_PROTEINS",
                           "GOMF_HEAT_SHOCK_PROTEIN_BINDING",
                           "GOMF_HSP90_PROTEIN_BINDING", 
                           "GOMF_HSP70_PROTEIN_BINDING")
fgseaRes_hsp <- fgseaRes[which(fgseaRes$pathway %in% selected_hsp_pathways),]
fgseaRes_hsp$leadingEdge <- as.character(fgseaRes_hsp$leadingEdge)

fn_fgsea <- paste0(hsp_output_dir, "fgsea_hsp_PTGES3_GAL4.xlsx")
write.xlsx(fgseaRes_hsp, fn_fgsea)

# Per Luke request, save all GSEA result
fgseaRes_all_pathways <- fgsea(pathways = c(custom_pathway, pathways.hallmark), stats = genes_rank, maxSize = 500, nproc = 6)
fgseaRes_all_pathways <- fgseaRes_all_pathways[order(fgseaRes_all_pathways$NES),]
View(fgseaRes_all_pathways)
fgseaRes_all_pathways$leadingEdge <- as.character(fgseaRes_all_pathways$leadingEdge)
fn_fgsea_all_pathways <- paste0(hsp_output_dir, "fgsea_all_pathways_sgPTGES3_GAL4.xlsx")
write.xlsx(fgseaRes_all_pathways, fn_fgsea_all_pathways)

###################################################################
# NG responses: 
# Plot the mean counts for sgGAL4 and sgPTGES3 for the HSF1 genes
###################################################################
fn_cor_plot <- paste0(hsp_output_dir, "cor_sgGAL4_sgPTGES3_REACTOME_HSF1_ACTIVATION_genes.pdf")

colnames(normalized_counts_mean)
hsf_genes <- custom_pathway[["REACTOME_HSF1_ACTIVATION"]]
r_hsf_genes <- normalized_counts_mean$geneName %in% hsf_genes
min <- min(c(normalized_counts_mean$sgGAL4_mean_log2[r_hsf_genes], normalized_counts_mean$sgPTGES3_mean_log2[r_hsf_genes]))
max <- max(c(normalized_counts_mean$sgGAL4_mean_log2[r_hsf_genes], normalized_counts_mean$sgPTGES3_mean_log2[r_hsf_genes]))

plot <- ggplot(data = normalized_counts_mean[normalized_counts_mean$geneName %in% hsf_genes,], 
               aes(x = sgGAL4_mean_log2, y = sgPTGES3_mean_log2)) + geom_point() + 
  geom_abline(slope = 1,intercept = 0, color = "blue", linewidth = 1) +
  xlim(min, max) + ylim(min, max) +
  xlab("sgGAL4 log2 normalized gene counts") +
  ylab("sgPTGES3 log2 normalized gene counts") +
  theme_bw(base_size = 12, base_rect_size = 1.5) + theme(text = element_text(size = 15),
                                                         axis.text = element_text(face = "bold"),
                                                         axis.line.x = element_line(),
                                                         axis.line.y = element_line(),
                                                         panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_rect(fill = "transparent",colour = NA),
                                                         plot.background = element_rect(fill = "transparent",colour = NA))

print(plot)

pdf(fn_cor_plot, width = 6, height = 6)
print(plot)
dev.off()

