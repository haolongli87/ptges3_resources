######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### LOAD LIBRARIES ---
library("stringr")
library("ggplot2")
library("gridExtra")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_masterdata.rds")
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_feature_counts.rds")
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks.rds")

###########################################################################################################################################################
### FUNCTION: getVolcanoPlot() ---
getVolcanoPlot <- function(df, genes, logfc.threshold, pvalue.threshold, filterby.pvalueadj, fc.upperbound, nlogp.upperbound){
	yval <- -log10(pvalue.threshold)
    xval <- logfc.threshold
	cbPalette1 <- c("#98a3a5","#e85748")
    
    # REMOVE NAs ---
    if(length(which(is.na(df$pvalue))) > 0){
        df <- df[-which(is.na(df$pvalue)),]
    }

    # PVALUE FILTER ---
    if(filterby.pvalueadj == TRUE){
        df$nlogp <- -log10(df$padj)
    }else{
        df$nlogp <- -log10(df$pvalue)
    }

    # GROUP GENES ---
    df$Group <- "STABLE"
    if(filterby.pvalueadj == TRUE){    
        df$Group[which((df$log2FoldChange >= logfc.threshold) & (df$padj <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$log2FoldChange <= -logfc.threshold) & (df$padj <= pvalue.threshold))] <- "VARIABLE"
    }else{
        df$Group[which((df$log2FoldChange >= logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$log2FoldChange <= -logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
    }

    # ADJUST FOLDCHANGE FOR VIZ ---
    if(max(df$log2FoldChange) > fc.upperbound){
        df$log2FoldChange[which(df$log2FoldChange > fc.upperbound)] <- fc.upperbound 
    }

    if(min(df$log2FoldChange) < -fc.upperbound){
        df$log2FoldChange[which(df$log2FoldChange < -fc.upperbound)] <- -fc.upperbound
    }

    # ADJUST nlogp FOR VIZ ---
    if(max(df$nlogp) > nlogp.upperbound){
        df$nlogp[which(df$nlogp > nlogp.upperbound)] <- nlogp.upperbound 
    }

    # ADD GENE LABELS BY SELECTED GENES ---
    df$Label <- NA
    df$Label[which(df$Gene %in% genes)] <- df$Gene[which(df$Gene %in% genes)] 
    df$Label[which(df$Group == "STABLE")] <- NA

    # FACTORIZE ---
    df$Group <- factor(df$Group, levels=c("STABLE","VARIABLE"))

    # PLOT ---
	p <- ggplot(df, aes(x=log2FoldChange, y=nlogp, label=Label)) + 
			geom_hline(yintercept = yval, colour="gray70", linetype=4, linewidth=0.5, alpha=0.7) +
			geom_vline(xintercept = xval, colour="gray70", linetype=4, linewidth=0.5, alpha=0.7) +		
			geom_vline(xintercept = -xval, colour="gray70", linetype=4, linewidth=0.5, alpha=0.7) +
			geom_point(aes(fill=Group, color=Group), shape = 21, stroke=0.3, size=0.5, alpha=0.8) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette1) +
            scale_x_continuous(breaks=seq(-fc.upperbound, fc.upperbound, 1)) +
            scale_y_continuous(breaks=seq(0, 10, 2)) +
			coord_cartesian(xlim=c(-fc.upperbound, fc.upperbound), ylim=c(0, nlogp.upperbound)) +
			geom_text(aes(label=Label, alpha=0.3, hjust=1, vjust=-1, lineheight=1.5), color="black", size=2, na.rm =TRUE) +
			theme(
				axis.text = element_text(size = 6, color="black"),
				axis.title = element_text(size = 8, color="black"),
				strip.text = element_text(size = 6, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 8, color="black", hjust=0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(linewidth=0.3, color="black"),	
				panel.background = element_rect(fill="white", color="black"),
				legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.3, "cm"),			
			    legend.position="none") +
			xlab("log2 Fold Change") + 
			ylab("-log10(adj. pvalue)") + 
			ggtitle("Differentially Accessible Genomic Regions") 	

    return(p)
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



###########################################################################################################################################################
################# VOLCANO PLOT #################
### DEFINE GENES FOR LABELLING ---
genes <- c("ADRB1","GALT","KLK3","MYT1L","PBX1","PTGES3","TMPRSS2","AR")

### VOLCANO PLOT ---
p <- getVolcanoPlot(df=dat, genes=genes,
        logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=TRUE, 
        fc.upperbound=5, nlogp.upperbound=10)

### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_plots, "extended_data_fig_9.pdf")
pdf(file.plot, height=3, width=3)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()	

