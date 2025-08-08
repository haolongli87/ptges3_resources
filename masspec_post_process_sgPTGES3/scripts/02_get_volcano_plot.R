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

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/masspec_post_process_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.dat <- file.path(dir.reproduce_data, "PTGES3ko_masspec_results_ttest.tsv")


#####################################################################################
### FUNCTION: getVolcanoPlot() ---
getVolcanoPlot <- function(df, compr.name, logfc.threshold, pvalue.threshold, filterby.pvalueadj, fc.upperbound, nlogp.upperbound){
	yval <- -log10(pvalue.threshold)
    xval <- logfc.threshold
	cbPalette1 <- c("#98a3a5","#e85748")

    
    # PVALUE FILTER ---
    if(filterby.pvalueadj == TRUE){
        df$nlogp <- -log10(df$fdr)
    }else{
        df$nlogp <- -log10(df$pvalue)
    }

    # GROUP GENES ---
    df$Group <- "STABLE"
    if(filterby.pvalueadj == TRUE){    
        df$Group[which((df$FoldChange >= logfc.threshold) & (df$padj <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$FoldChange <= -logfc.threshold) & (df$padj <= pvalue.threshold))] <- "VARIABLE"
    }else{
        df$Group[which((df$FoldChange >= logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
        df$Group[which((df$FoldChange <= -logfc.threshold) & (df$pvalue <= pvalue.threshold))] <- "VARIABLE"
    }

    # LABEL GROUP BY PATHWAY GENES TO HIGHLIGHT ---
    #df$Group[which((df$Group == "VARIABLE") & (df$Gene %in% genelist))] <- "ZHIGHLIGHT"

    # ADD GENE NAMES TO LABEL ---
    df$Gene_Label <- NA
    #df$Gene_Label[which(df$Group == "ZHIGHLIGHT")] <- df$Gene[which(df$Group == "ZHIGHLIGHT")]
    df$Gene_Label[which(df$Protein %in% c("DGKZ","PPM1E","SMIM12","ATP6AP1","PDCD4","FASN"))] <- df$Protein[which(df$Protein %in% c("DGKZ","PPM1E","SMIM12","ATP6AP1","PDCD4","FASN"))]
    df$Gene_Label[which(df$Protein %in% c("NF1","VIL1","SMC2","DPYSL2","ALDH7A1","HELLS"))] <- df$Protein[which(df$Protein %in% c("NF1","VIL1","SMC2","DPYSL2","ALDH7A1","HELLS"))]
    df$Gene_Label[which(df$Protein == "PTGES3")] <- "PTGES3"  
    df$Gene_Label[which(df$Protein == "AR")] <- "AR"  

    # ADJUST FOLDCHANGE FOR VIZ ---
    if(max(df$FoldChange) > fc.upperbound){
        df$FoldChange[which(df$FoldChange > fc.upperbound)] <- fc.upperbound 
    }

    if(min(df$FoldChange) < -fc.upperbound){
        df$FoldChange[which(df$FoldChange < -fc.upperbound)] <- -fc.upperbound
    }

    # ADJUST nlogp FOR VIZ ---
    if(max(df$nlogp) > nlogp.upperbound){
        df$nlogp[which(df$nlogp > nlogp.upperbound)] <- nlogp.upperbound 
    }

    # ORDER DATA ---
    df <- df[order(df$Group, decreasing=FALSE),]

    # FACTORIZE ---
    df$Group <- factor(df$Group, levels=c("STABLE","VARIABLE"))

    # PLOT ---
	p <- ggplot(df, aes(x=FoldChange, y=nlogp, label=Gene_Label)) + 
			geom_hline(yintercept = yval, colour="gray70", linetype=4, alpha=0.7) +
			geom_vline(xintercept = xval, colour="gray70", linetype=4, alpha=0.7) +		
			geom_vline(xintercept = -xval, colour="gray70", linetype=4, alpha=0.7) +
			geom_point(aes(fill=Group, color=Group), shape = 21, stroke=0.3, size=1.5, alpha=0.8) +
			scale_fill_manual(values=cbPalette1) +
			scale_color_manual(values=cbPalette1) +
            scale_x_continuous(breaks=seq(-fc.upperbound, fc.upperbound, 1)) +
            scale_y_continuous(breaks=seq(0, nlogp.upperbound, 1)) +
			coord_cartesian(xlim=c(-fc.upperbound, fc.upperbound), ylim=c(0, nlogp.upperbound)) +
			geom_text(aes(label=Gene_Label, alpha=0.3, hjust=1, vjust=-1, lineheight=1.5), color="black", size=2) +
			theme(
				axis.text = element_text(size = 6, color="black"),
				axis.title = element_text(size = 10, color="black"),
				strip.text = element_text(size = 8, color="black"),
				strip.background = element_rect(fill = "white", colour = "white"),
				plot.title = element_text(size = 10, color="black", hjust=0.5),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				axis.ticks = element_line(size=0.2, color="black"),	
				panel.background = element_rect(fill="white", color="black"),
				legend.text = element_text(size = 5, color="#000000"),
			    legend.title = element_blank(),
			    legend.key.size = unit(0.2, "cm"),			
			    legend.position="none") +
			xlab("Fold Change") + 
			ylab("-log10(pvalue)") + 
			ggtitle(compr.name) 	

    return(p)
}

#####################################################################################
### GET VOLCANO PLOT ---
p <- getVolcanoPlot(df=dat, compr.name="sgPTGES3-ko vs sgGAL4", logfc.threshold=1, pvalue.threshold=0.05, filterby.pvalueadj=FALSE, fc.upperbound=5, nlogp.upperbound=6)


# GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_plots, "PTGES3ko_masspec_results_ttest_volcanoplot.pdf")
pdf(file.plot, height=3, width=3)
    grid.arrange(p, nrow=1, ncol=1)
dev.off()

