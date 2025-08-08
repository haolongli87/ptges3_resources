######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("ggVennDiagram")
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
file.rds_chipseq_ar <- file.path(dir.reproduce_data, "AR_LNCaP_ChIPseq_Stelloo2018_PMID28925401.rds")


###########################################################################################################################################################
### LOAD CHIP-SEQ PEAKS DATA ---
gr_chipseq_ar <- readRDS(file=file.rds_chipseq_ar)

### LOAD ATAC-SEQ FEATURE COUNTS DATA ---
list.rds_atacseq_feature_counts <- readRDS(file=file.rds_atacseq_feature_counts)

### GET FEATURE ANNOTATION ---
dat_feature_annot <- list.rds_atacseq_feature_counts$feature_annotation

### LOAD ATAC-SEQ DIFF PEAKS DATA ---
list.rds_atacseq_diffpeaks <- readRDS(file=file.rds_atacseq_diffpeaks)
dat_diffpeaks <- list.rds_atacseq_diffpeaks$diffPeaks
dat_diffpeaks_pass <- list.rds_atacseq_diffpeaks$diffPeaks_pass

### MERGE DATA ---
#dat_diffpeaks_annot <- merge(x=dat_diffpeaks, y=dat_feature_annot, by="FeatureID", all.x=TRUE)


### GET CHR/START/END OF FEATUREID ---
dat_diffpeaks_pass$chr <- unlist(lapply(stringr::str_split(dat_diffpeaks_pass$FeatureID, "_"), function(x) x[1]))
dat_diffpeaks_pass$start <- unlist(lapply(stringr::str_split(dat_diffpeaks_pass$FeatureID, "_"), function(x) as.numeric(x[2]) ))
dat_diffpeaks_pass$end <- unlist(lapply(stringr::str_split(dat_diffpeaks_pass$FeatureID, "_"), function(x) as.numeric(x[3]) ))

### CONVERT TO GRANGE OBJ ---
gr_atacseq_diffpeaks <- GenomicRanges::makeGRangesFromDataFrame(df=dat_diffpeaks_pass, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field="chr", start.field="start", end.field="end")



#> length(gr_atacseq_diffpeaks)
#[1] 393
#> length(gr_chipseq_ar)
#[1] 22651


### GET ATAC-SEQ AND CHIPSEQ OVERLAP ---
gr_overlap_chipseq <- IRanges::subsetByOverlaps(gr_chipseq_ar, gr_atacseq_diffpeaks, ignore.strand=TRUE)
gr_overlap_atacseq <- IRanges::subsetByOverlaps(gr_atacseq_diffpeaks, gr_chipseq_ar, ignore.strand=TRUE)


percent_chipseq_overlap <- round((length(gr_overlap_chipseq) / length(gr_chipseq_ar)) * 100, 2)
percent_atacseq_overlap <- round((length(gr_overlap_atacseq) / length(gr_atacseq_diffpeaks)) * 100, 2)


#> length(gr_overlap_chipseq)
#[1] 323
#> length(gr_overlap_atacseq)
#[1] 320



df <- data.frame(Dataset=c("Stelloo2018", "sgPTGES3"),
                Analysis=c("ChIP-seq","ATAC-seq"), 
                TotalPeaks=c( length(gr_chipseq_ar), length(gr_atacseq_diffpeaks) ),
                OverlapPeaks=c( length(gr_overlap_chipseq), length(gr_overlap_atacseq) ),
                OverlapPercent=c( percent_chipseq_overlap, percent_atacseq_overlap ) )

#> df
#   Dataset Analysis TotalPeaks OverlapPeaks OverlapPercent
#1  Stelloo ChIP-seq      22651          323           1.43
#2 sgPTGES3 ATAC-seq        393          320          81.42



### PREPARE DATA FOR PLOT ---
dm1 <- data.frame(Dataset=df$Dataset[1],
                    Variable=c("ARE","NoARE"),
                    Value=c(df$OverlapPercent[1], (100 - df$OverlapPercent[1]) ))

dm2 <- data.frame(Dataset=df$Dataset[2],
                    Variable=c("ARE","NoARE"),
                    Value=c(df$OverlapPercent[2], (100 - df$OverlapPercent[2]) ))

dm1$Variable <- factor(dm1$Variable, levels=c("NoARE","ARE"))
dm2$Variable <- factor(dm2$Variable, levels=c("NoARE","ARE"))

### PLOT
p1 <- ggplot(dm1, aes(x=Value, y=Dataset)) + 
        geom_bar(aes(fill=Variable), stat="identity", color="#000000", width=0.8, size=0.5) + 
        scale_fill_manual(values=c("#bdbdbd","#e31a1c")) +
        theme(
          axis.text.x = element_text(size = 8, color="#000000"),
          axis.text.y = element_text(size = 8, color="#000000"),
          axis.title = element_text(size = 8, color="#000000"),
          plot.title = element_text(size = 8, color="#000000", hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(linewidth=0.4, color="#000000"), 
          strip.text = element_text(size=10, color="#000000"),
          strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
          panel.background = element_rect(fill="#FFFFFF", color="#000000"),
          legend.text = element_text(size = 10, color="#000000"),
          legend.title = element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "bottom") +
      ylab("") +
      xlab("AR ChIP-seq Peaks (in %)") + 
      ggtitle("") 

p2 <- ggplot(dm2, aes(x=Value, y=Dataset)) + 
        geom_bar(aes(fill=Variable), stat="identity", color="#000000", width=0.8, size=0.5) + 
        scale_fill_manual(values=c("#bdbdbd","#e31a1c")) +
        theme(
          axis.text.x = element_text(size = 8, color="#000000"),
          axis.text.y = element_text(size = 8, color="#000000"),
          axis.title = element_text(size = 8, color="#000000"),
          plot.title = element_text(size = 8, color="#000000", hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(linewidth=0.4, color="#000000"), 
          strip.text = element_text(size=10, color="#000000"),
          strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
          panel.background = element_rect(fill="#FFFFFF", color="#000000"),
          legend.text = element_text(size = 10, color="#000000"),
          legend.title = element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "bottom") +
      ylab("") +
      xlab("sgPTGES3 Differentially Accessible ATAC-seq Peaks (in %)") + 
      ggtitle("") 


### GENERATE PLOT ---
file.plot <- file.path(dir.reproduce_plots, "extended_data_fig_9b.pdf")
pdf(file.plot, height=3, width=5)
    grid.arrange(p2, p1, nrow=2, ncol=1)
dev.off()


