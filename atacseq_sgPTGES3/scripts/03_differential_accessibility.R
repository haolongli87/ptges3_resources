######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("GenomicRanges")
library("DESeq2")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)


### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_masterdata.rds")
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_feature_counts.rds")


###########################################################################################################################################################
### FUNCTION: getDiffPeaks() ---
getDiffPeaks <- function(obj, group_treat, group_ctrl){
    item <- c("CONDITION", group_treat, group_ctrl)
    res <- DESeq2::results(obj, contrast=item, independentFiltering=FALSE, pAdjustMethod="BH")

    index0 <- which(res$pvalue == 0)
    if(length(index0) != 0){
        res$pvalue[index0] <- .Machine$double.eps #smallest value
        res$padj <- p.adjust(res$pvalue, method = "BH", n = length(res$pvalue))
    }


    res <- as.data.frame(res)
    res <- res[order(res$pvalue, decreasing=FALSE),]

    dat_summary <- cbind(FeatureID=rownames(res), res)
    rownames(dat_summary) <- NULL

    return(dat_summary)
}

### FUNCTION: getBED() ---
getBED <- function(feature_ids){
    bed <- data.frame(chr=unlist( lapply( stringr::str_split(feature_ids, "_"), function(x) x[1]) ),
                        start= unlist( lapply( stringr::str_split(feature_ids, "_"), function(x) as.numeric(x[2]) ) ),
                        end= unlist( lapply( stringr::str_split(feature_ids, "_"), function(x) as.numeric(x[3]) ) ))
    return(bed)
}

###########################################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### EXTRACT METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata
rownames(metadata) <- metadata$SAMPLE_ID
metadata$CONDITION <- as.factor(metadata$CONDITION)


### LOAD ATAC-SEQ FEATURE COUNTS DATA ---
list.rds_atacseq_feature_counts <- readRDS(file=file.rds_atacseq_feature_counts)

### GET COUNT MATRIX ---
dat_feature_counts_raw <- list.rds_atacseq_feature_counts$feature_counts_raw

### GET DESEQ2 OBJECT ---
dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat_feature_counts_raw, colData=metadata, design= ~ CONDITION)
dds <- DESeq2::DESeq(dds)





###########################################################################################################################################################
### COMPUTE DIFFERENTIAL ACCESSIBILITY ---
dat_diffPeaks <- getDiffPeaks(obj=dds, group_treat="sgPTGES3", group_ctrl="sgGAL4")

#> nrow(dat_diffPeaks)
#[1] 82127


### FILTER BY FDR AND FOLDCHANGE ---
dat_diffPeaks_pass_up <- dat_diffPeaks[which( (dat_diffPeaks$padj <= 0.05) & (dat_diffPeaks$log2FoldChange > 1) ),]
dat_diffPeaks_pass_dn <- dat_diffPeaks[which( (dat_diffPeaks$padj <= 0.05) & (dat_diffPeaks$log2FoldChange < -1)),]

dat_diffPeaks_pass <- rbind(dat_diffPeaks_pass_up, dat_diffPeaks_pass_dn)

#> nrow(dat_diffPeaks_pass)
#[1] 393

#> nrow(dat_diffPeaks_pass_up)
#[1] 56

#> nrow(dat_diffPeaks_pass_dn)
#[1] 337

#> (nrow(dat_diffPeaks_pass)/nrow(dat_diffPeaks))* 100
#[1] 0.4785272

### PREPARE BED FILES ---
bed_diff <- getBED(feature_ids=dat_diffPeaks_pass$FeatureID)
bed_non_diff <- getBED(feature_ids= setdiff(dat_diffPeaks$FeatureID, dat_diffPeaks_pass$FeatureID) )

write.table(bed_diff, file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bed_non_diff, file.path(dir.reproduce_data, "sgPTGES3_atacseq_non_diffpeaks.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### ADD TO LIST ---
list.output <- list(diffPeaks=dat_diffPeaks,
                        diffPeaks_pass=dat_diffPeaks_pass,
                        diffPeaks_pass_up=dat_diffPeaks_pass_up,
                        diffPeaks_pass_dn=dat_diffPeaks_pass_dn)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_diffpeaks.rds")
saveRDS(object=list.output, file=file.rds_atacseq_diffpeaks)




###########################################################################################################################################################
### COMPUTE DIFFERENTIAL ACCESSIBILITY ---
dat_diffPeaks_2 <- getDiffPeaks(obj=dds, group_treat="sgAR", group_ctrl="sgGAL4")

### FILTER BY FDR AND FOLDCHANGE ---
dat_diffPeaks_2_pass_up <- dat_diffPeaks_2[which( (dat_diffPeaks_2$padj <= 0.05) & (dat_diffPeaks_2$log2FoldChange > 1) ),]
dat_diffPeaks_2_pass_dn <- dat_diffPeaks_2[which( (dat_diffPeaks_2$padj <= 0.05) & (dat_diffPeaks_2$log2FoldChange < -1)),]

dat_diffPeaks_2_pass <- rbind(dat_diffPeaks_2_pass_up, dat_diffPeaks_2_pass_dn)

#> nrow(dat_diffPeaks_2_pass_up)
#[1] 349
#> nrow(dat_diffPeaks_2_pass_dn)
#[1] 934

#> (nrow(dat_diffPeaks_pass)/nrow(dat_diffPeaks))* 100
#[1] 1.562215


### ADD TO LIST ---
list.output <- list(diffPeaks=dat_diffPeaks_2,
                        diffPeaks_pass=dat_diffPeaks_2_pass,
                        diffPeaks_pass_up=dat_diffPeaks_2_pass_up,
                        diffPeaks_pass_dn=dat_diffPeaks_2_pass_dn)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgAR_atacseq_diffpeaks.rds")
saveRDS(object=list.output, file=file.rds_atacseq_diffpeaks)




###########################################################################################################################################################
### COMPUTE DIFFERENTIAL ACCESSIBILITY ---
dat_diffPeaks_3 <- getDiffPeaks(obj=dds, group_treat="sgPTGES3", group_ctrl="sgAR")

### FILTER BY FDR AND FOLDCHANGE ---
dat_diffPeaks_3_pass_up <- dat_diffPeaks_3[which( (dat_diffPeaks_3$padj <= 0.05) & (dat_diffPeaks_3$log2FoldChange > 1) ),]
dat_diffPeaks_3_pass_dn <- dat_diffPeaks_3[which( (dat_diffPeaks_3$padj <= 0.05) & (dat_diffPeaks_3$log2FoldChange < -1)),]

dat_diffPeaks_3_pass <- rbind(dat_diffPeaks_3_pass_up, dat_diffPeaks_3_pass_dn)

#> nrow(dat_diffPeaks_3_pass_up)
#[1] 0
#> nrow(dat_diffPeaks_3_pass_dn)
#[1] 6

### ADD TO LIST ---
list.output <- list(diffPeaks=dat_diffPeaks_3,
                        diffPeaks_pass=dat_diffPeaks_3_pass,
                        diffPeaks_pass_up=dat_diffPeaks_3_pass_up,
                        diffPeaks_pass_dn=dat_diffPeaks_3_pass_dn)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_sgAR_atacseq_diffpeaks.rds")
saveRDS(object=list.output, file=file.rds_atacseq_diffpeaks)

#> merge(dat_diffPeaks_3_pass, annot, by = "FeatureID")
#                 FeatureID  baseMean log2FoldChange     lfcSE      stat
#1     chr1_9835764_9836663  72.11504      -1.813572 0.3416393 -5.308440
#2  chr12_56687311_56688901 237.50105      -1.234519 0.2182443 -5.656590
#3 chr5_140322993_140323885 117.93011      -1.415343 0.2772720 -5.104530
#4   chr5_90473466_90475254 333.54694      -1.067455 0.2279622 -4.682596
#5   chr6_80003994_80004957 224.23903      -1.381804 0.2232511 -6.189459
#6   chr9_34646294_34647284 154.50416      -1.962258 0.2724465 -7.202360
#        pvalue         padj           Feature distanceToTSS             geneId
#1 1.105675e-07 2.270145e-03 Distal Intergenic         -7493  ENSG00000280113.2
#2 1.544104e-08 4.227088e-04          Promoter             0 ENSG00000110958.15
#3 3.316178e-07 5.446955e-03 Distal Intergenic         11245  ENSG00000113070.7
#4 2.832641e-06 3.877272e-02          Promoter             0 ENSG00000113356.11
#5 6.037110e-10 2.479049e-05          Promoter             0  ENSG00000112742.9
#6 5.917905e-13 4.860197e-08          Promoter             0 ENSG00000213930.11
#       transcriptId         Gene       GeneType
#1 ENST00000639753.1 RP11-84A14.7            TEC
#2 ENST00000537473.2       PTGES3 protein_coding
#3 ENST00000482211.2        HBEGF protein_coding
#4 ENST00000505345.5       POLR3G protein_coding
#5 ENST00000230510.7          TTK protein_coding
#6 ENST00000450095.6         GALT protein_coding
