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
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3_inhibitor") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_masterdata.rds")
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_feature_counts.rds")


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
dat_diffPeaks <- getDiffPeaks(obj=dds, group_treat="RAC18", group_ctrl="DMSO")

#> nrow(dat_diffPeaks)
#[1] 95314

### FILTER BY FDR AND FOLDCHANGE ---
dat_diffPeaks_pass_up <- dat_diffPeaks[which( (dat_diffPeaks$padj <= 0.05) & (dat_diffPeaks$log2FoldChange > 1) ),]
dat_diffPeaks_pass_dn <- dat_diffPeaks[which( (dat_diffPeaks$padj <= 0.05) & (dat_diffPeaks$log2FoldChange < -1)),]

dat_diffPeaks_pass <- rbind(dat_diffPeaks_pass_up, dat_diffPeaks_pass_dn)

#> nrow(dat_diffPeaks_pass)
#[1] 471
#> nrow(dat_diffPeaks_pass_up)
#[1] 158
#> nrow(dat_diffPeaks_pass_dn)
#[1] 313


### PREPARE BED FILES ---
bed_diff <- getBED(feature_ids=dat_diffPeaks_pass$FeatureID)
bed_non_diff <- getBED(feature_ids= setdiff(dat_diffPeaks$FeatureID, dat_diffPeaks_pass$FeatureID) )

write.table(bed_diff, file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_diffpeaks.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bed_non_diff, file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_non_diffpeaks.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### ADD TO LIST ---
list.output <- list(diffPeaks=dat_diffPeaks,
                        diffPeaks_pass=dat_diffPeaks_pass,
                        diffPeaks_pass_up=dat_diffPeaks_pass_up,
                        diffPeaks_pass_dn=dat_diffPeaks_pass_dn)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_diffpeaks <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_diffpeaks.rds")
saveRDS(object=list.output, file=file.rds_atacseq_diffpeaks)

