######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARY ---
library("stringr")
library("data.table")
library("GenomicRanges")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

dir.chipseq <- file.path("/data1/home/rshrestha/projects/chip_seq/chipseq_haolong_20200706/data/chipseq_AR")

### DEFINE FILES ---
file.metadata <- file.path(dir.wrk, "metadata/metadata.tsv")
file.chipseq_ar <- file.path(dir.chipseq, "AR_LNCaP_Chipseq_Stelloo2018_PMID28925401.bed")

### CHROMOSOMES ---
chromosomes <- paste("chr", c(1:22,"X","Y"), sep="")



###########################################################################################################################################################
### LOAD METADATA ---
metadata <- data.table::fread(file=file.metadata, sep="\t", header=TRUE, nThread=1,  data.table=FALSE, verbose=FALSE)






###########################################################################################################################################################
######### EXTRACT ATAC-SEQ PEAKS #########
### FUNCTION: readPeakFile() ---
readPeakFile <- function(file.peak, chromosomes){
    # LOAD PEAK FILE ---
    peak.bed <- data.table::fread(file=file.peak, sep="\t", header=FALSE, nThread=1, data.table=FALSE, verbose=FALSE)
    peak.bed$V6 <- NULL
    colnames(peak.bed) <- c("chromosome", "start", "end", "peak_id", "int_qvalue", "fold_enrichment_summit", "pvalue", "qvalue", "summit_position_to_peak_start")  

    # SUBSET BY CHROMOSOME ---
    peak.bed <- peak.bed[which(peak.bed$chromosome %in% chromosomes),]

    return(peak.bed)
}

### FUNCTION: getPeakBEDList() ---
getPeakBEDList <- function(sampleids, files.peak, chromosomes){
    # REMOVE DUPLICATE PEAKS ---
    list.atacseq.peaks <- list()

    for(sampleid in sampleids){
        cat("START:", sampleid, "\n", sep="\t")

        # GET PEAKS ---
        dat  <- readPeakFile(file.peak=files.peak[sampleid], chromosomes)
        dat <- subset(dat, select=c("chromosome","start","end"))
        dat <- dat[!duplicated(dat),]
        dat$SampleID <- sampleid
        #dat$FeatureID <- apply(dat, 1, function(x) paste(x[4], x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))

        # CONVERT TO GENOMEIC RANGES ---
        gr <- GenomicRanges::makeGRangesFromDataFrame(df=dat, keep.extra.columns=TRUE)
        gr <- GenomeInfoDb::sortSeqlevels(gr)
        gr <- sort(gr, ignore.strand=TRUE)

        # CONVERT BACK TO DATAFRAME ---
        df <- as.data.frame(gr)
        df <- subset(df, select=c("seqnames","start","end","SampleID"))
        colnames(df) <- c("chr","start","end","SampleID")
        rownames(df) <- NULL

        # STORE DATA TO LIST ---
        list.atacseq.peaks[[sampleid]] <- df

        cat("DONE:", sampleid, "\n", sep="\t")
    }

    return(list.atacseq.peaks)
}

### GET ATAC-SEQ PEAKS ---
files.peak <- sapply(metadata$SAMPLE_ID, function(x) { paste(dir.wrk, "/processed/", sprintf("%s/peaks/%s.filtered.narrowPeak", x, x), sep="") })
list.peaks <- getPeakBEDList(sampleids=metadata$SAMPLE_ID, files.peak, chromosomes)




###########################################################################################################################################################
list.output <- list(metadata=metadata,
                    atacseq_peaks=list.peaks)


### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_atacseq_masterdata.rds")
saveRDS(object=list.output, file=file.rds_atacseq_masterdata)






###########################################################################################################################################################
################# AR ChIP-seq #################
### LOAD CHIP-SEQ FILE ---
chipseq_peak.bed <- data.table::fread(file=file.chipseq_ar, sep="\t", header=FALSE, nThread=1, data.table=FALSE, verbose=FALSE)
colnames(chipseq_peak.bed ) <- c("chr","start","end","name","score")

### SUBSET BY CHROMOSOMES ---
chipseq_peak.bed <- subset(chipseq_peak.bed, chipseq_peak.bed$chr %in% chromosomes)

#### CONVERT TO GR RANGES ---
gr_chipseq_ar <- GenomicRanges::makeGRangesFromDataFrame(df=chipseq_peak.bed, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field="chr", start.field="start", end.field="end")

### SAVE OBJECT TO RDATA FILE ---
file.rds_chipseq_ar <- file.path(dir.reproduce_data, "AR_LNCaP_ChIPseq_Stelloo2018_PMID28925401.rds")
saveRDS(object=gr_chipseq_ar, file=file.rds_chipseq_ar)
