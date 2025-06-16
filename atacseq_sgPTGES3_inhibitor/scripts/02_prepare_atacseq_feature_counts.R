######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

### DEFINE LIBRARIES ---
library("stringr")
library("data.table")
library("ChIPQC")
library("GenomicRanges")
library("Rsubread")
library("DESeq2")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

### DEFINE PATH ---
dir.project <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong/analysis/reproduce_2024/atacseq_sgPTGES3_inhibitor") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.project, "scripts/config_paths.R") 
source(file.config_path)

### DEFINE PATH ---
dir.wrk <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong") # CHANGE THE PATH AS REQUIRED
dir.reproduce <- file.path(dir.wrk, "analysis/reproduce_2024") # CHANGE THE PATH AS REQUIRED
file.config_path <- file.path(dir.reproduce, "scripts/config_paths.R") 
source(file.config_path)


### DEFINE FILES ---
file.rds_atacseq_masterdata <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_masterdata.rds")
file.reference_annotation <- file.path("/data1/projects/WCDT_atac_2020/reference/annotation_data/gene_annotation/gencode_v28/gencode.v28.annotation.gtf.gz")


###########################################################################################################################################################
### FUNCTION: parseGTF() --
parseGTF <- function(file.gtf, feature.type){
    # CHROMOSOMES ---
    chromosome <- paste("chr", c(1:22,"X","Y"), sep="")

    # LOAD GTF FILE ---
    #annot <- data.table::fread(file.gtf, data.table=FALSE, stringsAsFactors=FALSE, skip=5)
    annot <- data.table::fread(file=file.gtf, sep="\t", header=FALSE, nThread=50, data.table=TRUE, stringsAsFactors=FALSE, skip=5, verbose=FALSE)
    annot <- subset(annot, annot$V1 %in% chromosome)

    # PARSE GTF COLUMNS ---
    if(feature.type == "gene"){    
        annot <- subset(annot, annot$V3 == "gene")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$V3 == "exon")
    }else if(feature.type == "transcript"){
        annot <- subset(annot, annot$V3 == "transcript")
    }

    list.annot <- noquote(stringr::str_split(annot$V9, "; "))
    annot$EnsemblID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_id")]))), " "), function(x) x[2]))
    annot$Gene <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_name")]))), " "), function(x) x[2]))
    annot$GeneType <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "gene_type")]))), " "), function(x) x[2]))

    # TRIM DATA ---
    if(feature.type == "gene"){ 
        annot <- subset(annot, annot$GeneType == "protein_coding")   
        items <- c("V1","V4","V5","Gene","EnsemblID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblID","Strand")
    }else if(feature.type == "exon"){
        annot <- subset(annot, annot$GeneType == "protein_coding")
        list.annot <- stringr::str_split(annot$V9, "; ")
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        annot$ExonID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "exon_id")]))), " "), function(x) x[2]))  
        items <- c("V1","V4","V5","Gene","EnsemblID","TranscriptID","ExonID","V7")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","GeneID","TranscriptID","ExonID","Strand")
    }else if(feature.type == "transcript"){
        annot$TranscriptID <- unlist(lapply(stringr::str_split(unlist(lapply(list.annot, function(x) gsub('"', '', x[stringr::str_detect(x, "transcript_id")]))), " "), function(x) x[2]))
        items <- c("V1","V4","V5","Gene","EnsemblID","V7","TranscriptID","GeneType")
        df <- subset(annot, select=items)
        colnames(df) <- c("Chrom","Start","End","Gene","EnsemblGeneID","Strand","EnsemblTranscriptID","GeneType")
    }        

    return(df)
}

### FUNCTION: getNonOverlapPeaks() ---
# Compile all ATAC-seq peaks from all samples and then find a unique list of non-overlapping genomic regions
getNonOverlapPeaks <- function(list.peaks, sampleids){
    list.peaks_gr <- lapply(list.peaks, ChIPQC:::GetGRanges, simple = TRUE)
    names(list.peaks_gr) <- sampleids
    
    list.peaks_gr <- GenomicRanges::GRangesList(list.peaks_gr)   
    reduced <- GenomicRanges::reduce(unlist(list.peaks_gr))
    consensusIDs <- paste("consensus", seq(1, length(reduced)), sep="_")
    mcols(reduced) <- do.call(cbind, lapply(list.peaks_gr, function(x) (reduced %over% x) + 0))

    reducedConsensus <- reduced
    mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
    consensusIDs <- paste("consensus", seq(1, length(reducedConsensus)), sep="_")

    return(reducedConsensus)
}

### FUNCTION: getATACseqFeatureCount() ---
getATACseqFeatureCount <- function(files.bam, consensusToCount, sampleids, nthreads){
    # GET REGION TO COUNT ---
    regionsToCount <- data.frame(GeneID = paste(
                                            GenomicRanges::seqnames(consensusToCount), 
                                            GenomicRanges::start(consensusToCount), 
                                            GenomicRanges::end(consensusToCount), sep = "_"), 
                                Chr = GenomicRanges::seqnames(consensusToCount), 
                                Start = GenomicRanges::start(consensusToCount), 
                                End = GenomicRanges::end(consensusToCount), 
                                Strand = GenomicRanges::strand(consensusToCount))

    # COUNT READS
    fcResults <- Rsubread::featureCounts(
                                files=files.bam, 
                                annot.ext=regionsToCount, 
                                isPairedEnd=TRUE, 
                                countMultiMappingReads=FALSE, 
                                maxFragLength=100,
                                autosort=TRUE,
                                nthreads=nthreads,
                                verbose=FALSE)
  
    # COLLECT THE RAW COUNTS ---
    counts <- fcResults$counts
    colnames(counts) <- sampleids
    counts <- cbind.data.frame(feature_id=rownames(counts), counts)
    rownames(counts) <- NULL

    return(counts)
}

### FUNCTION: annotPeaksbyRange() ---
annotPeaksbyRange <- function(gr_peaks, genecode.txdb){
    # ANNOTATE PEAKS ---
    annot_peaks <- ChIPseeker::annotatePeak(peak=gr_peaks, TxDb=genecode.txdb, verbose=TRUE)
    annot_peaks <- as.data.frame(annot_peaks)

    # ADD FEATURES ---
    annot_peaks$Feature <- ""
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (<=1kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (1-2kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "Promoter (2-3kb)")] <- "Promoter"
    annot_peaks$Feature[which(annot_peaks$annotation == "5' UTR")] <- "5' UTR"
    annot_peaks$Feature[which(annot_peaks$annotation == "3' UTR")] <- "3' UTR"

    exon <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Exon") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == TRUE)])] <- "Exon"
    annot_peaks$Feature[which(annot_peaks$annotation %in% exon[which(str_detect(exon, "exon 1 of") == FALSE)])] <- "Exon"

    intn <- annot_peaks$annotation[which(stringr::str_detect(annot_peaks$annotation, "Intron") == TRUE)]
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == TRUE)])] <- "Intron"
    annot_peaks$Feature[which(annot_peaks$annotation %in% intn[which(str_detect(intn, "intron 1 of") == FALSE)])] <- "Intron"

    annot_peaks$Feature[which(stringr::str_detect(annot_peaks$annotation, "Downstream") == TRUE)] <- "Downstream (<=300)"
    annot_peaks$Feature[which(annot_peaks$annotation == "Distal Intergenic")] <- "Distal Intergenic"

    # ADD FEATURE IDS ---
    annot_peaks$FeatureID <- apply(annot_peaks, 1, function(x) paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep="_"))

    return(annot_peaks)
}

### FUNCTION: addGenes() ---
addGenes <- function(annot, dat){
    pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
    dat$Gene <- ""
    dat$GeneType <- ""
    n <- nrow(dat)

    # ADD GENE SYBBOLS ---
    for(i in 1:n){
        index <- which(annot$EnsemblTranscriptID == dat$transcriptId[i])
        if(length(index) != 0){
            dat$Gene[i] <- annot$Gene[index]
            dat$GeneType[i] <- annot$GeneType[index]

            # UPDATE PROGRESS BAR ---
            setTxtProgressBar(pb, i)
        }
    }

    return(dat)
}

###########################################################################################################################################################
### LOAD ATAC-SEQ MASTER DATA ---
list.rds_atacseq_masterdata <- readRDS(file=file.rds_atacseq_masterdata)

### EXTRACT METADATA ---
metadata <- list.rds_atacseq_masterdata$metadata
rownames(metadata) <- metadata$SAMPLE_ID
metadata$CONDITION <- as.factor(metadata$CONDITION)

### DEFINE ATAC-SEQ BAM FILES ---
files.bam <- sapply(metadata$SAMPLE_ID, function(x) { paste(dir.wrk, "/2023_PTGES3i_atacseq/processed/", sprintf("%s/alignment/%s", x, x), ".bowtie.sorted.nodup.shifted.bam", sep="") })

### EXTRACT ATAC-SEQ PEAKS ---
list.atacseq_peaks <- list.rds_atacseq_masterdata$atacseq_peaks


###########################################################################################################################################################
### GET NON-OVERLAPPING CONSENSUS PEAKS ---
gr_consensus_peaks <- getNonOverlapPeaks(list.peaks=list.atacseq_peaks, sampleids=metadata$SAMPLE_ID)
gr_consensus_peaks <- GenomeInfoDb::sortSeqlevels(gr_consensus_peaks)
gr_consensus_peaks <- sort(gr_consensus_peaks, ignore.strand=TRUE)

### RETAIN PEAKS THAT ARE DETEDCTED IN AT LEAST 2 SAMPLES ---
occurrences <- GenomicRanges::elementMetadata(gr_consensus_peaks) %>% 
                                    as.data.frame %>% 
                                    dplyr::select(-consensusIDs) %>% 
                                    rowSums

gr_consensus_peaks <- gr_consensus_peaks[occurrences >= 2,]

### GET ATAC-SEQ FEATURE READ COUNTS ---
dat_feat_counts <- getATACseqFeatureCount(files.bam, consensusToCount=gr_consensus_peaks, sampleids=metadata$SAMPLE_ID, nthreads=50)
rownames(dat_feat_counts) <- dat_feat_counts$feature_id
dat_feat_counts$feature_id <- NULL

### GET DESEQ2 OBJECT ---
dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat_feat_counts, colData=metadata, design= ~ CONDITION)
dds <- DESeq2::DESeq(dds)

### NORMALIZE DATA ---
d_rlog <- DESeq2::rlog(dds)
dat_feat_counts_norm <- SummarizedExperiment::assay(d_rlog)





###########################################################################################################################################################
################ FEATURE ANNOTATION ################ 
### GET ANNOTATION DATA ---
dat_reference_annotation <- parseGTF(file.gtf=file.reference_annotation, feature.type="transcript")

### LOAD GENE ANNOTATION DATABASE FROM GFF/GTF ---
genecode.txdb <- GenomicFeatures::makeTxDbFromGFF(file=file.reference_annotation, format="gtf", organism="Homo sapiens", chrominfo=seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene))

### ANNOTATE PEAKS ---
dat_feature_annotation <- annotPeaksbyRange(gr_peaks=gr_consensus_peaks, genecode.txdb)
dat_feature_annotation <- subset(dat_feature_annotation, select=c("FeatureID","Feature","distanceToTSS","geneId","transcriptId"))
dat_feature_annotation <- addGenes(annot=dat_reference_annotation, dat=dat_feature_annotation)




###########################################################################################################################################################
### ADD TO LIST ---
list.output <- list(feature_counts_raw=dat_feat_counts,
                        feature_counts_norm=dat_feat_counts_norm,
                        feature_annotation=dat_feature_annotation)

### SAVE OBJECT TO RDATA FILE ---
file.rds_atacseq_feature_counts <- file.path(dir.reproduce_data, "sgPTGES3_inhibitor_atacseq_feature_counts.rds")
saveRDS(object=list.output, file=file.rds_atacseq_feature_counts)
