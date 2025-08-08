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
file.dat <- file.path(dir.reproduce_data, "22rv1_MS_data_normalized.tsv")
file.annot <- file.path(dir.reproduce_data, "protein_gene_conversion_uniprot.tsv")

#####################################################################################
### FUNCTION: get.normalizeQuantile() ---
get.normalizeQuantile <- function(dat){
    require("aroma.light")
	dQNorm <- data.frame(aroma.light::normalizeQuantile(as.matrix(dat)))
    colnames(dQNorm) <- colnames(dat)
	return(dQNorm)
}	

### FUNCTION: get_ttest() ----
get_ttest <- function(dat, class1, class2){
	cat("COMPUTE T-TEST ...", "\n", sep="")
	
	#Get class indices
	index1 <- which(colnames(dat) %in% class1)
	index2 <- which(colnames(dat) %in% class2)
	
	# Check if the data points are constant ----
	exp <- dat[,c(index1, index2)]
	len <- apply(exp, 1, function(x) length(unique(as.numeric(x))))
	del.index <- which(len == 1)

	if(length(del.index) != 0){
		dat <- dat[-del.index,]
	}

	# Check if each group has at least 2 samples represented ---
	length_not_na.A <- as.numeric(apply(dat[,index1], 1, function(x) length(which(!is.na(as.numeric(x))))))
	length_not_na.B <- as.numeric(apply(dat[,index2], 1, function(x) length(which(!is.na(as.numeric(x))))))
	del.index_A <- which(length_not_na.A < 2)
	del.index_B <- which(length_not_na.B < 2)

	del.index_AB <- unique(c(del.index_A, del.index_B))
	
	if(length(del.index_AB) != 0){
		dat <- dat[-del.index_AB,]
	}
	
	# Compute Median and StdDev for the Sample Class
	median1 <- apply(as.matrix(dat)[,index1], 1, median, na.rm=TRUE)
	median2 <- apply(as.matrix(dat)[,index2], 1, median, na.rm=TRUE)
	cat("MEDIAN COMPUTED ...", "\n", sep="")
	
	stdev1 <- apply(as.matrix(dat)[,index1], 1, sd, na.rm=TRUE)
	stdev2 <- apply(as.matrix(dat)[,index2], 1, sd, na.rm=TRUE)
	cat("STANDARD DEVIATION COMPUTED ...", "\n", sep="")
	
	# Compute FoldChange
	foldchange <- ifelse(median2 > median1, 2^median2/2^median1, -2^median1/2^median2)
	cat("FOLD CHANGE COMPUTED ...", "\n", sep="")
	
	cat("PERFORMING T-TEST ...", "\n", sep="")
	# FUNCTION: func_ttest() ---
	func_ttest <- function(x, index1, index2){
		res_ttest <- t.test(x[index1], x[index2], alternative ="two.sided", conf.level = 0.95)
		return(res_ttest)
	}

	ttest.list <- suppressWarnings(apply(as.matrix(dat), 1, function(x) func_ttest(x, index1, index2)))
	cat("T-TEST COMPUTED ...", "\n", sep="")
	
	pvalue <- suppressWarnings(do.call(c,lapply(ttest.list,function(x) x$p.value)))

	index0 <- which(pvalue == 0)
	if(length(index0) != 0){
		pvalue[index0] <- .Machine$double.eps #smallest value
	}

	fdr <- suppressWarnings(p.adjust(pvalue, method = "BH", n = length(pvalue)))
	cat("P-VALUE COMPUTED ...", "\n", sep="")

	dat.summary <- data.frame(Gene=rownames(dat), 
								MedianA=median1, StdevA=stdev1, 
								MedianB=median2, StdevB=stdev2, 
								FoldChange=foldchange, 
								pvalue=pvalue, fdr=fdr)
	rownames(dat.summary) <- c(1:nrow(dat))
	dat.summary$Gene <- as.character(dat.summary$Gene)
	dat.summary <- dat.summary[order(dat.summary$pvalue, decreasing=FALSE),]
	rownames(dat.summary) <- NULL
	cat("DATA SUMMARY COMPUTED ...", "\n", sep="")
		
	return(dat.summary)
}

#####################################################################################
### LOAD DESIGN TABLE ---
annot <- data.table::fread(file=file.annot, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)

### LOAD DESIGN TABLE ---
dat <- data.table::fread(file=file.dat, sep="\t", header=TRUE, nThread=1, data.table=FALSE, verbose=FALSE)
rownames(dat) <- dat$ProteinID
dat$ProteinID <- NULL


#####################################################################################
### QUANTILE NORMALIZATION ---
dQnorm <- get.normalizeQuantile(dat)

### SAVE OBJECT TO RDATA FILE ---
list.output <- list(expr_norm=dQnorm, annotation=annot)
file.rds_masspec <- file.path(dir.reproduce_data, "sgPTGES3_masspec_expression.rds")
saveRDS(object=list.output, file=file.rds_masspec)



#####################################################################################
### T-TEST ----
df <- get_ttest(dat=dQnorm, class1=c("sgGAL4_1","sgGAL4_2","sgGAL4_3"), class2=c("sgPTGES3_1","sgPTGES3_2","sgPTGES3_3") )
colnames(df)[1] <- "UniProtID" 
df$Protein <- NA
df <- subset(df, select=c("UniProtID","Protein","MedianA","StdevA","MedianB","StdevB","FoldChange","pvalue","fdr"))

### ADD PROTEIN NAME ---
for(i in 1:nrow(df)){
	index <- which(annot$UniprotID == df$UniProtID[i])

	if(length(index) != 0){
		df$Protein[i] <- annot$Protein[index[1]]  
	}
}


### WRITE OUTPUT ---
file.output <- file.path(dir.reproduce_data, "PTGES3ko_masspec_results_ttest.tsv")
write.table(df, file.output, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
