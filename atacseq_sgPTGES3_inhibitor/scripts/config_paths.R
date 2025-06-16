######################################################################################
## COPYRIGHT UNIVERSITY OF CALIFORNIA, SAN FRANCISCO (UCSF) 2025                    ##
## Li H. et al. Genome-scale CRISPR screens identify PTGES3 as a direct modulator   ##
##              of AR function in advanced prostate cancer. Nature Genetics. 2025   ##
######################################################################################

cat("###############################################################################################################\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "LOADAING PATHS ... ", "\n", sep=" ")

### DEFINE PATHS ---
dir.wrk <- file.path("/data1/projects/feng_cells_ATAC/2021_PTGES3_haolong")
dir.reproduce <- file.path(dir.wrk, "analysis/reproduce_2024")  
dir.project <- file.path(dir.reproduce, "atacseq_sgPTGES3_inhibitor") # CHANGE THIS AS REQUIRED: PATH TO RE-PRODUCTION DIRECTORY
dir.reproduce_data <- file.path(dir.project, "data") # PATH TO DIRECTORY CONTAINING ALL REQUIRED DATA 
dir.reproduce_plots <- file.path(dir.project, "plots") # PATH TO DIRECTORY CONTAINING ALL OUTPUT PLOTS 
dir.reproduce_scripts <- file.path(dir.reproduce, "scripts")  # PATH TO DIRECTORY CONTAINING SCRIPTS

#### SET WORKING DIRECTORY TO PRODUCTION DIRECTORY ---------------------------------------
setwd(dir.project)
cat(format(Sys.time(), "%b %d %X"), "WORKING DIRECTORY IS SET TO:", getwd(), "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_data:", dir.reproduce_data, "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_plots:", dir.reproduce_plots, "\n", sep=" ")
cat(format(Sys.time(), "%b %d %X"), "dir.reproduce_scripts:", dir.reproduce_scripts, "\n", sep=" ")

cat(format(Sys.time(), "%b %d %X"), "DONE!", "\n", sep=" ")
cat("###############################################################################################################\n", sep=" ")
