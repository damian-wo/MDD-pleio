##### Create read file ##### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] #Tissue 

setwd("./12_TWAS/5_all_combined/1_cvd")

# List files 
files <- list(paste0("AF_", tissue, "_ALL.dat"), paste0("CAD_", tissue, "_ALL.dat"), paste0("HF_", tissue, "_ALL.dat"), paste0("MDD_", tissue, "_ALL.dat"), paste0("STK_", tissue, "_ALL.dat"))

# List abbreviations 
trait.names = c("AF", "CAD", "HF", "MDD", "STK")

# List effective sample size 
N = c(222050.36,
530185.43,
162461.31,
783651,
117522.00)

# List binary outcomes 
binary=c(T, T, T, T, T)

# Perform read fusion 
CVD_genes <- read_fusion(files,trait.names=trait.names,N=N,binary=binary, perm = FALSE)

output = paste0("./12_TWAS/6_read_fusion/1_cvd/CVD_", tissue, ".rds")  

saveRDS(CVD_genes, file = output)

