##### Create read file ##### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] #Tissue 

setwd("./12_TWAS/5_all_combined/4_im")

# List files 
files <- list(paste0("AD_", tissue, "_ALL.dat"), paste0("ASTH_", tissue, "_ALL.dat"), paste0("IBD_", tissue, "_ALL.dat"), paste0("MDD_", tissue, "_ALL.dat"), paste0("MS_", tissue, "_ALL.dat"))

# List abbreviations 
trait.names = c("AD", "ASTH", "IBD", "MDD", "MS")

# List effective sample size 
N = c(203311.45,
203321.92,
57637.70,
783651,
35828.01)

# List binary outcomes 
binary=c(T, T, T, T, T)

# Perform read fusion 
IM_genes <- read_fusion(files,trait.names=trait.names,N=N,binary=binary,perm = FALSE)

output = paste0("./6_read_fusion/4_im/IM_", tissue, ".rds")  

saveRDS(IM_genes, file = output)

