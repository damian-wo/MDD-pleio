##### Create read file ##### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] #Tissue 

setwd("./12_TWAS/5_all_combined/3_gastro")

# List files 
files <- list(paste0("GORD_", tissue, "_ALL.dat"), paste0("GSD_", tissue, "_ALL.dat"), paste0("IBS_", tissue, "_ALL.dat"), paste0("MDD_", tissue, "_ALL.dat"), paste0("PUD_", tissue, "_ALL.dat"))

# List abbreviations 
trait.names = c("GORD", "GSD", "IBS", "MDD", "PUD")

# List effective sample size 
N = c(193040.52,
160663.48,
188382.95,
783651,
64229.29
)

# List binary outcomes 
binary=c(T, T, T, T, T)

# Perform read fusion 
Gastro_genes <- read_fusion(files,trait.names=trait.names,N=N,binary=binary,perm = FALSE)

output = paste0("./12_TWAS/6_read_fusion/3_gastro/Gastro_", tissue, ".rds")  

saveRDS(Gastro_genes, file = output)

