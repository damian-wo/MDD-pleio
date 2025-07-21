##### Create read file ##### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] #Tissue 

setwd(".12_TWAS/5_all_combined/2_metab/")

# List files 
files <- list(paste0("GSD_", tissue, "_ALL.dat"), paste0("HDL_", tissue, "_ALL.dat"), paste0("MDD_", tissue, "_ALL.dat"), paste0("MSYN_", tissue, "_ALL.dat"), paste0("T2D_", tissue, "_ALL.dat"), paste0("TG_", tissue, "_ALL.dat")) 

# List abbreviations 
trait.names = c("GSD", "HDL", "MDD", "MSYN", "T2D", "TG")

# List effective sample size 
N = c(160663.48,
1320016,
783651,
189772.81,
751754.69,
1320016)

# List binary outcomes 
binary=c(T, F, T, T, T, F)

# Perform read fusion 
Metab_genes <- read_fusion(files, trait.names=trait.names, N=N, binary=binary, perm = FALSE)

output = paste0("./6_read_fusion/2_metab/Metab_", tissue, ".rds")  

saveRDS(Metab_genes, file = output)
