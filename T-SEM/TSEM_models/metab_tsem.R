##### T_SEM Script ###### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] 

setwd("./12_TWAS/T_SEM")

load("LDSCoutput_metab.RData")

genes <-  readRDS(paste0("./6_read_fusion/2_metab/Metab_", tissue, ".rds"))

model1 <- "
F_metab =~ MSYN + T2D + TG + HDL + GSD 
F_metab_MDD =~ a*MDD + a*F_metab

b >0.001
MSYN ~~ b*MSYN 

F_metab_MDD ~ Gene
"

results1 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model1, sub=c("F_metab_MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out1 <- paste0("./2_metab/Metab_", tissue, "TSEM")

write.csv(results1, file = out1, row.names = FALSE, quote = TRUE)

model2 <- "
F_metab =~ MSYN + T2D + TG + HDL + GSD 
F_metab_MDD =~ a*MDD + a*F_metab

b >0.001
MSYN ~~ b*MSYN 

F_metab + MDD ~ Gene
"

results2 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model2, sub=c("MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out2 <- paste0("./2_metab/Metab_", tissue, "_independent_TSEM")

write.csv(results2, file = out2, row.names = FALSE, quote = TRUE)

