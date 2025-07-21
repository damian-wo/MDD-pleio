##### T_SEM Script ###### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] 

setwd("./12_TWAS/T_SEM")

load("LDSCoutput_CVD.RData")

genes <-  readRDS(paste0("./6_read_fusion/1_cvd/CVD_", tissue, ".rds"))

model1 <- "
F_CVD =~ HF + CAD + STK + AF
F_CVD_MDD =~ a*MDD + a*F_CVD
F_CVD_MDD ~ Gene
"

results1 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model1, sub=c("F_CVD_MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out1 <- paste0("./1_cvd/CVD_", tissue, "TSEM")

write.csv(results1, file = out1, row.names = FALSE, quote = TRUE)

model2 <- "
F_CVD =~ HF + CAD + STK + AF
F_CVD_MDD =~ a*MDD + a*F_CVD
F_CVD + MDD ~ Gene
"

results2 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model2, sub=c("MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out2 <- paste0("./1_cvd/CVD_", tissue, "_independent_TSEM")

write.csv(results2, file = out2, row.names = FALSE, quote = TRUE)

