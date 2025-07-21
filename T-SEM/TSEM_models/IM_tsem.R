##### T_SEM Script ###### 

require(GenomicSEM)
library(devtools) 

args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] 

setwd("./12_TWAS/T-SEM")

load("LDSCoutput_IM.RData")

genes <-  readRDS(paste0("./6_read_fusion/4_im/IM_", tissue, ".rds"))

model1 <- "
F_IM =~ ASTH + AD + IBD + MS 
F_IM_MDD =~ a*F_IM + a*MDD
ASTH ~~ AD
F_IM ~~ b*F_IM
b > 0.001
F_IM_MDD ~ Gene
"

results1 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model1, sub=c("F_IM_MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out1 <- paste0("./4_im/IM_", tissue, "TSEM")

write.csv(results1, file = out1, row.names = FALSE, quote = TRUE)

model2 <- "
F_IM =~ ASTH + AD + IBD + MS 
F_IM_MDD =~ a*F_IM + a*MDD
ASTH ~~ AD
F_IM ~~ b*F_IM
b > 0.001
F_IM + MDD ~ Gene 
"

results2 <- userGWAS(covstruc = LDSCoutput, SNPs = genes, model = model2, sub=c("MDD ~ Gene"), estimation = "DWLS", parallel = FALSE, TWAS = TRUE)

out2 <- paste0("./4_im/IM_", tissue, "_independent_TSEM")

write.csv(results2, file = out2, row.names = FALSE, quote = TRUE)
