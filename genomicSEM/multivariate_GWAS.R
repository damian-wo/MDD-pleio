########################################
##### Multivariate GWAS models CVD #####
########################################

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1
output2 <- args[4] #Output 2

## Run multivariate GWAS for CVD factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

## Set Working Directory 
setwd("./1_cvd/mdd_cvd_factor")

## Load in LDSC results and sumstats results
load("./1_cvd/LDSCoutput_CVD.RData")
load("./1_cvd/CVD_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

## Create double factor for shared liability between MDD and CVD factor 

model1 <- "
F_CVD =~ HF + CAD + STK + AF
F_CVD_MDD =~ a*MDD + a*F_CVD
F_CVD_MDD ~ SNP
"
results <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, sub=c("F_CVD_MDD ~ SNP"), estimation = "DWLS", parallel = FALSE, smooth_check = FALSE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

## Run independent pathway

model2 <- "
F_CVD =~ HF + CAD + STK + AF
F_CVD_MDD =~ a*MDD + a*F_CVD
F_CVD + MDD ~ SNP
"
results2 <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model2, estimation = "DWLS", sub=c("MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE)

write.csv(results2, file = paste0(output2), row.names = FALSE, quote = TRUE)

##########################################
##### Multivariate GWAS models Metab #####
##########################################

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1
output2 <- args[4] #Output 2

## Run multivariate GWAS for metab factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

## Set Working Directory 
setwd("./2_metab/mdd_metab_factor/")

## Load in LDSC results and sumstats results
load("./2_metab/LDSCoutput_metab.RData")
load("./2_metab/metab_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

## Create a double factor for shared liability between MDD and metab factor 

model1 <- "
F_metab =~ MSYN + T2D + TG + HDL + GSD 
F_metab_MDD =~ a*MDD + a*F_metab

b >0.001
MSYN ~~ b*MSYN 

F_metab_MDD ~ SNP
"
results <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, estimation = "DWLS", sub=c("F_metab_MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

## Create model for independent pathway 

model2 <- "
F_metab =~ MSYN + T2D + TG + HDL + GSD 
F_metab_MDD =~ a*MDD + a*F_metab

b >0.001
MSYN ~~ b*MSYN 

F_metab + MDD ~ SNP
"
results2 <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model2, estimation = "DWLS", sub=c("MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE)

write.csv(results2, file = paste0(output2), row.names = FALSE, quote = TRUE)

###########################################
##### Multivariate GWAS models Gastro #####
###########################################

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1
output2 <- args[4] #Output 2

## Run multivariate GWAS for gastro factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

## Set Working Directory 
setwd("./3_gastro/mdd_gastro_factor")

## Load in LDSC results and sumstats results
load("./3_gastro/LDSCoutput_gastro.RData")
load("./3_gastro/gastro_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

## Create double factor for shared liability between MDD and gastro factor 

model1 <- '
F_gastro =~ GORD + GSD + IBS + PUD
F_gastro_MDD =~ a*MDD + a*F_gastro
F_gastro_MDD ~ SNP 
'
results <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, estimation = "DWLS", sub=c("F_gastro_MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

## Create model for independent pathway 

model2 <- '
F_gastro =~ GORD + GSD + IBS + PUD
F_gastro_MDD =~ a*MDD + a*F_gastro
F_gastro + MDD ~ SNP 
'
results2 <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model2, estimation = "DWLS", sub=c("MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE)

write.csv(results2, file = paste0(output2), row.names = FALSE, quote = TRUE)

############################################
##### Multivariate GWAS models Immune ######
############################################

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1
output2 <- args[4] #Output 2

## Run multivariate GWAS for IM factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

## Set Working Directory 
setwd("./4_im/mdd_im_factor")

## Load in LDSC results and sumstats results
load("./4_im/LDSCoutput_IM.RData")
load("./4_im/IM_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

## Create double factor for shared liability between MDD and IM factor 

model1 <- "
F_IM =~ ASTH + AD + IBD + MS 
F_IM_MDD =~ a*F_IM + a*MDD
ASTH ~~ AD
F_IM ~~ b*F_IM
b > 0.001
F_IM_MDD ~ SNP
"
results <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, sub=c("F_IM_MDD ~ SNP"), estimation = "DWLS", parallel = FALSE, smooth_check = FALSE, toler = 1e-60)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

model2 <- "
F_IM =~ ASTH + AD + IBD + MS 
F_IM_MDD =~ a*F_IM + a*MDD
ASTH ~~ AD
F_IM ~~ b*F_IM
b > 0.001
F_IM + MDD ~ SNP 
"
results2 <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model2, estimation = "DWLS", sub=c("MDD ~ SNP"), parallel = FALSE, smooth_check = FALSE, toler = 1e-60)

write.csv(results2, file = paste0(output2), row.names = FALSE, quote = TRUE)