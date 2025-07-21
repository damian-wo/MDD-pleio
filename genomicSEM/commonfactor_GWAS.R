########################################################
###### Run Common Factor GWAS on Disease Groups ########
########################################################

#########
## CVD ##
#########

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1

## Run multivariate GWAS for CVD only factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

setwd("./3_common_factor_models/1_cvd/cvd_factor/")

## Load in LDSC results and sumstats results
load("./3_common_factor_models/1_cvd/LDSCoutput_CVD.RData")
load("./3_common_factor_models/1_cvd/CVD_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

model1 <- "
F_CVD =~ 1*HF + CAD + AF + STK
F_CVD ~~ b*F_CVD

HF ~~ a*HF
a > 0.001

b < 1

F_CVD ~ SNP
"

results  <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, estimation = "DWLS", sub = c("F_CVD ~ SNP"), parallel = FALSE, Q_SNP = TRUE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

###########
## Metab ##
###########

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1

## Run multivariate GWAS for Metab only factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

setwd("./3_common_factor_models/2_metab/metab_factor/")

## Load in LDSC results and sumstats results
load("./3_common_factor_models/2_metab/LDSCoutput_metab.RData")
load("./3_common_factor_models/2_metab/metab_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

model1 <- "
F_metab =~ MSYN + TG + HDL + GSD + T2D
F_metab ~~ b*F_metab

MSYN ~~ a*MSYN
a > 0.001

b < 1

F_metab ~ SNP
"

results  <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, sub = c("F_metab ~ SNP"), estimation = "DWLS", parallel = FALSE, Q_SNP = TRUE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

############
## Gastro ##
############

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1

## Run multivariate GWAS for CVD only factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

setwd("./3_common_factor_models/3_gastro/gastro_factor")

## Load in LDSC results and sumstats results
load("./3_common_factor_models/3_gastro/LDSCoutput_gastro.RData")
load("./3_common_factor_models/3_gastro/gastro_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]


model1 <- "
F_gastro =~ GORD + GSD + IBS + PUD 

F_gastro ~~ NA*F_gastro

F_gastro ~ SNP
"

results  <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, sub = c("F_gastro ~ SNP"), estimation = "DWLS", parallel = FALSE, Q_SNP = TRUE)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)

############
## Immune ##
############

### Set arguments from pbs script ###

args = commandArgs(trailingOnly=TRUE)
n_start <- args[1] #Nstart
n_stop <- args[2] #Nend
output1 <- args[3] #Output 1

## Run multivariate GWAS for CVD only factor ##

library(devtools)
library(data.table)
require(GenomicSEM)

setwd("./3_common_factor_models/4_im/im_factor/")

## Load in LDSC results and sumstats results
load("./3_common_factor_models/4_im/LDSCoutput_IM.RData")
load("./3_common_factor_models/4_im/IM_linearsumstats.RData")
subsetSNPs <- linear_sumstats[paste0(n_start):paste0(n_stop),]

model1 <- "
F_IM =~ IBD + AD + ASTH + MS 
F_IM ~~ NA*F_IM

AD ~~ ASTH

F_IM ~ SNP
"

results  <- userGWAS(covstruc = LDSCoutput, SNPs = subsetSNPs, model = model1, sub = c("F_IM ~ SNP"), estimation = "DWLS", parallel = FALSE, Q_SNP = TRUE, toler = 1e-60)

write.csv(results, file = paste0(output1), row.names = FALSE, quote = TRUE)



