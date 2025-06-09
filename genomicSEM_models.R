### Genomic SEM Models for MDD and Disease Groups ### 

library(devtools)
require(GenomicSEM)

## Load in LDSC results 
load("LDSCoutput_CVD.RData")

## CVD-MDD Single Factor Model ## 

model_CVD <- '
F_CVD=~ NA*HF + CAD + AF + STK

F_CVD ~~ 1*F_CVD

MDD ~ F_CVD
'

model_CVD <-usermodel(LDSCoutput, estimation = "DWLS", model = modelCVD, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)


## Metabolic-MDD Single Factor Model ## 

## Load in LDSC results 
load("LDSCoutput_metab.RData")

modelmetab <- '
F_metab =~ NA*T2D + b*MSYN + TG + HDL 
F_metab ~~ 1*F_metab

MSYN ~~ a*MSYN
a > 0.001
b < 1

MDD ~ F_metab
'

model_metab <-usermodel(LDSCoutput, estimation = "DWLS", model = modelmetab, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)


## Gastro-MDD Single Factor Model ## 

## Load in LDSC results 
load("LDSCoutput_gastro.RData")

modelgastro <- '
F_gastro =~ NA*GORD + GSD + IBS + PUD 

F_gastro ~~ 1*F_gastro

MDD ~ F_gastro
'

model_gastro <-usermodel(LDSCoutput, estimation = "DWLS", model = modelgastro, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)


## Gastro-MDD Single Factor Model ## 

## Load in LDSC results 
load("LDSCoutput_IM.RData")

modelIM <- '
F_IM =~ NA*AD + ASTH + IBD + MS 
F_IM ~~ 1*F_IM 

ASTH ~~ AD 

MDD ~ F_IM
'

model_IM <-usermodel(LDSCoutput, estimation = "DWLS", model = modelIM, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)


## Four Factor Model with MDD ## 

## Load in LDSC results 
load("LDSCoutput_all.traits.RData")

modelall <- '
F_CVD =~ NA*CAD + HF + STK + AF 
F_metab =~ NA*T2D + b*MSYN + TG + HDL + GSD 
F_gastro =~ NA*GORD + GSD + IBS + PUD 
F_IM =~ NA*ASTH + AD + IBD + MS

F_metab ~~ 1*F_metab
F_CVD ~~ 1*F_CVD
F_gastro ~~ 1*F_gastro
F_IM ~~ 1*F_IM

F_CVD ~~ F_metab 
F_CVD ~~ F_gastro
F_CVD ~~ F_IM
F_metab ~~ F_gastro
F_metab ~~ F_IM
F_gastro ~~ F_IM

ASTH ~~ AD

MSYN ~~ a*MSYN
a > 0.001
b < 1

MDD ~ F_CVD 
MDD ~ F_metab
MDD ~ F_gastro 
MDD ~ F_IM 
'

model_all <- usermodel(LDSCoutput, estimation = "DWLS", model = modelall, CFIcalc = TRUE, std.lv = FALSE, imp_cov = TRUE)

