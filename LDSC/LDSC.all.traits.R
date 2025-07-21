#### Multivairable LDSC #####

library(devtools)
require(GenomicSEM)
library(data.table)
library(dplyr)

## Set working directory
setwd("/mnt/lustre/working/lab_esked/damianWo/Chapter2/1_genomicSEM")

## import all sumstat files 
traits <- list.files(pattern = "*sumstats.gz$", recursive = TRUE)

## Enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev <- c(0.5,
                 0.5,
                 0.5,
                 0.5,
                 0.5,
                 NA,
                 0.5,
                 0.5,
                 NA,
                 0.5,
                 0.5,
                 0.5,
                 0.5,
                 0.5,
                 0.5,
                 0.5,
                 0.5)

## Vector of population prevalences
population.prev <- c(
  0.075,
  0.015,
  0.035,
  0.015,
  0.015,
  NA,
  0.24,
  0.07,
  NA,
  0.15,
  0.125,
  0.125,
  0.075,
  0.1,
  0.08,
  0.008130,
  0.001)

## The folder of LD scores
ld <- "/mnt/lustre/working/lab_esked/damianWo/reference_files/eur_w_ld_chr/"

## The folder of LD weights
wld <- "/mnt/lustre/working/lab_esked/damianWo/reference_files/eur_w_ld_chr/"

## Create vector for trait names 
trait.names <- sub(".sumstats\\.gz$", "", basename(traits))

## Run LDSC
LDSCoutput <- ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

save(LDSCoutput,file="LDSCoutput_all.traits5.RData")
