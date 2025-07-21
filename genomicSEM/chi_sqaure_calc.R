############################################################
##### Filter SNPs based on Chi2 and in LD with Chi Sq ######
############################################################

setwd("./2_multivariateGWAS/")

library(data.table) 
library(dplyr)
library(tidyr)

## Load in file for common and independent pathways ## 
common_file <- fread("./GSEM_CVD_common_combined.tsv.gz")
independent_file <- fread("./GSEM_CVD__independent_combined.tsv.gz")

##############################################################################################

#### Calculate Effective Sample Size of Factor #####

## Estimate Factor Expected Sample Size 
subset_df <- subset(common_file, common_file$MAF <= .4 & common_file$MAF >= .1 & !is.na(SE))

# To calculate N_Hat remove any NA or infinite values 
na_maf <- any(is.na(subset_df$MAF)) # Check NA
na_se <- any(is.na(subset_df$SE)) # Check NA
inf_se <- any(is.infinite(subset_df$SE)) # Check infinite values 

N_hat <- mean(1/((2*subset_df$MAF*(1-subset_df$MAF))*subset_df$SE^2))

##############################################################################################

#### Data Checking and Cleaning ####

## Check for errors 
subset_errors <- subset(common_file, common_file$error != 0)

## Check for warnings 
warningSNPS <- common_file %>% filter(warning != 0) # Subset warning SNPS 
num <- unique(warningSNPS$warning) # Identify number of unique warnings 

# Define the specific warning message to be replaced
warning_message <- "lavaan WARNING:\n    the optimizer (NLMINB) claimed the model converged, but not all\n    elements of the gradient are (near) zero; the optimizer may not\n    have found a local solution use check.gradient = FALSE to skip\n    this check."

# Replace the specific warning message with 0
common_file$warning <- ifelse(common_file$warning == warning_message, 0, common_file$warning)
warningSNPS <- common_file %>% filter(warning != 0) ## Check warnings
SNPs <- warningSNPS %>% select(SNP) ## Select SNPs if needed 

## Filter file to no warning 
common_file <- common_file %>% filter(warning == 0)

## Check for errors in independent pathway
subset_errors <- subset(independent_file, independent_file$error != 0)

## Save Common file No QSNP removal for comparison 

## Write output path 
output_file_path <- "./6_sumstats/MDD_CVD_Factor_NoQSNP.tsv"

## Remove unnecessary columns to reduce file size 
no_QSNP <- common_file[,1:15]

## Save as .tsv table 
write.table(no_QSNP, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

##############################################################################################

### Calculating Chi-square for each SNP in the GWAS ###

## Rename chisq columns and independent pathways fo chisquare pvalue
common_sub <- common_file %>% 
  rename(
    chisq_common = chisq,
    chisq_df_common = chisq_df,
    chisq_pval_common = chisq_pval
  )

independent_sub <- independent_file %>% 
  select(SNP, chisq, chisq_df, chisq_pval) %>%
  rename(
    chisq_independent = chisq,
    chisq_df_independent = chisq_df,
    chisq_pval_independent = chisq_pval
  )

## Join independent pathway chi square 
chisquare_df <- left_join(common_sub, independent_sub, by = "SNP")

## Create Q_snp columns
Q_df <- chisquare_df %>% 
  mutate(Q_chisq = chisq_common - chisq_independent) %>%
  mutate(df = chisq_df_common - chisq_df_independent) %>% 
  mutate(Q_chisq_pvalue = pchisq(Q_chisq, df, lower.tail = FALSE))

##############################################################################################

#### Removing QSNPs and SNPs in LD

## Rename columns
Q_df <- Q_df %>% 
  rename( 
    P = Q_chisq_pvalue,
    Pval_est = P)

### Create file path for clumping 
output_file_path <- "./cvd_clump/cvd_QSNP_Pfile.tsv"

## Save as .tsv table 
write.table(Q_df, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

#### PERFORM CLUMPING IN PLINK #####

## Read in clumping file 
clump <- fread("./cvd_clump/GSEM_cvd_common_clump.clumped")

## Count SNPS to check if matched with snps to remove ##

## Select columns 
clump <- clump %>% select(SNP, SP2)

## Remove (1) from columns 
clump$SP2 <- gsub("\\(1\\)", "", clump$SP2)

## Ensure SP2 column is split into individual SNPs
clump$SP2 <- strsplit(as.character(clump$SP2), split = ",")

## Create vector of unique SNPs to remove - the QSNPs and SNPs in LD 
split_snps <- lapply(clump$SP2, function(x) unlist(strsplit(as.character(x), ",")))
split_snps <- unlist(split_snps)
snps_to_remove <- c(clump$SNP, split_snps)
snps_to_remove <- unique(snps_to_remove)

## Remove QSNP SNPs and LD SNPs
GWAS <- common_file %>% filter(!SNP %in% snps_to_remove)

## Remove unnecessary columns to reduce file size 
GWAS <- GWAS[,1:15]

##############################################################################################

## Write output path of QSNP removed GWAS 
output_file_path <- "./MDD_CVD_Factor.tsv"

## Save as .tsv table 
write.table(GWAS, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)

## Count number of SNPs 
total_SNPs_no_removal <- nrow(common_file)
total_SNPs_removal <- nrow(GWAS)
genomewide_sig_no_removal <- common_file %>% filter(P < 5e-8) %>% nrow()
genomewide_sig_removal <- GWAS %>% filter(P < 5e-8) %>% nrow()