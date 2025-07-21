######################################
##### Filter Genes based on Chi2 #####
######################################

setwd("/mnt/lustre/working/lab_esked/damianWo/Chapter2/12_TWAS/T_SEM")

library(data.table) 
library(dplyr)
library(tidyr)

## Combine Files ### 

## Name directory and file pattern 
directory_path <- "/mnt/lustre/working/lab_esked/damianWo/Chapter2/12_TWAS/T_SEM/2_metab"

# List all files matching the pattern
all_files <- list.files(path = directory_path, pattern = "*TSEM$", full.names = TRUE)

# Select common pathway files only 
common_files <- all_files[!grepl("independent", all_files)]

# Function to standardise warning columns 
standardize_types <- function(df) {
  # Convert specific columns to character
  df$warning <- as.character(df$warning)  
  return(df)
}

# Read in files, standardise warnings, and combine 
common_genes <- common_files %>%
  lapply(read.csv) %>%
  lapply(standardize_types) %>%
  bind_rows()

## Create tissue column 
common_genes$tissue <- str_extract(common_genes$Panel, "(?<=GTExv8\\.EUR\\.)[^/]+")

## Remove tissue that are irrelevant 
common_genes <- common_genes %>% filter(!tissue %in% c("Kidney_Cortex", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Testis"))

## Check removal
unique(common_genes$tissue)

# Select independent files 
independent_files <- all_files[grepl("independent", all_files)]

# Read in files, standardise warnings, and combine 
independent_genes <- independent_files %>%
  lapply(read.csv) %>%
  lapply(standardize_types) %>%
  bind_rows()

##############################################################################################

#### Data Checking and Cleaning ####

## Check for errors 
subset_errors <- subset(common_genes, common_genes$error != 0)

## Check for warnings 
warningGENE <- common_genes %>% filter(warning != 0) # Subset warning SNPS 
num <- unique(warningGENE$warning) # Identify number of unique warnings 

# Define the specific warning message to be replaced
warning_message <- "lavaan WARNING:\n    the optimizer (NLMINB) claimed the model converged, but not all\n    elements of the gradient are (near) zero; the optimizer may not\n    have found a local solution use check.gradient = FALSE to skip\n    this check."

# Replace the specific warning message with 0
common_genes$warning <- ifelse(common_genes$warning == warning_message, 0, common_genes$warning)
warningGENE <- common_genes %>% filter(warning != 0) ## Check warnings
Gene <- warningGENE %>% select(Gene) ## Select SNPs if needed 

## Filter file to no warning 
common_genes <- common_genes %>% filter(warning == 0)

## Check for errors in independent pathway
subset_errors <- subset(independent_genes, independent_genes$error != 0)

##############################################################################################

### Calculating Chi-square for each GENE  ###

## Rename chisq columns and independent pathways fo chisquare pvalue
common_sub <- common_genes %>% 
  rename(
    chisq_common = chisq,
    chisq_df_common = chisq_df,
    chisq_pval_common = chisq_pval
  )

independent_sub <- independent_genes %>% 
  select(Gene, Panel, chisq, chisq_df, chisq_pval) %>%
  rename(
    chisq_independent = chisq,
    chisq_df_independent = chisq_df,
    chisq_pval_independent = chisq_pval
  )

## Join independent pathway chi square 
chisquare_df <- left_join(common_sub, independent_sub, by = c("Gene","Panel"))

## Create Q_GENE columns
Q_df <- chisquare_df %>% 
  mutate(Q_chisq = chisq_common - chisq_independent) %>%
  mutate(df = chisq_df_common - chisq_df_independent) %>% 
  mutate(Q_chisq_pvalue = pchisq(Q_chisq, df, lower.tail = FALSE))

##############################################################################################

## Write output path of QGENE removed GWAS 
output_file_path <- "/mnt/lustre/working/lab_esked/damianWo/Chapter2/12_TWAS/7_genes/Metab_combined.genes.csv"

## Save as .tsv table 
write.table(Q_df, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
