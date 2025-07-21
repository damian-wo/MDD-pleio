### Conduct ACAT on TSEM Genes ###

library(devtools)
library(dplyr)
library(data.table)
devtools::install_github("yaowuliu/ACAT")

setwd("./ACAT")

##################################################

CVD <- fread("CVD_pvalue_matrix.tsv")

## Change rownames to gene names 
CVD <- as.data.frame(CVD)
rownames(CVD) <- CVD$Gene
CVD$Gene <- NULL

# Perform ACAT test for each gene(row) in the df 
acat_results <- apply(CVD, 1, function(row) {
  
  # Remove NA values
  pvals <- row[!is.na(row)]
  
  # Run ACAT only if there are valid P-values
  if (length(pvals) > 0) {
    return(list(
      p_value = ACAT::ACAT(pvals),  # ACAT p-value
      n_tissues = length(pvals)    # Number of non-NA tissues
    ))
  } else {
    return(list(
      p_value = NA,
      n_tissues = 0
    ))
  }
})

# Extract the results into separate vectors
acat_pvals <- sapply(acat_results, function(res) res$p_value)
n_tissues <- sapply(acat_results, function(res) res$n_tissues)

# Create dataframe 
data <- data.frame(
  Gene = rownames(CVD),      # Gene names from the input matrix
  ACAT_pval = acat_pvals,    # Extracted ACAT p-values
  n_tissues = n_tissues      # Extracted tissue counts
)

# Apply Multiple Testing Correction 
data <- data %>% filter(ACAT_pval < 0.05/nrow(data))

# Save file 
fwrite(data, "CVD_ACAT.tsv", sep = "\t") 

##############################################################

Metab <- fread("Metab_pvalue_matrix.tsv")

## Change rownames to gene names 
Metab <- as.data.frame(Metab)
rownames(Metab) <- Metab$Gene
Metab$Gene <- NULL

# Perform ACAT test for each gene(row) in the df 
acat_results <- apply(Metab, 1, function(row) {
  
  # Remove NA values
  pvals <- row[!is.na(row)]
  
  # Run ACAT only if there are valid P-values
  if (length(pvals) > 0) {
    return(list(
      p_value = ACAT::ACAT(pvals),  # ACAT p-value
      n_tissues = length(pvals)    # Number of non-NA tissues
    ))
  } else {
    return(list(
      p_value = NA,
      n_tissues = 0
    ))
  }
})

# Extract the results into separate vectors
acat_pvals <- sapply(acat_results, function(res) res$p_value)
n_tissues <- sapply(acat_results, function(res) res$n_tissues)

# Create dataframe 
data <- data.frame(
  Gene = rownames(Metab),      # Gene names from the input matrix
  ACAT_pval = acat_pvals,    # Extracted ACAT p-values
  n_tissues = n_tissues      # Extracted tissue counts
)

# Apply Multiple Testing Correction 
data <- data %>% filter(ACAT_pval < 0.05/nrow(data))

# Save file 
fwrite(data, "Metab_ACAT.tsv", sep = "\t") 

#############################################

Gastro <- fread("Gastro_pvalue_matrix.tsv")

## Change rownames to gene names 
Gastro <- as.data.frame(Gastro)
rownames(Gastro) <- Gastro$Gene
Gastro$Gene <- NULL

# Perform ACAT test for each gene(row) in the df 
acat_results <- apply(Gastro, 1, function(row) {
  
  # Remove NA values
  pvals <- row[!is.na(row)]
  
  # Run ACAT only if there are valid P-values
  if (length(pvals) > 0) {
    return(list(
      p_value = ACAT::ACAT(pvals),  # ACAT p-value
      n_tissues = length(pvals)    # Number of non-NA tissues
    ))
  } else {
    return(list(
      p_value = NA,
      n_tissues = 0
    ))
  }
})

# Extract the results into separate vectors
acat_pvals <- sapply(acat_results, function(res) res$p_value)
n_tissues <- sapply(acat_results, function(res) res$n_tissues)

# Create dataframe 
data <- data.frame(
  Gene = rownames(Gastro),      # Gene names from the input matrix
  ACAT_pval = acat_pvals,    # Extracted ACAT p-values
  n_tissues = n_tissues      # Extracted tissue counts
)

# Apply Multiple Testing Correction 
data <- data %>% filter(ACAT_pval < 0.05/nrow(data))

# Save file 
fwrite(data, "Gastro_ACAT.tsv", sep = "\t") 

###################################################

Immune <- fread("Immune_pvalue_matrix.tsv")

## Change rownames to gene names 
Immune <- as.data.frame(Immune)
rownames(Immune) <- Immune$Gene
Immune$Gene <- NULL

# Perform ACAT test for each gene(row) in the df 
acat_results <- apply(Immune, 1, function(row) {
  
  # Remove NA values
  pvals <- row[!is.na(row)]
  
  # Run ACAT only if there are valid P-values
  if (length(pvals) > 0) {
    return(list(
      p_value = ACAT::ACAT(pvals),  # ACAT p-value
      n_tissues = length(pvals)    # Number of non-NA tissues
    ))
  } else {
    return(list(
      p_value = NA,
      n_tissues = 0
    ))
  }
})

# Extract the results into separate vectors
acat_pvals <- sapply(acat_results, function(res) res$p_value)
n_tissues <- sapply(acat_results, function(res) res$n_tissues)

# Create dataframe 
data <- data.frame(
  Gene = rownames(Immune),      # Gene names from the input matrix
  ACAT_pval = acat_pvals,    # Extracted ACAT p-values
  n_tissues = n_tissues      # Extracted tissue counts
)

# Apply Multiple Testing Correction 
data <- data %>% filter(ACAT_pval < 0.05/nrow(data))

# Save file 
fwrite(data, "Immune_ACAT.tsv", sep = "\t") 

###################################################

MDD <- fread("MDD_pvalue_matrix.tsv")

## Change rownames to gene names 
MDD <- as.data.frame(MDD)
rownames(MDD) <- MDD$Gene
MDD$Gene <- NULL

# Perform ACAT test for each gene(row) in the df 
acat_results <- apply(MDD, 1, function(row) {
  
  # Remove NA values and values of 1 
  pvals <- row[!is.na(row) & row < 1]
  
  # Run ACAT only if there are valid P-values
  if (length(pvals) > 0) {
    return(list(
      p_value = ACAT::ACAT(pvals),  # ACAT p-value
      n_tissues = length(pvals)    # Number of non-NA tissues
    ))
  } else {
    return(list(
      p_value = NA,
      n_tissues = 0
    ))
  }
})

# Extract the results into separate vectors
acat_pvals <- sapply(acat_results, function(res) res$p_value)
n_tissues <- sapply(acat_results, function(res) res$n_tissues)

# Create dataframe 
data <- data.frame(
  Gene = rownames(MDD),      # Gene names from the input matrix
  ACAT_pval = acat_pvals,    # Extracted ACAT p-values
  n_tissues = n_tissues      # Extracted tissue counts
)

# Apply Multiple Testing Correction 
data <- data %>% filter(ACAT_pval < 0.05/nrow(data))

# Save file 
fwrite(data, "MDD_ACAT.tsv", sep = "\t") 

