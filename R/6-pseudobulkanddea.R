# Author: Rui Xiang Yu
# Date: 2024 November 5th
# This script creates the pseudobulk files from the 0.RaW files, as well as normalizes them. It also performs DEA.
# Usage: R/6-pseudobulkanddea.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)

main <- function() {
  
  
  
}


create_pseudo_bulk <- function(expr, meta) {
  
  library(dplyr)
  # https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
  
  mm <- model.matrix(~ 0 + patientID:disorder, data = meta)
  mat_mm <- expr %*% mm
  
  mat_mm = as.matrix(mat_mm)
  keep_genes <- rowSums(mat_mm > 0) > 0
  mat_mm <- mat_mm[keep_genes, ] %>% as.matrix() %>% as.data.frame()
  
  # Clean up column names
  colnames(mat_mm) <- gsub("replicate|label", "", colnames(mat_mm))
  
  # Drop empty columns
  keep_samples <- colSums(mat_mm) > 0
  mat_mm <- mat_mm[, keep_samples]
  mat_mm <- mat_mm[, ]
  
  # Create the targets data frame
  targets <- data.frame(group_sample = colnames(mat_mm)) %>%
    mutate(group = gsub(".*\\:", "", group_sample))
  
  # Adjust the group factor
  targets$group <- as.factor(targets$group)
  targets$group <- relevel(targets$group, ref = "disorderno")
  
  # Create the PB list
  PB <- list(meta = targets, expr = mat_mm)
  
  return(PB)
}

perform_DGE <- function(PB) {
  # Check if all sample names in the expression matrix match group sample names
  if (!all(colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }
  
  # Load necessary libraries
  library(edgeR)
  library(dplyr)
  
  ### make it return tmm normalized matrix
  design <- model.matrix(~ group, data = PB$meta)
  y <- DGEList(counts = PB$expr, group = PB$meta$group)
  y <- calcNormFactors(y, method = "TMM")  # Normalization
  x <- estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)   # Estimate dispersion
  fit <- glmFit(x, design = design)   # Fit GLM (edgeR)
  test <- glmLRT(fit)   # Likelihood ratio test
  
  # Extract results
  res <- topTags(test, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    mutate(test = 'pseudobulk_edgeR')
  
  return(res)
}

do_cpm_log <- function(mtx, log) {
  colsums <- colSums(mtx)
  cpm_result <- t(t(mtx) / colsums * 1e6)
  
  if (log) {
    cpm_result <- log1p(cpm_result)
  }
  
  return(cpm_result)
}