# Author: Rui Xiang Yu
# Date: 2024 November 5th
# This script creates the pseudobulk files from the 0.Raw files, as well as normalizes them. It also performs DEA.
# Usage: R/6-pseudobulkanddea.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)
library(edgeR)

main <- function() {
  ### create pseudobulks and cpm them. also cpm filtered one
  cohort_list <- c("CMC", "MultiomeBrain", "SZBDMulti-Seq")
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
    
    for (cell in cell_list) {
      full_path <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell)
      unfiltered <- readRDS(full_path)
      pseudobulked <- create_pseudo_bulk(unfiltered$expr, unfiltered$meta)
      final_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)
      saveRDS(pseudobulked, final_path)
      
      cpm_matrix <- do_cpm_log(unfiltered$expr, FALSE)
      final_cpmed <- list(expr = cpm_matrix, meta = unfiltered$meta)
      
      pseudobulk_cpm <- do_cpm_log(pseudobulked$expr, FALSE)
      final_cpmpb <- list(expr = pseudobulk_cpm, meta = pseudobulked$meta)
      
      cpm_path <- paste0("data/data_processed/", cohort, "/SingleCell/CPM/", cell)
      saveRDS(final_cpmed, cpm_path)
      
      pb_cpm_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkCPM/", cell)
      saveRDS(final_cpmpb, pb_cpm_path)
    }
  }

  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
    
    for (cell in cell_list) {
      full_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)
      pseudobulk <- readRDS(full_path)
      dea_res <- perform_DGE(pseudobulk)
      
      dea_path <- paste0("results/DEA/", cohort, "/TMMDEAresults_", cell)
      saveRDS(dea_res$dge_results, dea_path)
      
      tmm_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkTMM/", cell)
      saveRDS(dea_res$tmmfile, tmm_path)
      }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
    
    for (cell in cell_list) {
      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkCPM/", cell)
      pseudobulk <- readRDS(pb_path)
      dea_res <- perform_DGE_on_CPM(pseudobulk)
      
      dea_path <- paste0("results/DEA/", cohort, "/CPMDEAresults_", cell)
      saveRDS(dea_res$dge_results, dea_path)
    }
  }
  
}

create_pseudo_bulk <- function(expr, meta) {
 
   # https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
  
  mm <- model.matrix(~ 0 + patientID:disorder, data = meta)
  mat_mm <- expr %*% mm
  
  mat_mm = as.matrix(mat_mm)
  #keep_genes <- rowSums(mat_mm > 0) > 0
  #mat_mm <- mat_mm[keep_genes, ] %>% as.matrix() %>% as.data.frame()
  mat_mm <- mat_mm |> as.matrix() |> as.data.frame()
  
  # Clean up column names
  colnames(mat_mm) <- gsub("replicate|label", "", colnames(mat_mm))
  
  # Drop empty columns
  keep_samples <- colSums(mat_mm) > 0
  mat_mm <- mat_mm[, keep_samples]
  mat_mm <- mat_mm[, ]
  
  # Create the targets data frame
  targets <- data.frame(group_sample = colnames(mat_mm)) %>%
    mutate(group = gsub(".*\\:", "", group_sample),
           patientID = sub(".*patientID(.*)\\:.*", "\\1", group_sample))
  
  targets <- targets %>%
    left_join(meta %>% select(patientID, disorder, sex, age) %>% distinct(),
              by = "patientID")
  
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
  
  #save tmm matrix
  tmm_normalized_matrix <- y$counts
  final_tmmfile <- list(expr = tmm_normalized_matrix, meta = PB$meta)
  
  return(list(dge_results = res, tmmfile = final_tmmfile))
}

do_cpm_log <- function(mtx, log) {
  colsums <- colSums(mtx)
  cpm_result <- t(t(mtx) / colsums * 1e6)
  
  if (log) {
    cpm_result <- log1p(cpm_result)
  }
  
  return(cpm_result)
}

perform_DGE_on_CPM <- function(PB) {
  # Check if all sample names in the expression matrix match group sample names
  if (!all(colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }
  
  design <- model.matrix(~ group, data = PB$meta)
  y <- DGEList(counts = PB$expr, group = PB$meta$group)
  y$samples$norm.factors <- 1
 
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

main()