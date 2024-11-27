# Author: Rui Xiang Yu
# Date: 2024 November 5th
# This script creates the pseudobulk files from the 0.Raw files, as well as normalizes them. It also performs DEA.
# Usage: R/6-pseudobulkanddea.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)
library(edgeR)
library(EnhancedVolcano)
library(limma)

main <- function() {
  ### create pseudobulks and cpm them. also cpm filtered one
  cohort_list <- c("CMC", "MultiomeBrain", "SZBDMulti-Seq", "Batiuk")

  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      full_path <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell)

      if (!file.exists(full_path)) {
        message("File ", full_path, " does not exist. Skipping to next.")
        next
      }
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

  ### perform DEA and save TMM PB

  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      full_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)

      if (!file.exists(full_path)) {
        message("File ", full_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(full_path)
      dea_res <- perform_DGE(pseudobulk)

      dea_path <- paste0("results/DEA/", cohort, "/TMM/DEAresults_", cell)
      saveRDS(dea_res$dge_results, dea_path)

      tmm_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkTMM/", cell)
      saveRDS(dea_res$tmmfile, tmm_path)
      }
  }

  ### perform DEA on CPM PB
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkCPM/", cell)

      if (!file.exists(pb_path)) {
        message("File ", pb_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(pb_path)
      dea_res <- perform_DGE_on_CPM(pseudobulk)

      dea_path <- paste0("results/DEA/", cohort, "/CPM/DEAresults_", cell)
      saveRDS(dea_res, dea_path)
    }
  }

  ### CPMLOG PBs and perform DEA
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      full_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)

      if (!file.exists(full_path)) {
        message("File ", full_path, " does not exist. Skipping to next.")
        next
      }
      unfiltered <- readRDS(full_path)

      cpm_matrix <- do_cpm_log(unfiltered$expr, TRUE)
      final_cpmed <- list(expr = cpm_matrix, meta = unfiltered$meta)

      cpm_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkCPMLog/", cell)
      saveRDS(final_cpmed, cpm_path)

      results <- perform_DGE_on_CPM(final_cpmed)

      cpmlogres_path <- paste0("results/DEA/", cohort, "/CPMLog/DEAresults_", cell)
      saveRDS(results, cpmlogres_path)
    }
  }

  ### DEA on CPMLog with other covariates
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkCPMLog/", cell)

      if (!file.exists(pb_path)) {
        message("File ", pb_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(pb_path)
      dea_res <- perform_DGE_with_covs(pseudobulk)

      dea_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/DEAresults_", cell)
      saveRDS(dea_res, dea_path)
    }
  }

  ### DEA on TMM with other covariates
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)

      if (!file.exists(pb_path)) {
        message("File ", pb_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(pb_path)
      dea_res <- perform_DGE_with_covs_on_tmm(pseudobulk)

      dea_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/DEAresults_", cell)
      saveRDS(dea_res, dea_path)
    }
  }

  ### including on TMM single cell number as a covariate
  celltype_summary_CMC <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )

  celltype_summary_SZBD <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )

  celltype_summary_MB <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )

  celltype_summary_Batiuk <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )

  celltype_summary_CMC <- countcells_pertype("CMC", celltype_summary_CMC)
  celltype_summary_SZBD <- countcells_pertype("SZBDMulti-Seq", celltype_summary_SZBD)
  celltype_summary_MB <- countcells_pertype("MultiomeBrain", celltype_summary_MB)
  celltype_summary_Batiuk <- countcells_pertype("Batiuk", celltype_summary_Batiuk)

  all_celltype <- rbind(celltype_summary_Batiuk, celltype_summary_CMC, celltype_summary_MB, celltype_summary_SZBD)

  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      celltype <- strsplit(cell, "_")[[1]][1]

      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)

      if (!file.exists(pb_path)) {
        message("File ", pb_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(pb_path)

      filt_celltype <- filter(all_celltype, cell_type == celltype)

      pseudobulk$meta <- merge(pseudobulk$meta, filt_celltype[, c("patientID", "cells")],
                               by = "patientID", all.x = TRUE)

      dea_res <- perform_DGE_with_cells(pseudobulk)

      dea_path <- paste0("results/DEA/", cohort, "/TMMWithCells/DEAresults_", cell)
      saveRDS(dea_res, dea_path)
    }
  }

  ### limma cpmlog
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")

    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

    for (cell in cell_list) {
      celltype <- strsplit(cell, "_")[[1]][1]

      pb_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", cell)

      if (!file.exists(pb_path)) {
        message("File ", pb_path, " does not exist. Skipping to next.")
        next
      }

      pseudobulk <- readRDS(pb_path)
      filt_celltype <- filter(all_celltype, cell_type == celltype)
      pseudobulk$meta <- merge(pseudobulk$meta, filt_celltype[, c("patientID", "cells")],
                               by = "patientID", all.x = TRUE)
      pseudobulk$meta$age <- gsub("\\+", "", pseudobulk$meta$age)
      pseudobulk$meta$age <- as.numeric(pseudobulk$meta$age)
      pseudobulk$meta$sex <- as.factor(pseudobulk$meta$sex)
      pseudobulk$meta$cells <- as.numeric(pseudobulk$meta$cells)

      design = model.matrix(~ group, data = pseudobulk$meta)
      dea_res <- limma_dge(PB = pseudobulk, design = design)
      dea_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/DEAresults_", cell)
      saveRDS(dea_res, dea_path)

      design_agesex = model.matrix(~ group + age + sex, data = pseudobulk$meta)
      dea_res_age_sex <- limma_dge(PB = pseudobulk, design = design_agesex)
      dea_path_age_sex <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/DEAresults_", cell)
      saveRDS(dea_res_age_sex, dea_path_age_sex)

      design_agesexcell = model.matrix(~ group + age + sex + cells, data = pseudobulk$meta)
      dea_res_age_sex_cell <- limma_dge(PB = pseudobulk, design = design_agesexcell)
      dea_path_age_sex_cell <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/DEAresults_", cell)
      saveRDS(dea_res_age_sex_cell, dea_path_age_sex_cell)
    }
  }

  ### pseudobulk batiuk subtypes
  exc_neurons <- c("L2_3_CUX2_FREM3", "L2_CUX2_LAMP5", "L5_6_THEMIS", "L5_6_FEZF2_TLE4", "L3_CUX2_PRSS12",
                   "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A")

  inh_neurons <- c("ID2_LAMP5", "VIP", "ID2_PAX6", "ID2_NCKAP5", "PVALB", "SST")

  for (excitatory in exc_neurons){
    initial_path <- paste0("data/data_processed/Batiuk/Excitatory/FilteredV1/", excitatory, "_Batiuk_SZ.rds")
    list_file <- readRDS(initial_path)

    pseudobulked <- create_pseudo_bulk(expr = list_file$expr, meta = list_file$meta)
    pb_path <- paste0("data/data_processed/Batiuk/Excitatory/PseudobulkRaw/", excitatory, "_Batiuk_SZ.rds")
    saveRDS(pseudobulked, pb_path)

    pseudobulked$meta$age <- as.numeric(pseudobulked$meta$age)
    pseudobulked$meta$sex <- as.factor(pseudobulked$meta$sex)

    design = model.matrix(~ group, data = pseudobulked$meta)
    regular <- limma_dge(design, pseudobulked)
    reg_path <- paste0("results/DEA/Batiuk/LimmaCPMLog/DEAresults_", excitatory, "_Batiuk_SZ.rds")
    saveRDS(regular, reg_path)

    design_covs = model.matrix(~ group + age + sex, data = pseudobulked$meta)
    withcovariates <- limma_dge(design_covs, pseudobulked)
    cov_path <- paste0("results/DEA/Batiuk/LimmaCPMLogAgeSex/DEAresults_", excitatory, "_Batiuk_SZ.rds")
    saveRDS(withcovariates, cov_path)
  }

  for (inhibitory in inh_neurons){
    initial_path <- paste0("data/data_processed/Batiuk/Inhibitory/FilteredV1/", inhibitory, "_Batiuk_SZ.rds")
    list_file <- readRDS(initial_path)

    pseudobulked <- create_pseudo_bulk(expr = list_file$expr, meta = list_file$meta)
    pb_path <- paste0("data/data_processed/Batiuk/Inhibitory/PseudobulkRaw/", inhibitory, "_Batiuk_SZ.rds")
    saveRDS(pseudobulked, pb_path)

    pseudobulked$meta$age <- as.numeric(pseudobulked$meta$age)
    pseudobulked$meta$sex <- as.factor(pseudobulked$meta$sex)

    design = model.matrix(~ group, data = pseudobulked$meta)
    regular <- limma_dge(design, pseudobulked)
    reg_path <- paste0("results/DEA/Batiuk/LimmaCPMLog/DEAresults_", inhibitory, "_Batiuk_SZ.rds")
    saveRDS(regular, reg_path)

    design_covs = model.matrix(~ group + age + sex, data = pseudobulked$meta)
    withcovariates <- limma_dge(design_covs, pseudobulked)
    cov_path <- paste0("results/DEA/Batiuk/LimmaCPMLogAgeSex/DEAresults_", inhibitory, "_Batiuk_SZ.rds")
    saveRDS(withcovariates, cov_path)
  }
}

### helper functions

create_pseudo_bulk <- function(expr, meta) {

   # https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/

  mm <- model.matrix(~ 0 + patientID:disorder, data = meta)
  mat_mm <- expr %*% mm

  mat_mm = as.matrix(mat_mm)
  keep_genes <- rowSums(mat_mm > 0) > 0
  mat_mm <- mat_mm[keep_genes, ] |> as.matrix() |> as.data.frame()
  mat_mm <- mat_mm |> as.matrix() |> as.data.frame()

  # Clean up column names
  colnames(mat_mm) <- gsub("replicate|label", "", colnames(mat_mm))

  # Drop empty columns
  keep_samples <- colSums(mat_mm) > 0
  mat_mm <- mat_mm[, keep_samples]
  mat_mm <- mat_mm[, ]

  # Create the targets data frame
  targets <- data.frame(group_sample = colnames(mat_mm))|>
    mutate(group = gsub(".*\\:", "", group_sample),
           patientID = sub(".*patientID(.*)\\:.*", "\\1", group_sample))

  targets <- targets |>
    left_join(meta |> dplyr::select(patientID, disorder, sex, age) |> distinct(),
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
  res <- topTags(test, n = Inf) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
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

  x <- estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)   # Estimate dispersion
  fit <- glmFit(x, design = design)   # Fit GLM (edgeR)
  test <- glmLRT(fit)   # Likelihood ratio test

  # Extract results
  res <- topTags(test, n = Inf) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
    mutate(test = 'pseudobulk_edgeR')

  return(res)
}

perform_DGE_with_covs <- function(PB) {
  # Check if all sample names in the expression matrix match group sample names
  if (!all(colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }

  PB$meta$age <- gsub("\\+", "", PB$meta$age)
  PB$meta$age <- as.numeric(PB$meta$age)
  PB$meta$sex <- as.factor(PB$meta$sex)

  design <- model.matrix(~ group + sex + age, data = PB$meta)
  y <- DGEList(counts = PB$expr, group = PB$meta$group)

  x <- estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)   # Estimate dispersion
  fit <- glmFit(x, design = design)   # Fit GLM (edgeR)
  test <- glmLRT(fit)   # Likelihood ratio test

  # Extract results
  res <- topTags(test, n = Inf) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
    mutate(test = 'pseudobulk_edgeR')

  return(res)
}

perform_DGE_with_covs_on_tmm <- function(PB) {
  # Check if all sample names in the expression matrix match group sample names
  if (!all(colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }

  PB$meta$age <- gsub("\\+", "", PB$meta$age)
  PB$meta$age <- as.numeric(PB$meta$age)
  PB$meta$sex <- as.factor(PB$meta$sex)

  ### make it return tmm normalized matrix
  design <- model.matrix(~ group + sex + age, data = PB$meta)
  y <- DGEList(counts = PB$expr, group = PB$meta$group)
  y <- calcNormFactors(y, method = "TMM")  # Normalization
  x <- estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)   # Estimate dispersion
  fit <- glmFit(x, design = design)   # Fit GLM (edgeR)
  test <- glmLRT(fit)   # Likelihood ratio test

  # Extract results
  res <- topTags(test, n = Inf) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
    mutate(test = 'pseudobulk_edgeR')

  return(res)
}

countcells_pertype <- function(cohort, results_df) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")
  gli <- paste0("Gli_", cohort, "_SZ.rds")

  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)

  for (cell_type_file in cell_list) {
    cell_type <- strsplit(cell_type_file, "_")[[1]][1]

    sample_file <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell_type_file)

    if (!file.exists(sample_file)) {
      message("File ", sample_file, " does not exist. Skipping to next.")
      next
    }

    sample_data <- readRDS(sample_file)

    n_cells <- sample_data$meta |>
      group_by(patientID, disorder) |>
      summarize(cells = n()) |>
      mutate(cell_type = cell_type) |>
      dplyr::select(patientID, cells, disorder, cell_type)
    results_df <- bind_rows(results_df, n_cells)
  }

  return(results_df)
}

perform_DGE_with_cells <- function(PB) {
  # Check if all sample names in the expression matrix match group sample names
  if (!all(colnames(PB$expr) %in% PB$meta$group_sample)) {
    stop("Sample names in expression matrix do not match group sample names in meta data")
  }

  PB$meta$age <- gsub("\\+", "", PB$meta$age)
  PB$meta$age <- as.numeric(PB$meta$age)
  PB$meta$sex <- as.factor(PB$meta$sex)
  PB$meta$cells <- as.numeric(PB$meta$cells)

  ### make it return tmm normalized matrix
  design <- model.matrix(~ group + sex + age + cells, data = PB$meta)
  y <- DGEList(counts = PB$expr, group = PB$meta$group)
  y <- calcNormFactors(y, method = "TMM")  # Normalization
  x <- estimateDisp(y, design, trend.method = "locfit", tagwise = TRUE, prior.df = NULL)   # Estimate dispersion
  fit <- glmFit(x, design = design)   # Fit GLM (edgeR)
  test <- glmLRT(fit)   # Likelihood ratio test

  # Extract results
  res <- topTags(test, n = Inf) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
    mutate(test = 'pseudobulk_edgeR')

  return(res)
}

limma_dge <- function(design, PB){
  x = voom(as.matrix(PB$expr), design)
  fit = lmFit(x, design) |> eBayes()
  res = fit |> topTable(number = Inf)
  return(res)
}


main()
