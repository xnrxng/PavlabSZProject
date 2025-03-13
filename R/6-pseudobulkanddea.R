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

  ### ling
  ling_meta <- read_csv("data/data_processed/Ling/all-patient-metaData.csv")

  cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc")

  for (celltype in cell_types){
    initial_path <- paste0("data/data_processed/Ling/PseudobulkRaw/", celltype, "_ling-2024_meta.csv")
    cell_meta <- read_csv(initial_path)

    merged_meta <- left_join(cell_meta, ling_meta, by = "patientID") |>
      dplyr::select(patientID, sex, age, diagnosis)

    final_path <- paste0("data/data_processed/Ling/PseudobulkRaw/", celltype, "_Ling_MetaFiltered.rds")
    saveRDS(merged_meta, final_path)
  }
    for (celltype in cell_types){
      meta_path <- paste0("data/data_processed/Ling/PseudobulkRaw/", celltype, "_Ling_MetaFiltered.rds")
      cell_meta <- readRDS(meta_path) |> as.data.frame()
      rownames(cell_meta) <- cell_meta$patientID
      cell_meta$diagnosis <- as.factor(cell_meta$diagnosis)
      cell_meta$diagnosis <- relevel(cell_meta$diagnosis, ref = "no")
      cell_meta$age <- as.numeric(cell_meta$age)
      cell_meta$sex <- as.factor(cell_meta$sex)

      cell_type_path <- paste0("data/data_processed/CMC/FilteredV1/", celltype, "_CMC_SZ.rds")
      cell_type_genes <- readRDS(cell_type_path)
      cell_type_genes_CMC <- rownames(cell_type_genes$expr)

      expr_path <- paste0("data/data_processed/Ling/PseudobulkRaw/", celltype, "_ling-2024_expr.csv")
      cell_expr <- read_csv(expr_path) |> as.data.frame()
      rownames(cell_expr) <- cell_expr$gene_names
      cell_expr <- cell_expr[, -1]
      cell_expr <- cell_expr[rownames(cell_expr) %in% cell_type_genes_CMC, ]
      cell_expr <- as.matrix(cell_expr)
      cell_expr <- Matrix(cell_expr, sparse = TRUE)

      PB <- list(expr = cell_expr, meta = cell_meta)

      design <- model.matrix(~ diagnosis, data = PB$meta)
      results <- limma_dge(design = design, PB = PB)
      cpmlogpath <- paste0("results/DEA/Ling/LimmaCPMLog/DEAresults_", celltype, "_Ling_SZ.rds")
      saveRDS(results, cpmlogpath)

      designagesex <- model.matrix(~ diagnosis+age+sex, data = PB$meta)
      resultsagesex <- limma_dge(design = designagesex, PB = PB)
      agesexpath <- paste0("results/DEA/Ling/LimmaCPMLogAgeSex/DEAresults_", celltype, "_Ling_SZ.rds")
      saveRDS(resultsagesex, agesexpath)
    }

  ### nairuz microglia
  all_meta <- fread("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
  studies <- unique(all_meta$study)
  studies <- studies[studies != "ROSMAP-microGlia"]

  for (study in studies){
    meta_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/Mic_", study, "_meta.csv")
    cell_meta <- read_csv(meta_path) |> as.data.frame()  |> na.omit()
    rownames(cell_meta) <- cell_meta$patientID
    cell_meta$diagnosis <- as.factor(cell_meta$diagnosis)
    cell_meta$diagnosis <- relevel(cell_meta$diagnosis, ref = "no")
    #cell_meta$age <- as.numeric(cell_meta$age)
    #cell_meta$sex <- as.factor(cell_meta$sex)

    expr_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/Mic_", study, "_expr.csv")
    cell_expr <- read_csv(expr_path) |> as.data.frame()
    rownames(cell_expr) <- cell_expr[,1]
    cell_expr <- cell_expr[, -1]
    cell_expr <- as.matrix(cell_expr)
    cell_expr <- Matrix(cell_expr, sparse = TRUE)

    if (length(colnames(cell_expr)) > length(rownames(cell_meta))) {
      cell_expr <- cell_expr[, colnames(cell_expr) %in% rownames(cell_meta)]
    }

    if (length(colnames(cell_expr)) < length(rownames(cell_meta))) {
      cell_meta <- cell_meta[rownames(cell_meta) %in% colnames(cell_expr), ]
    }

    PB <- list(expr = cell_expr, meta = cell_meta)

    design <- model.matrix(~ diagnosis, data = PB$meta)
    results <- limma_dge(design = design, PB = PB)
    micpath <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/nairuz-microglia/DEA/DEAres_Mic_", study, ".rds")
    saveRDS(results, micpath)
  }

  all_res <- list()
  for (study in studies){
    micpath <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/nairuz-microglia/DEA/DEAres_Mic_", study, ".rds")
    deares <- readRDS(micpath)
    all_res[[study]] <- deares
  }

  saveRDS(all_res, "data/data_processed/nairuz-microglia/DEA/all_res.rds")

  ### DGE for SCZ
  all_meta <- fread("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
  studies <- unique(all_meta$study)
  studies <- studies[studies != "ROSMAP-microGlia"]
  cell_types <- c("Ast", "Mic", "Opc", "Oli", "Exc", "Inh")

  for (study in studies){
    for (cell_type in cell_types){
    meta_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_meta.csv")

    if (!file.exists(meta_path)) {
      message("Skipping: Meta file not found - ", meta_path)
      next
    }

    cell_meta <- read_csv(meta_path) |> as.data.frame()  |> na.omit()
    cell_meta <- cell_meta |> filter(disease == "CTL" | disease == "SCZ")
    rownames(cell_meta) <- cell_meta$patientID
    cell_meta$diagnosis <- as.factor(cell_meta$diagnosis)
    cell_meta$diagnosis <- relevel(cell_meta$diagnosis, ref = "no")
    #cell_meta$age <- as.numeric(cell_meta$age)
    #cell_meta$sex <- as.factor(cell_meta$sex)

    if (nlevels(cell_meta$diagnosis) == 1) next

    expr_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_expr.csv")

    if (!file.exists(expr_path)) {
      message("Skipping: Expression file not found - ", expr_path)
      next
    }

    cell_expr <- read_csv(expr_path) |> as.data.frame()
    rownames(cell_expr) <- cell_expr[,1]
    cell_expr <- cell_expr[, -1]
    cell_expr <- as.matrix(cell_expr)
    cell_expr <- Matrix(cell_expr, sparse = TRUE)

    if (length(colnames(cell_expr)) > length(rownames(cell_meta))) {
      cell_expr <- cell_expr[, colnames(cell_expr) %in% rownames(cell_meta)]
    }

    if (length(colnames(cell_expr)) < length(rownames(cell_meta))) {
      cell_meta <- cell_meta[rownames(cell_meta) %in% colnames(cell_expr), ]
    }

    cell_meta <- cell_meta[colnames(cell_expr), , drop = FALSE]

    PB <- list(expr = cell_expr, meta = cell_meta)

    design <- model.matrix(~ diagnosis, data = PB$meta)
    results <- limma_dge(design = design, PB = PB)
    respath <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/results/DEA/SCZ_18/", cell_type, "_", study, ".rds")
    saveRDS(results, respath)
    gc()
  }}

  all_res <- list()
  for (study in studies){
    for (cell_type in cell_types){
    sczpath <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/results/DEA/SCZ_18/", cell_type, "_", study, ".rds")

    if (!file.exists(sczpath)) {
      message("Skipping: File not found - ", sczpath)
      next
    }

    deares <- readRDS(sczpath)
    if (!is.list(all_res[[study]])) {
      all_res[[study]] <- list()
    }

    all_res[[study]][[cell_type]] <- deares
  }}

  saveRDS(all_res, "results/DEA/SCZ_18/all_res.rds")

  ### DGE for AD
  all_meta <- fread("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
  studies <- unique(all_meta$study)
  studies <- studies[studies != "ROSMAP-microGlia"]
  cell_types <- c("Ast", "Mic", "Opc", "Oli", "Exc", "Inh")

  for (study in studies){
    for (cell_type in cell_types){
      meta_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_meta.csv")

      if (!file.exists(meta_path)) {
        message("Skipping: Meta file not found - ", meta_path)
        next
      }

      cell_meta <- read_csv(meta_path) |> as.data.frame()  |> na.omit()
      cell_meta <- cell_meta |> filter(disease == "CTL" | disease == "AD")
      rownames(cell_meta) <- cell_meta$patientID
      cell_meta$diagnosis <- as.factor(cell_meta$diagnosis)
      cell_meta$diagnosis <- relevel(cell_meta$diagnosis, ref = "no")
      #cell_meta$age <- as.numeric(cell_meta$age)
      #cell_meta$sex <- as.factor(cell_meta$sex)

      if (nlevels(cell_meta$diagnosis) == 1) next

      expr_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_expr.csv")

      if (!file.exists(expr_path)) {
        message("Skipping: Expression file not found - ", expr_path)
        next
      }

      cell_expr <- read_csv(expr_path) |> as.data.frame()
      rownames(cell_expr) <- cell_expr[,1]
      cell_expr <- cell_expr[, -1]
      cell_expr <- as.matrix(cell_expr)
      cell_expr <- Matrix(cell_expr, sparse = TRUE)

      if (length(colnames(cell_expr)) > length(rownames(cell_meta))) {
        cell_expr <- cell_expr[, colnames(cell_expr) %in% rownames(cell_meta)]
      }

      if (length(colnames(cell_expr)) < length(rownames(cell_meta))) {
        cell_meta <- cell_meta[rownames(cell_meta) %in% colnames(cell_expr), ]
      }

      cell_meta <- cell_meta[colnames(cell_expr), , drop = FALSE]

      PB <- list(expr = cell_expr, meta = cell_meta)

      design <- model.matrix(~ diagnosis, data = PB$meta)
      results <- limma_dge(design = design, PB = PB)
      respath <- paste0("results/DEA/AD_18/", cell_type, "_", study, ".rds")
      saveRDS(results, respath)
      gc()
    }}

  all_res <- list()
  for (study in studies){
    for (cell_type in cell_types){
      sczpath <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/results/DEA/AD_18/", cell_type, "_", study, ".rds")

      if (!file.exists(sczpath)) {
        message("Skipping: File not found - ", sczpath)
        next
      }

      deares <- readRDS(sczpath)
      if (!is.list(all_res[[study]])) {
        all_res[[study]] <- list()
      }

      all_res[[study]][[cell_type]] <- deares
    }}

  saveRDS(all_res, "results/DEA/AD_18/all_res.rds")

  ### save as lists
  all_meta <- fread("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
  studies <- unique(all_meta$study)
  studies <- studies[studies != "ROSMAP-microGlia"]
  cell_types <- c("Ast", "Mic", "Opc", "Oli", "Exc", "Inh")

  for (study in studies){
    for (cell_type in cell_types){
      meta_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_meta.csv")

      if (!file.exists(meta_path)) {
        message("Skipping: Meta file not found - ", meta_path)
        next
      }

      cell_meta <- read_csv(meta_path) |> as.data.frame()  |> na.omit()
      rownames(cell_meta) <- cell_meta$patientID
      cell_meta$diagnosis <- as.factor(cell_meta$diagnosis)
      cell_meta$disease <- as.factor(cell_meta$disease)
      cell_meta$diagnosis <- relevel(cell_meta$disease, ref = "CTL")
      #cell_meta$age <- as.numeric(cell_meta$age)
      #cell_meta$sex <- as.factor(cell_meta$sex)

      expr_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_expr.csv")

      if (!file.exists(expr_path)) {
        message("Skipping: Expression file not found - ", expr_path)
        next
      }

      cell_expr <- read_csv(expr_path) |> as.data.frame()
      rownames(cell_expr) <- cell_expr[,1]
      cell_expr <- cell_expr[, -1]
      cell_expr <- as.matrix(cell_expr)
      cell_expr <- Matrix(cell_expr, sparse = TRUE)

      if (length(colnames(cell_expr)) > length(rownames(cell_meta))) {
        cell_expr <- cell_expr[, colnames(cell_expr) %in% rownames(cell_meta)]
      }

      if (length(colnames(cell_expr)) < length(rownames(cell_meta))) {
        cell_meta <- cell_meta[rownames(cell_meta) %in% colnames(cell_expr), ]
      }

      cell_meta <- cell_meta[colnames(cell_expr), , drop = FALSE]

      cpmed <- do_cpm_log(cell_expr)

      PB <- list(expr = cpmed, meta = cell_meta)
      cpmpath <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/cpm/", cell_type, "_", study, ".rds")
      saveRDS(PB, cpmpath)
      gc()

      cpmloged <- do_cpm_log(cell_expr, log = TRUE)

      PB_log <- list(expr = cpmloged, meta = cell_meta)
      cpmlogpath <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/cpm-log/", cell_type, "_", study, ".rds")
      saveRDS(PB_log, cpmlogpath)
      gc()
    }}

### summmaries
  all_meta <- fread("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
  studies <- unique(all_meta$study)
  studies <- studies[studies != "ROSMAP-microGlia"]
  cell_types <- c("Ast", "Mic", "Opc", "Oli", "Exc", "Inh")

  for (cell_type in cell_types){
    summary <- data.frame()

    for (study in studies){
      meta_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_meta.csv")

      if (!file.exists(meta_path)) {
        message("Skipping: Meta file not found - ", meta_path)
        next
      }

      cell_meta <- read_csv(meta_path) |> as.data.frame()  |> na.omit()
      rownames(cell_meta) <- cell_meta$patientID
      n_ad <- sum(cell_meta$disease == "AD")
      n_ctl <- sum(cell_meta$disease == "CTL")

      summary <- rbind(summary, data.frame(
        N_AD = n_ad,
        N_CTL = n_ctl,
        study = study,
        celltype = cell_type
      ))
    }
    write_tsv(summary, paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/summaries/patient_AD_", cell_type, ".tsv"))
  }

  for (cell_type in cell_types){
    summary <- data.frame()

    for (study in studies){
      expr_path <- paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/", cell_type, "_", study, "_expr.csv")

      if (!file.exists(expr_path)) {
        message("Skipping: Expression file not found - ", expr_path)
        next
      }

      cell_expr <- read_csv(expr_path) |> as.data.frame()
      rownames(cell_expr) <- cell_expr[,1]
      cell_expr <- cell_expr[, -1]
      cell_expr <- as.matrix(cell_expr)
      cell_expr <- Matrix(cell_expr, sparse = TRUE)

      summary <- rbind(summary, data.frame(
        N_genes = nrow(cell_expr),
        study = study,
        celltype = cell_type
      ))
    }
    write_tsv(summary, paste0("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/summaries/n_genes_", cell_type, ".tsv"))
  }


  ### subtypes
  excitatory_cell_types <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 CT",
                             "L6 IT Car3", "L5 ET", "L5/6 NP", "L6b")

  inhibitory_cell_types <- c("Sst", "Sst Chodl", "Pvalb", "Chandelier", "Pax6",
                             "Lamp5 Lhx6", "Lamp5", "Sncg", "Vip")

  all_subtypes <-c(excitatory_cell_types, inhibitory_cell_types)
  subtypes_safe <- gsub("[ /]", "_", all_subtypes)
  for (cohort in cohort_list){
    for (subtype in subtypes_safe){
      full_path <- paste0("data/data_processed/", cohort, "/FilteredV1/", subtype, "_", cohort, "_SZ.rds")

      if (!file.exists(full_path)) {
        message("File ", full_path, " does not exist. Skipping to next.")
        next
      }
      unfiltered <- readRDS(full_path)
      pseudobulked <- create_pseudo_bulk(unfiltered$expr, unfiltered$meta)
      final_path <- paste0("data/data_processed/", cohort, "/Pseudobulk/PseudobulkRaw/", subtype, "_", cohort, "_SZ.rds")
      saveRDS(pseudobulked, final_path)

      pseudobulked$meta$age <- gsub("\\+", "", pseudobulked$meta$age)
      pseudobulked$meta$age <- as.numeric(pseudobulked$meta$age)
      pseudobulked$meta$sex <- as.factor(pseudobulked$meta$sex)

      design_agesex = model.matrix(~ group + age + sex, data = pseudobulked$meta)
      dea_res_age_sex <- limma_dge(PB = pseudobulked, design = design_agesex)
      dea_path_age_sex <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/DEAresults_", subtype, "_", cohort, "_SZ.rds")
      saveRDS(dea_res_age_sex, dea_path_age_sex)
    }

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

do_cpm_log <- function(mtx, log = FALSE) {
  colsums <- Matrix::colSums(mtx)
  cpm_result <- Matrix::t(Matrix::t(mtx) / colsums * 1e6)

  if (log) {
    cpm_result <- base::log1p(cpm_result)
  }

  return(cpm_result)
}

main()
