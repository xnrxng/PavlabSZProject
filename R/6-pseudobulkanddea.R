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
  
  ### create p-value hists
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPM/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
    results <- readRDS(res_path)
    
    res_his <- results |>
      ggplot(aes(x = PValue)) +
      geom_histogram() +
      theme_classic()+
      labs(
        title = cell_type,
        x = "P-value",
        y = NULL
      )
    
    plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPM")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPMLog/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = PValue)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPMLog")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLog/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLog/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/TMM/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = PValue)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": TMM")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMM/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMM/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = PValue)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPMLog with Age and Sex")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = PValue)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": TMM with Age and Sex")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  ###create hists for logFC
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPM/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPM")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPMLog/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPMLog")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLog/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLog/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/TMM/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": TMM")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMM/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMM/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": CPMLog with Age and Sex")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPMLogWithCovariates/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )
      
      plot_list[[cell_type]] <- res_his
    }
    
    cohort_title <- paste0(cohort, ": TMM with Age and Sex")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCovariates/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
  ### create enhanced volcano plots
  for (cohort in cohort_list){
    astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
    excitatory <- paste0("Exc_", cohort, "_SZ.rds")
    inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
    microglia <- paste0("Mic_", cohort, "_SZ.rds")
    oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
    opc <- paste0("Opc_", cohort, "_SZ.rds")
    gli <- paste0("Gli_", cohort, "_SZ.rds")
    
    cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc, gli)
    plot_list <- list()
    
    for (cell in cell_list) {
      res_path <- paste0("results/DEA/", cohort, "/CPM/DEAresults_", cell)
      
      cell_type <- strsplit(cell, "_")[[1]][1]
      
      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)
      
      res_volcano <- EnhancedVolcano(results, lab = results$gene, 
                      x = 'logFC', y = 'FDR', title = cell_type, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                      axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$FDR))+1))
      
      plot_list[[cell_type]] <- res_volcano
    }
    
    cohort_title <- paste0(cohort, ": CPM")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/volcanoplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/CPM/volcanoplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
  }
  
}

### helper functions

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
  res <- topTags(test, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
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
  res <- topTags(test, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    mutate(test = 'pseudobulk_edgeR')
  
  return(res)
}

main()