# Author: Rui Xiang Yu
# Date: 2024 November 22nd
# This script creates visualizations from the DEA results.
# Usage: R/7-deavisualizations.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)
library(edgeR)
library(EnhancedVolcano)
library(limma)

main <- function(){
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
      res_path <- paste0("results/DEA/", cohort, "/TMMWithCells/DEAresults_", cell)
      
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
    
    cohort_title <- paste0(cohort, ": TMM with Age, Sex and Cell Number")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCells/pvalueplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCells/pvalueplot_", cohort, ".png")
      
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
      res_path <- paste0("results/DEA/", cohort, "/TMMWithCells/DEAresults_", cell)
      
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
    
    cohort_title <- paste0(cohort, ": TMM with Age, Sex, and Cell Number")
    
    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
    
    title_plot <- ggdraw() + 
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    
    if (cohort == "Batiuk") {
      
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCells/logFCplot_", cohort, ".png")
      
      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
      
    }
    
    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
      
      final_path <- paste0("results/DEA/", cohort, "/TMMWithCells/logFCplot_", cohort, ".png")
      
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

main()