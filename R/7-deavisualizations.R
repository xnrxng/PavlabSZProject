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
  cohort_list <- c("Batiuk", "CMC", "SZBDMulti-Seq", "MultiomeBrain")

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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = P.Value)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )

      plot_list[[cell_type]] <- res_his
    }

    cohort_title <- paste0(cohort, ": CPMLog (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/pvalueplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/pvalueplot_", cohort, ".png")

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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = P.Value)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )

      plot_list[[cell_type]] <- res_his
    }

    cohort_title <- paste0(cohort, ": CPMLog with Age and Sex (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/pvalueplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/pvalueplot_", cohort, ".png")

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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = P.Value)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "P-value",
          y = NULL
        )

      plot_list[[cell_type]] <- res_his
    }

    cohort_title <- paste0(cohort, ": CPMLog with Age, Sex, and Cell number (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/pvalueplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/pvalueplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }
  }

  exc_neurons <- c("L2_3_CUX2_FREM3", "L2_CUX2_LAMP5", "L5_6_THEMIS", "L5_6_FEZF2_TLE4", "L3_CUX2_PRSS12",
                   "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A")

  inh_neurons <- c("ID2_LAMP5", "VIP", "ID2_PAX6", "ID2_NCKAP5", "PVALB", "SST")

  method_list <- c("LimmaCPMLog", "LimmaCPMLogAgeSex")
  for (method in method_list){
    plot_list <- list()
    for (excitatory in exc_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", excitatory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = P.Value)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = excitatory,
          x = "P-value",
          y = NULL
        )

      plot_list[[excitatory]] <- res_his}

      cohort_title <- paste0("Batiuk: ", get_title(method))

      all_plots <- plot_grid(plotlist = plot_list, ncol = 4, align = "hv")

      title_plot <- ggdraw() +
        draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
        all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
        final_path <- paste0("results/DEA/Batiuk/", method, "/excitatorypvalueplot_Batiuk.png")
        ggsave(plot = all_plots_with_title, filename = final_path, width = 12)
  }

  for (method in method_list){
    plot_list <- list()
    for (inhibitory in inh_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", inhibitory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = P.Value)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = inhibitory,
          x = "P-value",
          y = NULL
        )

      plot_list[[inhibitory]] <- res_his}

    cohort_title <- paste0("Batiuk: ", get_title(method))

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
    final_path <- paste0("results/DEA/Batiuk/", method, "/inhibitorypvalueplot_Batiuk.png")
    cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/DEAresults_", cell)

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

    cohort_title <- paste0(cohort, ": CPMLog (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/logFCplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLog/logFCplot_", cohort, ".png")

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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = groupdisorderyes)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )

      plot_list[[cell_type]] <- res_his
    }

    cohort_title <- paste0(cohort, ": CPMLog with Age and Sex (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/logFCplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/logFCplot_", cohort, ".png")

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
      res_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      res_his <- results |>
        ggplot(aes(x = groupdisorderyes)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )

      plot_list[[cell_type]] <- res_his
    }

    cohort_title <- paste0(cohort, ": CPMLog with Age, Sex, and Cell Number (Limma)")

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/logFCplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSexCell/logFCplot_", cohort, ".png")

      cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

    }
  }

  exc_neurons <- c("L2_3_CUX2_FREM3", "L2_CUX2_LAMP5", "L5_6_THEMIS", "L5_6_FEZF2_TLE4", "L3_CUX2_PRSS12",
                   "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A")

  inh_neurons <- c("ID2_LAMP5", "VIP", "ID2_PAX6", "ID2_NCKAP5", "PVALB", "SST")

  method_list <- c("LimmaCPMLog", "LimmaCPMLogAgeSex")
  for (method in method_list){
    plot_list <- list()
    for (excitatory in exc_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", excitatory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)


      if (method == "LimmaCPMLog"){

      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = excitatory,
          x = "logFC",
          y = NULL
        )
      }

      else{
        res_his <- results |>
          ggplot(aes(x = groupdisorderyes)) +
          geom_histogram() +
          theme_classic()+
          labs(
            title = excitatory,
            x = "logFC",
            y = NULL
          )
      }
      plot_list[[excitatory]] <- res_his}

    cohort_title <- paste0("Batiuk: ", get_title(method))

    all_plots <- plot_grid(plotlist = plot_list, ncol = 4, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
    final_path <- paste0("results/DEA/Batiuk/", method, "/excitatorylogFCplot_Batiuk.png")
    ggsave(all_plots_with_title, filename = final_path, width = 12)
  }

  for (method in method_list){
    plot_list <- list()
    for (inhibitory in inh_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", inhibitory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)

      if (method == "LimmaCPMLog"){
      res_his <- results |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = inhibitory,
          x = "logFC",
          y = NULL
        )}

      else{
        res_his <- results |>
          ggplot(aes(x = groupdisorderyes)) +
          geom_histogram() +
          theme_classic()+
          labs(
            title = inhibitory,
            x = "logFC",
            y = NULL
          )
      }

      plot_list[[inhibitory]] <- res_his}

    cohort_title <- paste0("Batiuk: ", get_title(method))

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
    final_path <- paste0("results/DEA/Batiuk/", method, "/inhibitorylogFCplot_Batiuk.png")
    cowplot::save_plot(plot = all_plots_with_title, filename = final_path)
  }


  ### create enhanced volcano plots

  method_list <- c("LimmaCPMLog", "LimmaCPMLogAgeSex", "LimmaCPMLogAgeSexCell")

  for (method in method_list){
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
      res_path <- paste0("results/DEA/", cohort, "/", method, "/DEAresults_", cell)

      cell_type <- strsplit(cell, "_")[[1]][1]

      if (!file.exists(res_path)) {
        message("File ", res_path, " does not exist. Skipping to next.")
        next
      }
      results <- readRDS(res_path)

      if(method == "LimmaCPMLog"){

      res_volcano <- EnhancedVolcano(results, lab = rownames(results),
                                     x = 'logFC', y = 'adj.P.Val', title = cell_type, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                                     axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))

      plot_list[[cell_type]] <- res_volcano
      }
      else{
        res_volcano <- EnhancedVolcano(results, lab = rownames(results),
                                       x = 'groupdisorderyes', y = 'adj.P.Val', title = cell_type, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                                       axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))

        plot_list[[cell_type]] <- res_volcano
      }
    }

    cohort_title <- paste0(cohort, ": ", get_title(method))

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

    if (cohort == "Batiuk") {

      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/", method, "/volcanoplot_", cohort, ".png")
      ggsave(all_plots_with_title, filename = final_path, width = 12)

    }

    else{
      all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

      final_path <- paste0("results/DEA/", cohort, "/", method, "/volcanoplot_", cohort, ".png")
      ggsave(all_plots_with_title, filename = final_path, width = 12, height = 9)
    }
  }}

  non_limma <- c("CPM", "CPMLog", "CPMLogWithCovariates", "TMM", "TMMWithCells", "TMMWithCovariates")

  for (method in non_limma){
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
        res_path <- paste0("results/DEA/", cohort, "/", method, "/DEAresults_", cell)

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

      cohort_title <- paste0(cohort, ": ", get_title(method))

      all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

      title_plot <- ggdraw() +
        draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)

      if (cohort == "Batiuk") {

        all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

        final_path <- paste0("results/DEA/", cohort, "/", method, "/volcanoplot_", cohort, ".png")
        ggsave(all_plots_with_title, filename = final_path, width = 12)

      }

      else{
        all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

        final_path <- paste0("results/DEA/", cohort, "/", method, "/volcanoplot_", cohort, ".png")
        ggsave(all_plots_with_title, filename = final_path, width = 12, height = 9)
      }
    }}
  ### batiuk subtypes

  method_list <- c("LimmaCPMLog", "LimmaCPMLogAgeSex")

  for (method in method_list){
    plot_list <- list()
    for (excitatory in exc_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", excitatory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)

      if(method == "LimmaCPMLog"){

        res_volcano <- EnhancedVolcano(results, lab = rownames(results),
            x = 'logFC', y = 'adj.P.Val', title = excitatory, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
            axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))

        plot_list[[excitatory]] <- res_volcano
      }

        else{
          res_volcano <- EnhancedVolcano(results, lab = rownames(results),
            x = 'groupdisorderyes', y = 'adj.P.Val', title = excitatory, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
            axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))
          plot_list[[excitatory]] <- res_volcano
        }
      }

      cohort_title <- paste0("Batiuk: ", get_title(method))

      all_plots <- plot_grid(plotlist = plot_list, ncol = 2, align = "hv")

      title_plot <- ggdraw() +
        draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
     all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
     final_path <- paste0("results/DEA/Batiuk/", method, "/excitatoryvolcanoplot_Batiuk.png")
    ggsave(all_plots_with_title, filename = final_path, height = 12)
  }

  for (method in method_list){
    plot_list <- list()
    for (inhibitory in inh_neurons){
      res_path <- paste0("results/DEA/Batiuk/", method, "/DEAresults_", inhibitory, "_Batiuk_SZ.rds")
      results <- readRDS(res_path)

      if(method == "LimmaCPMLog"){

        res_volcano <- EnhancedVolcano(results, lab = rownames(results),
                                       x = 'logFC', y = 'adj.P.Val', title = inhibitory, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                                       axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))

        plot_list[[inhibitory]] <- res_volcano
      }

      else{
        res_volcano <- EnhancedVolcano(results, lab = rownames(results),
                                       x = 'groupdisorderyes', y = 'adj.P.Val', title = inhibitory, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                                       axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1, ylim = c(0, max(-log10(results$adj.P.Val))+1))
        plot_list[[inhibitory]] <- res_volcano
      }
    }

    cohort_title <- paste0("Batiuk: ", get_title(method))

    all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

    title_plot <- ggdraw() +
      draw_label(cohort_title, x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))
    final_path <- paste0("results/DEA/Batiuk/", method, "/inhibitoryvolcanoplot_Batiuk.png")
    ggsave(all_plots_with_title, filename = final_path, width = 12, height = 9)
  }
}
### helper functions
get_title <- function(method){
  if (method == "LimmaCPMLog"){
    return("CPMLog (Limma)")
  }

  if (method == "LimmaCPMLogAgeSex"){
    return("CPMLog with Age and Sex (Limma)")
  }

  if (method == "LimmaCPMLogAgeSexCell"){
    return("CPMLog with Age, Sex, and Cell Number (Limma)")
  }
  if (method == "CPM"){
    return("CPM")
  }
  if (method == "CPMLog"){
    return("CPMLog")
  }
  if (method == "CPMLogWithCovariates"){
    return("CPMLog with Age and Sex")
  }
  if (method == "TMM"){
    return("TMM")
  }
  if (method == "TMMWithCells"){
    return("TMM with Age, Sex, and Cell Number")
  }
  if (method == "TMMWithCovariates"){
    return("TMM with Age and Sex")
  }
}

main()
