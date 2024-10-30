# Author: Rui Xiang Yu
# Date: 2024 October 21st
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/4-szagedistributionandumap.R

library(tidyverse)
library(ggplot2)
library(cowplot)
library(Seurat)
library(future.apply)
library(Matrix)
library(magick)

main <- function() {
  ### CMC umap
  set.seed(11)

  CMC_ast <- readRDS("data/data_processed/CMC/Filtered/Filt_Ast_CMC_SZ.rds")
  CMC_ast$meta$cell_type <- "Ast"

  CMC_exc <- readRDS("data/data_processed/CMC/Filtered/Filt_Exc_CMC_SZ.rds")
  CMC_exc$meta$cell_type <- "Exc"

  CMC_inh <- readRDS("data/data_processed/CMC/Filtered/Filt_Inh_CMC_SZ.rds")
  CMC_inh$meta$cell_type <- "Inh"

  CMC_mic <- readRDS("data/data_processed/CMC/Filtered/Filt_Mic_CMC_SZ.rds")
  CMC_mic$meta$cell_type <- "Mic"

  CMC_oli <- readRDS("data/data_processed/CMC/Filtered/Filt_Oli_CMC_SZ.rds")
  CMC_oli$meta$cell_type <- "Oli"

  CMC_opc <- readRDS("data/data_processed/CMC/Filtered/Filt_Opc_CMC_SZ.rds")
  CMC_opc$meta$cell_type <- "Opc"

  CMCexpr_list <- list(CMC_ast$expr, CMC_opc$expr, CMC_exc$expr, CMC_inh$expr, CMC_mic$expr, CMC_oli$expr)

  CMCall_genes <- Reduce(union, lapply(CMCexpr_list, rownames))
  CMCpad_genes <- function(expr, CMCall_genes) {
    CMCmissing_genes <- setdiff(CMCall_genes, rownames(expr))
    CMCzero_mat <- matrix(0, nrow = length(CMCmissing_genes), ncol = ncol(expr))
    rownames(CMCzero_mat) <- CMCmissing_genes
    CMCexpr_padded <- rbind(expr, CMCzero_mat)
    return(CMCexpr_padded[CMCall_genes, , drop = FALSE])
  }

  CMCexpr_list_padded <- lapply(CMCexpr_list, CMCpad_genes, CMCall_genes = CMCall_genes)
  CMCcombined_expr <- do.call(cbind, CMCexpr_list_padded)

  CMCmetadata_list <- list(CMC_ast$meta, CMC_opc$meta, CMC_exc$meta, CMC_inh$meta, CMC_mic$meta, CMC_oli$meta)
  CMCcombined_meta <- do.call(rbind, CMCmetadata_list)

  CMC_seurat <- CreateSeuratObject(counts = CMCcombined_expr, meta.data = CMCcombined_meta)
  CMC_seurat <- NormalizeData(CMC_seurat)
  CMC_seurat <- FindVariableFeatures(CMC_seurat, selection.method = "vst", nfeatures = 2000)
  
  CMCall.genes <- rownames(CMC_seurat)
  CMC_seurat <- ScaleData(CMC_seurat, features = CMCall.genes)
  
  CMC_seurat <- RunPCA(CMC_seurat, features = VariableFeatures(object = CMC_seurat))
  
  CMC_seurat <- RunUMAP(CMC_seurat, dims = 1:17)
  CMC_umapplot <- DimPlot(CMC_seurat, reduction = "umap", group.by = "cell_type") +
    scale_color_manual(values = c("Ast" = "turquoise",
                                 "Opc" = "coral",
                                 "Exc" = "pink2",
                                 "Inh" = "olivedrab3",
                                 "Mic" = "navy",
                                 "Oli" = "yellow2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("CMC")
  
  CMC_disorderumapplot <- DimPlot(CMC_seurat, reduction = "umap", group.by = "disorder") +
    scale_color_manual(values = c("no" = "grey21", "yes" = "firebrick2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("CMC")
  
  CMC_sexumapplot <- DimPlot(CMC_seurat, reduction = "umap", group.by = "sex") +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("CMC")
  
  ggsave(file.path("results/9.1-CMC_umap.png"), CMC_umapplot)
  ggsave(file.path("results/9.2-CMC_umap_disorder.png"), CMC_disorderumapplot)
  ggsave(file.path("results/9.3-CMC_umap_sex.png"), CMC_sexumapplot)

  ### SZBDMulti-Seq umap
  SZBD_ast <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Ast_SZBDMulti-Seq_SZ.rds")
  SZBD_ast$meta$cell_type <- "Ast"

  SZBD_exc <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Exc_SZBDMulti-Seq_SZ.rds")
  SZBD_exc$meta$cell_type <- "Exc"

  SZBD_inh <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Inh_SZBDMulti-Seq_SZ.rds")
  SZBD_inh$meta$cell_type <- "Inh"

  SZBD_mic <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Mic_SZBDMulti-Seq_SZ.rds")
  SZBD_mic$meta$cell_type <- "Mic"

  SZBD_oli <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Oli_SZBDMulti-Seq_SZ.rds")
  SZBD_oli$meta$cell_type <- "Oli"

  SZBD_opc <- readRDS("data/data_processed/SZBDMulti-Seq/Filtered/Filt_Opc_SZBDMulti-Seq_SZ.rds")
  SZBD_opc$meta$cell_type <- "Opc"

  SZBDexpr_list <- list(SZBD_ast$expr, SZBD_opc$expr, SZBD_exc$expr, SZBD_inh$expr, SZBD_mic$expr, SZBD_oli$expr)

  SZBDall_genes <- Reduce(union, lapply(SZBDexpr_list, rownames))
  SZBDpad_genes <- function(expr, SZBDall_genes) {
    SZBDmissing_genes <- setdiff(SZBDall_genes, rownames(expr))
    SZBDzero_mat <- matrix(0, nrow = length(SZBDmissing_genes), ncol = ncol(expr))
    rownames(SZBDzero_mat) <- SZBDmissing_genes
    SZBDexpr_padded <- rbind(expr, SZBDzero_mat)
    return(SZBDexpr_padded[SZBDall_genes, , drop = FALSE])
  }

  SZBDexpr_list_padded <- lapply(SZBDexpr_list, SZBDpad_genes, SZBDall_genes = SZBDall_genes)
  SZBDcombined_expr <- do.call(cbind, SZBDexpr_list_padded)

  SZBDmetadata_list <- list(SZBD_ast$meta, SZBD_opc$meta, SZBD_exc$meta, SZBD_inh$meta, SZBD_mic$meta, SZBD_oli$meta)
  SZBDcombined_meta <- do.call(rbind, SZBDmetadata_list)

  SZBD_seurat <- CreateSeuratObject(counts = SZBDcombined_expr, meta.data = SZBDcombined_meta)
  SZBD_seurat <- NormalizeData(SZBD_seurat)
  SZBD_seurat <- FindVariableFeatures(SZBD_seurat, selection.method = "vst", nfeatures = 2000)
 
  SZBDall.genes <- rownames(SZBD_seurat)
  SZBD_seurat <- ScaleData(SZBD_seurat, features = SZBDall.genes)
  
  SZBD_seurat <- RunPCA(SZBD_seurat, features = VariableFeatures(object = SZBD_seurat))
  

  SZBD_seurat <- RunUMAP(SZBD_seurat, dims = 1:17)
  
  SZBD_umapplot <- DimPlot(SZBD_seurat, reduction = "umap", group.by = "cell_type") +
    scale_color_manual(values = c("Ast" = "turquoise",
                                  "Opc" = "coral",
                                  "Exc" = "pink2",
                                  "Inh" = "olivedrab3",
                                  "Mic" = "navy",
                                  "Oli" = "yellow2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("SZBD")
  
  SZBD_disorderumapplot <- DimPlot(SZBD_seurat, reduction = "umap", group.by = "disorder") +
    scale_color_manual(values = c("no" = "grey21", "yes" = "firebrick2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("SZBD")
  
  SZBD_sexumapplot <- DimPlot(SZBD_seurat, reduction = "umap", group.by = "sex") +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("SZBD")
  ggsave(file.path("results/9.4-SZBD_umap.png"), SZBD_umapplot)
  ggsave(file.path("results/9.5-SZBD_umap_disorder.png"), SZBD_disorderumapplot)
  ggsave(file.path("results/9.6-SZBD_umap_sex.png"), SZBD_sexumapplot)

  ### MultiomeBrain umap
  MB_ast <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Ast_MultiomeBrain_SZ.rds")
  MB_ast$meta$cell_type <- "Ast"

  MB_exc <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Exc_MultiomeBrain_SZ.rds")
  MB_exc$meta$cell_type <- "Exc"

  MB_inh <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Inh_MultiomeBrain_SZ.rds")
  MB_inh$meta$cell_type <- "Inh"

  MB_mic <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Mic_MultiomeBrain_SZ.rds")
  MB_mic$meta$cell_type <- "Mic"

  MB_oli <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Oli_MultiomeBrain_SZ.rds")
  MB_oli$meta$cell_type <- "Oli"

  MB_opc <- readRDS("data/data_processed/MultiomeBrain/Filtered/Filt_Opc_MultiomeBrain_SZ.rds")
  MB_opc$meta$cell_type <- "Opc"

  MBexpr_list <- list(MB_ast$expr, MB_opc$expr, MB_exc$expr, MB_inh$expr, MB_mic$expr, MB_oli$expr)

  MBall_genes <- Reduce(union, lapply(MBexpr_list, rownames))
  MBpad_genes <- function(expr, MBall_genes) {
    MBmissing_genes <- setdiff(MBall_genes, rownames(expr))
    MBzero_mat <- matrix(0, nrow = length(MBmissing_genes), ncol = ncol(expr))
    rownames(MBzero_mat) <- MBmissing_genes
    MBexpr_padded <- rbind(expr, MBzero_mat)
    return(MBexpr_padded[MBall_genes, , drop = FALSE])
  }

  MBexpr_list_padded <- lapply(MBexpr_list, MBpad_genes, MBall_genes = MBall_genes)
  MBcombined_expr <- do.call(cbind, MBexpr_list_padded)

  MBmetadata_list <- list(MB_ast$meta, MB_opc$meta, MB_exc$meta, MB_inh$meta, MB_mic$meta, MB_oli$meta)
  MBcombined_meta <- do.call(rbind, MBmetadata_list)

  MB_seurat <- CreateSeuratObject(counts = MBcombined_expr, meta.data = MBcombined_meta)
  MB_seurat <- NormalizeData(MB_seurat)
  MB_seurat <- FindVariableFeatures(MB_seurat, selection.method = "vst", nfeatures = 2000)

  MBall.genes <- rownames(MB_seurat)
  MB_seurat <- ScaleData(MB_seurat, features = MBall.genes)
  
  MB_seurat <- RunPCA(MB_seurat, features = VariableFeatures(object = MB_seurat))

  MB_seurat <- RunUMAP(MB_seurat, dims = 1:17)
  
  MB_umapplot <- DimPlot(MB_seurat, reduction = "umap", group.by = "cell_type") +
    scale_color_manual(values = c("Ast" = "turquoise",
                                  "Opc" = "coral",
                                  "Exc" = "pink2",
                                  "Inh" = "olivedrab3",
                                  "Mic" = "navy",
                                  "Oli" = "yellow2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("MB")
  
  MB_disorderumapplot <- DimPlot(MB_seurat, reduction = "umap", group.by = "disorder") +
    scale_color_manual(values = c("no" = "grey21", "yes" = "firebrick2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("MB")
  
  MB_sexumapplot <- DimPlot(MB_seurat, reduction = "umap", group.by = "sex") +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    xlab(NULL) +
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("MB")
  
  ggsave(file.path("results/9.7-MB_umap.png"), MB_umapplot)
  ggsave(file.path("results/9.8-MB_umap_disorder.png"), MB_disorderumapplot)
  ggsave(file.path("results/9.9-MB_umap_sex.png"), MB_sexumapplot)

  ### combine umaps in one plot
  CMC_cellumap <- ggdraw() + draw_image("results/9.1-CMC_umap.png")
  SZBD_cellumap <- ggdraw() + draw_image("results/9.4-SZBD_umap.png")
  MB_cellumap <- ggdraw() + draw_image("results/9.7-MB_umap.png")
  
  celltype_dummy <- data.frame(x = 1:6, y = 1:6, cell_type = c("Ast", "Opc", "Oli", "Exc", "Inh", "Mic"))
  celltype_dummyplot <- ggplot(celltype_dummy) +
    geom_point(aes(x = x, y = y, color = cell_type), size = 3) +
    labs(color = "Cell type") +
    scale_color_manual(values = c("Ast" = "turquoise",
                                  "Opc" = "coral",
                                  "Exc" = "pink2",
                                  "Inh" = "olivedrab3",
                                  "Mic" = "navy",
                                  "Oli" = "yellow2")) +
    theme_minimal()
  
  celltypelegend <- get_legend(celltype_dummyplot)
  cell_combined_plot <- plot_grid(CMC_cellumap, SZBD_cellumap, MB_cellumap, celltypelegend, ncol = 4)

  ggsave(file.path("results/9.10-cell_type_umap.png"), cell_combined_plot, dpi = 300, width = 8, height = 6)
  
  CMC_disorderumap <- ggdraw() + draw_image("results/9.2-CMC_umap_disorder.png")
  SZBD_disorderumap <- ggdraw() + draw_image("results/9.5-SZBD_umap_disorder.png")
  MB_disorderumap <- ggdraw() + draw_image("results/9.8-MB_umap_disorder.png")
  
  disorder_dummy <- data.frame(x = 1:2, y = 1:2, disorder = c("Schizophrenia", "Control"))
  disorder_dummyplot <- ggplot(disorder_dummy) +
    geom_point(aes(x = x, y = y, color = disorder), size = 3) +
    labs(color = "Disorder") +
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    theme_minimal()
  
  disorderlegend <- get_legend(disorder_dummyplot)
  disorder_combined_plot <- plot_grid(CMC_disorderumap, SZBD_disorderumap, MB_disorderumap, disorderlegend, ncol = 4)
  
  ggsave(file.path("results/9.11-disorder_umap.png"), disorder_combined_plot, dpi = 300, width = 8, height = 6)
  
  CMC_sexumap <- ggdraw() + draw_image("results/9.3-CMC_umap_sex.png")
  SZBD_sexumap <- ggdraw() + draw_image("results/9.6-SZBD_umap_sex.png")
  MB_sexumap <- ggdraw() + draw_image("results/9.9-MB_umap_sex.png")
  
  sex_dummy <- data.frame(x = 1:2, y = 1:2, sex = c("Male", "Female"))
  sex_dummyplot <- ggplot(sex_dummy) +
    geom_point(aes(x = x, y = y, color = sex), size = 3) +
    labs(color = "Biological \n sex") +
    scale_color_manual(values = c("Male" = "hotpink1", "Female" = "seagreen2")) +
    theme_minimal()
  
  sexlegend <- get_legend(sex_dummyplot)
  sex_combined_plot <- plot_grid(CMC_sexumap, SZBD_sexumap, MB_sexumap, sexlegend, ncol = 4)
  
  ggsave(file.path("results/9.12-sex_umap.png"), sex_combined_plot, dpi = 300, width = 8, height = 6)
  
  ### age distributions
  clean_meta <- readRDS("data/data_processed/clean_metadata.rds") |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "CMC" | Cohort == "SZBDMulti-Seq" | Cohort == "MultiomeBrain") |>
    filter(Disorder == "Schizophrenia" | Disorder == "Control") |>
    select(Cohort, Biological_Sex, Disorder, Age)
  
  median_age <- median(clean_meta$Age, na.rm = TRUE)
  
  disorder_median_values <- clean_meta %>%
    group_by(Disorder) %>%
    summarize(median_age = median(Age, na.rm = TRUE))
  
  sex_median_values <- clean_meta %>%
    group_by(Biological_Sex) %>%
    summarize(median_age = median(Age, na.rm = TRUE))
  
  age_hist <- clean_meta |>
    ggplot(aes(x = Age)) +
    geom_histogram(fill = "#948dd1") +
    geom_vline(xintercept = median_age, linetype = "dashed", color = "navy", size = 1) +
    labs(x = "Age at death (years)",
         y = "Number of patients") +
    annotate("text", x = 15, y = 10, label = paste("Median =", median_age, "years"), 
             color = "black", hjust = 0) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(10, 10, 10, 10)) +
    scale_x_continuous(breaks= c(10, 20, 30, 40, 50, 60, 70, 80, 90))
  
  age_hist_disorder <- clean_meta |>
    ggplot(aes(x = Age)) +
    geom_histogram(aes(fill = Disorder, color = Disorder), alpha = 0.6) +
    labs(x = "Age at death (years)",
         y = "Number of patients") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_x_continuous(breaks= c(10, 20, 30, 40, 50, 60, 70, 80, 90)) +
    scale_fill_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    geom_vline(data = disorder_median_values, aes(xintercept = median_age, color = Disorder), linetype = "dashed", size = 1) +
    annotate("text", x = 40, y = 15, 
             label = paste("Median:", round(disorder_median_values$median_age[1], 1), "years"), 
             color = "black") +
    annotate("text", x = 40, y = 12, 
             label = paste("Median:", round(disorder_median_values$median_age[2], 1), "years"), 
             color = "firebrick")
  
  age_hist_sex <- clean_meta |>
    ggplot(aes(x = Age)) +
    geom_histogram(aes(fill = Biological_Sex, color = Biological_Sex), alpha = 0.6) +
    labs(x = "Age at death (years)",
         y = "Number of patients",
         fill = "Biological Sex",
         color = "Biological Sex") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_x_continuous(breaks= c(10, 20, 30, 40, 50, 60, 70, 80, 90)) +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    scale_fill_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    geom_vline(data = sex_median_values, aes(xintercept = median_age, color = Biological_Sex), linetype = "dashed", size = 1) +
    annotate("text", x = 40, y = 15, 
             label = paste("Median:", round(sex_median_values$median_age[1], 1), "years"), 
             color = "seagreen") +
    annotate("text", x = 40, y = 12, 
             label = paste("Median:", round(sex_median_values$median_age[2], 1), "years"), 
             color = "hotpink3")
  
  control_quantiles <- quantile(clean_meta$Age[clean_meta$Disorder == "Control"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  schizophrenia_quantiles <- quantile(clean_meta$Age[clean_meta$Disorder == "Schizophrenia"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  
  qq_data <- data.frame(
    Control_Ages = control_quantiles,
    Schizophrenia_Ages = schizophrenia_quantiles
  )
  qq_plot <- ggplot(qq_data, aes(x = Control_Ages, y = Schizophrenia_Ages)) +
    geom_point(color = "#948dd1", size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "navy", linetype = "dashed") +
    labs(x = "Quantiles of Control Ages",
         y = "Quantiles of Schizophrenia Ages") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(10, 10, 10, 10))
  
  disorder_violin <- clean_meta |>
    ggplot(aes(x = Disorder, y = Age)) +
    geom_violin(aes(color = Disorder), trim = FALSE) +
    geom_jitter(aes(color = Disorder), alpha = 0.6, width = 0.2) +
    geom_boxplot(width = 0.1, aes(color = Disorder)) +
    labs(x = NULL, y = "Age") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_fill_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2"))
  
  sex_violin <- clean_meta |>
    ggplot(aes(x = Biological_Sex, y = Age)) +
    geom_violin(aes(color = Biological_Sex), trim = FALSE) +
    geom_jitter(aes(color = Biological_Sex), alpha = 0.6, width = 0.2) +
    geom_boxplot(width = 0.1, aes(color = Biological_Sex)) +
    labs(x = NULL, y = "Age") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "mediumseagreen")) +
    scale_fill_manual(values = c("male" = "hotpink1", "female" = "mediumseagreen"))
  
  age_dist_combined_plot <- plot_grid(age_hist, qq_plot, age_hist_disorder, disorder_violin, age_hist_sex, sex_violin, ncol = 2, rel_heights = c(1, 1.2, 1))
  
  ggsave(file.path("results/10-age_distributions.png"), age_dist_combined_plot, dpi = 300, width = 10, height = 12)
  
  ### cell violin and qq plots
  CMC_meta_count <- CMCcombined_meta |>
    select(patientID, disorder, sex) |>
    group_by(patientID, disorder, sex) |>
    summarize(total_cells = n()) |>
    mutate(cohort = "CMC")
  
  SZBD_meta_count <- SZBDcombined_meta |>
    select(patientID, disorder, sex) |>
    group_by(patientID, disorder, sex) |>
    summarize(total_cells = n())|>
    mutate(cohort = "SZBDMulti-Seq")
  
  MB_meta_count <- MBcombined_meta |>
    select(patientID, disorder, sex) |>
    group_by(patientID, disorder, sex) |>
    summarize(total_cells = n())|>
    mutate(cohort = "MultiomeBrain")
  
  meta_count <- rbind(CMC_meta_count, SZBD_meta_count, MB_meta_count)
  
  disorder_cellcount_violin <- meta_count |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = disorder, y = total_cells)) +
    geom_violin(aes(color = disorder), trim = FALSE) +
    geom_jitter(aes(color = disorder), alpha = 0.6, width = 0.2) +
    geom_boxplot(width = 0.1, aes(color = disorder)) +
    labs(x = NULL, y = "Number of cells per patient", color = "Disorder", fill = "Disorder") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.spacing = unit(0, "lines"),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_fill_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    facet_wrap(~cohort)
  
  sex_cellcount_violin <- meta_count |>
    ggplot(aes(x = sex, y = total_cells)) +
    geom_violin(aes(color = sex), trim = FALSE) +
    geom_jitter(aes(color = sex), alpha = 0.6, width = 0.2) +
    geom_boxplot(width = 0.1, aes(color = sex)) +
    labs(x = NULL, y = "Number of cells per patient", color = "Biological Sex", fill = "Biological Sex") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.spacing = unit(0, "lines"),
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10)) +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "mediumseagreen")) +
    scale_fill_manual(values = c("male" = "hotpink1", "female" = "mediumseagreen")) +
    facet_wrap(~cohort)
  
  control_cell_quantiles <- quantile(meta_count$total_cells[meta_count$disorder == "no"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  schizophrenia_cell_quantiles <- quantile(meta_count$total_cells[meta_count$disorder == "yes"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  
  qq_data_cellcount <- data.frame(
    Control_cell_counts = control_cell_quantiles,
    Schizophrenia_cell_counts = schizophrenia_cell_quantiles
  )
  
  qq_plot_disorder <- ggplot(qq_data_cellcount, aes(x = Control_cell_counts, y = Schizophrenia_cell_counts)) +
    geom_point(color = "#948dd1", size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "navy", linetype = "dashed") +
    labs(x = "Quantiles of Control Cell Counts",
         y = "Quantiles of Schizophrenia Cell Counts") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(10, 10, 10, 10))
  
  female_cell_quantiles <- quantile(meta_count$total_cells[meta_count$sex == "female"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  male_cell_quantiles <- quantile(meta_count$total_cells[meta_count$sex == "male"], probs = seq(0, 1, 0.01), na.rm = TRUE)
  
  qq_sex_cellcount <- data.frame(
    female_cell_counts = female_cell_quantiles,
    male_cell_counts = male_cell_quantiles
  )
  
  qq_plot_sex <- ggplot(qq_sex_cellcount, aes(x = female_cell_counts, y = male_cell_counts)) +
    geom_point(color = "#948dd1", size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "navy", linetype = "dashed") +
    labs(x = "Quantiles of Female Cell Counts",
         y = "Quantiles of Male Cell Counts") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = margin(10, 10, 10, 10))
  
  cellcount_combined_plot <- plot_grid(disorder_cellcount_violin, qq_plot_disorder, sex_cellcount_violin, qq_plot_sex, ncol = 2, rel_heights = c(1, 1.2, 1))
  print(cellcount_combined_plot)
  
  ggsave(file.path("results/11-cellcount_plots.png"), cellcount_combined_plot, dpi = 300, width = 12, height = 12)
}

main()
