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

main <- function() {
 
  ### age distrbutions
  #  CMC_meta <- readRDS("data/data_processed/CMC/CMC-patients.rds") |>
  #   mutate(age = as.numeric(str_replace(age, "\\+", "")))
  # MB_meta <- readRDS("data/data_processed/MultiomeBrain/MultiomeBrain-patients.rds") |>
  #   mutate(age = as.numeric(str_replace(age, "\\+", "")))
  # SZBD_meta <- readRDS("data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patients.rds") |>
  #   mutate(age = as.numeric(str_replace(age, "\\+", "")))
  # 
  # CMC_age_hist <- CMC_meta |>
  #   mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
  #   ggplot(aes(x = age)) +
  #   geom_histogram(fill = "blue4") +
  #   labs(x = "Age at death (years)",
  #        y = "Number of patients",
  #        title = "CMC") +
  #   theme_minimal() +
  #   theme(
  #     panel.border = element_blank(),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank(),
  #     axis.line = element_line(colour = "black"))
  # 
  # print(CMC_age_hist)
  
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
  
  #CMCtop10 <- head(VariableFeatures(CMC_seurat), 10)
  #CMCplot1 <- VariableFeaturePlot(CMC_seurat, cols = c("grey42", "purple3")) + labs(title = "CMC")
  #CMCplot2 <- LabelPoints(plot = CMCplot1, points = CMCtop10, repel = TRUE)
  
  CMCall.genes <- rownames(CMC_seurat)
  CMC_seurat <- ScaleData(CMC_seurat, features = CMCall.genes)
  CMC_seurat <- RunPCA(CMC_seurat, features = VariableFeatures(object = CMC_seurat))
  
  #CMC_dimloadings <- VizDimLoadings(CMC_seurat, dims = 1:2, reduction = "pca")
  #CMC_elbow <- ElbowPlot(CMC_seurat)
  
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
  
  #SZBDtop10 <- head(VariableFeatures(SZBD_seurat), 10)
  #SZBDplot1 <- VariableFeaturePlot(SZBD_seurat, cols = c("grey42", "purple3")) + labs(title = "SZBD")
  #SZBDplot2 <- LabelPoints(plot = SZBDplot1, points = SZBDtop10, repel = TRUE)
  
  SZBDall.genes <- rownames(SZBD_seurat)
  SZBD_seurat <- ScaleData(SZBD_seurat, features = SZBDall.genes)
  SZBD_seurat <- RunPCA(SZBD_seurat, features = VariableFeatures(object = SZBD_seurat))
  
  #SZBD_dimloadings <- VizDimLoadings(SZBD_seurat, dims = 1:2, reduction = "pca")
  SZBD_elbow <- ElbowPlot(SZBD_seurat)
  
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
  
  #MBtop10 <- head(VariableFeatures(MB_seurat), 10)
  #MBplot1 <- VariableFeaturePlot(MB_seurat, cols = c("grey42", "purple3")) + labs(title = "MB")
  #MBplot2 <- LabelPoints(plot = MBplot1, points = MBtop10, repel = TRUE)
  
  MBall.genes <- rownames(MB_seurat)
  MB_seurat <- ScaleData(MB_seurat, features = MBall.genes)
  MB_seurat <- RunPCA(MB_seurat, features = VariableFeatures(object = MB_seurat))
  
  #MB_dimloadings <- VizDimLoadings(MB_seurat, dims = 1:2, reduction = "pca")
  MB_elbow <- ElbowPlot(MB_seurat)
  
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
  
  ### Overall umap
  expr_list <- list(CMCcombined_expr, SZBDcombined_expr, MBcombined_expr)
  
  all_genes <- Reduce(union, lapply(expr_list, rownames))
  pad_genes <- function(expr, all_genes) {
    missing_genes <- setdiff(all_genes, rownames(expr))
    zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(expr))
    rownames(zero_mat) <- missing_genes
    expr_padded <- rbind(expr, zero_mat)
    return(expr_padded[all_genes, , drop = FALSE])
  }
  
  expr_list_padded <- lapply(expr_list, pad_genes, all_genes = all_genes)
  combined_expr <- do.call(cbind, expr_list_padded)
  
  metadata_list <- list(CMCcombined_meta, SZBDcombined_meta, MBcombined_meta)
  combined_meta <- do.call(rbind, metadata_list)
  
  all_seurat <- CreateSeuratObject(counts = combined_expr, meta.data = combined_meta)
  all_seurat <- NormalizeData(all_seurat)
  all_seurat <- FindVariableFeatures(all_seurat, selection.method = "vst", nfeatures = 2000)
  
  allall.genes <- rownames(all_seurat)
  all_seurat <- ScaleData(all_seurat, features = allall.genes)
  all_seurat <- RunPCA(all_seurat, features = VariableFeatures(object = all_seurat))
  all_elbow <- ElbowPlot(all_seurat)
  ggsave(file.path("results/elbow.png"), all_elbow)
  
  all_seurat <- RunUMAP(all_seurat, dims = 1:17)
  all_umapplot <- DimPlot(all_seurat, reduction = "umap", group.by = "cell_type") +
    scale_color_manual(values = c("Ast" = "turquoise", 
                                  "Opc" = "coral", 
                                  "Exc" = "pink2", 
                                  "Inh" = "olivedrab3", 
                                  "Mic" = "navy",
                                  "Oli" = "yellow2")) +
    xlab(NULL) + 
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("all")
  
  all_disorderumapplot <- DimPlot(all_seurat, reduction = "umap", group.by = "disorder") +
    scale_color_manual(values = c("no" = "grey21", "yes" = "firebrick2")) +
    xlab(NULL) + 
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("all")
  
  all_sexumapplot <- DimPlot(all_seurat, reduction = "umap", group.by = "sex") +
    scale_color_manual(values = c("male" = "hotpink1", "female" = "seagreen2")) +
    xlab(NULL) + 
    ylab(NULL) +
    theme(legend.position = "none") +
    ggtitle("all")
  
  ggsave(file.path("results/9.10-all_umap.png"), all_umapplot)
  ggsave(file.path("results/9.11-all_umap_disorder.png"), all_disorderumapplot)
  ggsave(file.path("results/9.12-all_umap_sex.png"), all_sexumapplot)
  
}

main()
