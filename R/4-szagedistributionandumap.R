# Author: Rui Xiang Yu
# Date: 2024 October 21st
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/4-szagedistributionandumap.R

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(Seurat)
library(future.apply)

main <- function() {
 
  ### age distrbutions
   CMC_meta <- readRDS("data/data_processed/CMC/CMC-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  MB_meta <- readRDS("data/data_processed/MultiomeBrain/MultiomeBrain-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  SZBD_meta <- readRDS("data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  
  CMC_age_hist <- CMC_meta |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = age)) +
    geom_histogram(fill = "blue4") +
    labs(x = "Age at death (years)",
         y = "Number of patients",
         title = "CMC") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))
  
  print(CMC_age_hist)
  
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
}

main()
