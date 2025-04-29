# Author: Rui Xiang Yu
# Date: 2025 March 4th
# This script compares Ling results.
# Usage: R/9-ling.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)

main <- function() {

studies <- c("Batiuk", "CMC", "Ling", "MultiomeBrain", "SZBDMulti-Seq")

ll<-read.table("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt", header=T, sep='\t')

ll.ex<-ll[ grepl(x=ll$id,pattern=".+_glutamatergic") ,]
ll.ex[,1]<-gsub(ll.ex[,1], pattern="_glutamatergic", r="")
ll.ex.lf4<-ll.ex[,c(1,5)] # isolate factor 4
row.names(ll.ex.lf4)<-ll.ex.lf4[,1]
ll.ex.lf4<-ll.ex.lf4[,2, drop=F]

ll.in<-ll[ grepl(x=ll$id,pattern=".+_gabaergic") ,]
ll.in[,1]<-gsub(ll.in[,1], pattern="_gabaergic", r="")
ll.in.lf4<-ll.in[,c(1,5)] # isolate factor 4
row.names(ll.in.lf4)<-ll.in.lf4[,1]
ll.in.lf4<-ll.in.lf4[,2, drop=F]

ll.as<-ll[grepl(x=ll$id,pattern=".+_astrocyte"),]
ll.as[,1]<-gsub(ll.as[,1], pattern="_astrocyte", r="")
ll.as.lf4<-ll.as[,c(1,5)]
row.names(ll.as.lf4)<-ll.as.lf4[,1]
ll.as.lf4<-ll.as.lf4[,2, drop=F]

cell_types <- c("Exc", "Inh", "Ast")

for (study in studies){
  all_plots <- list()
for (cell_type in cell_types){
  res_path <- paste0("results/DEA/", study, "/LimmaCPMLogAgeSex/DEAresults_", cell_type, "_", study, "_SZ.rds")

  if (!file.exists(res_path)) {
    message("File ", res_path, " does not exist. Skipping to next.")
    next
  }

  des_res <- readRDS(res_path)
  des_res <- des_res[des_res$adj.P.Val < 0.05, ]

  deg_up <- des_res[des_res$groupdisorderyes > 0.5, ] |> rownames() |> unique()
  deg_down <- des_res[des_res$groupdisorderyes < -0.5, ] |> rownames() |> unique()

  if (cell_type == "Ast"){
    p_up <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="Astrocytes",
           subtitle = "Rug: top up-regulated genes")

    p_up <-p_up + geom_rug(data=ll.as.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

    p_down <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="Astrocytes",
           subtitle = "Rug: top down-regulated genes")

    p_down <-p_down + geom_rug(data=ll.as.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
  }

  if (cell_type == "Exc"){
    p_up <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="Ex. Neurons",
           subtitle = "Rug: top up-regulated genes")

    p_up <-p_up + geom_rug(data=ll.ex.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

    p_down <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="Ex. Neurons",
           subtitle = "Rug: top down-regulated genes")

    p_down <-p_down + geom_rug(data=ll.ex.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
  }

  if (cell_type == "Inh"){
    p_up <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="In. Neurons",
           subtitle = "Rug: top up-regulated genes")

    p_up <-p_up + geom_rug(data=ll.in.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

    p_down <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
      geom_histogram(bins=100) +
      theme_classic() +
      labs(title="In. Neurons",
           subtitle = "Rug: top down-regulated genes")

    p_down <-p_down + geom_rug(data=ll.in.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
  }

  all_plots[[paste0(cell_type, "_up")]] <- p_up
  all_plots[[paste0(cell_type, "_down")]] <- p_down
  }

  plot_list <- plot_grid(plotlist = all_plots, ncol = 2, align = "hv")

  title_plot <- ggdraw() +
    draw_label(paste0(study, ": disease contrast"), x = 0.5, y = 0.5, size = 12, hjust = 0.5)
  all_plots_with_title <- plot_grid(title_plot, plot_list, ncol = 1, rel_heights = c(0.1, 1))

  final_path <- paste0("results/Ling/", study, "/initialrug_disease.png")
  ggsave(plot = all_plots_with_title, filename = final_path, height = 12, width = 10)
}

for (study in studies){
  all_plots <- list()
  for (cell_type in cell_types){
    res_path <- paste0("results/DEA/", study, "/LimmaCPMLogAgeSex/DEAresults_", cell_type, "_", study, "_SZ.rds")

    if (!file.exists(res_path)) {
      message("File ", res_path, " does not exist. Skipping to next.")
      next
    }

    des_res <- readRDS(res_path)
    des_res <- des_res[des_res$adj.P.Val < 0.05, ]

    deg_up <- des_res[des_res$age > 0.05, ] |> rownames() |> unique()
    deg_down <- des_res[des_res$age < -0.05, ] |> rownames() |> unique()

    if (cell_type == "Ast"){
      p_up <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Astrocytes",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.as.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Astrocytes",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.as.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }

    if (cell_type == "Exc"){
      p_up <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Ex. Neurons",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.ex.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Ex. Neurons",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.ex.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }

    if (cell_type == "Inh"){
      p_up <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="In. Neurons",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.in.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="In. Neurons",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.in.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }

    all_plots[[paste0(cell_type, "_up")]] <- p_up
    all_plots[[paste0(cell_type, "_down")]] <- p_down
  }

  plot_list <- plot_grid(plotlist = all_plots, ncol = 2, align = "hv")

  title_plot <- ggdraw() +
    draw_label(paste0(study, ": age contrast"), x = 0.5, y = 0.5, size = 12, hjust = 0.5)
  all_plots_with_title <- plot_grid(title_plot, plot_list, ncol = 1, rel_heights = c(0.1, 1))

  final_path <- paste0("results/Ling/", study, "/initialrug_age.png")
  ggsave(plot = all_plots_with_title, filename = final_path, height = 12, width = 10)
}

### save ling as pbs
ling_meta <- read_csv("data/data_processed/Ling/all-patient-metaData.csv") |>
  rename(disorder = diagnosis) |>
  filter(study == "ling-2024") |>
  dplyr::select(patientID, sex, age, disorder, brainRegion, study)

saveRDS(ling_meta, "data/data_processed/Ling/Ling-patient.rds")

cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc")

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

  cell_meta <- cell_meta[colnames(cell_expr), , drop = FALSE]

  PB <- list(expr = cell_expr, meta = cell_meta)

  pb_path <- paste0("data/data_processed/Ling/Pseudobulk/PseudobulkRaw/", celltype, "_Ling_SZ.rds")
  saveRDS(PB, pb_path)
}

### combine matrices
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc", "Gli")
for (study in studies){

  all_matrix <- list()

  for (cell_type in cell_types){
    pb_path <- paste0("data/data_processed/", study, "/Pseudobulk/PseudobulkRaw/", cell_type, "_", study, "_SZ.rds")

    if (!file.exists(pb_path)) {
      message("Skipping missing file: ", pb_path)
      next
    }

    PB <- readRDS(pb_path)

    matrix_genes <- rownames(PB$expr)

    matrix_append <- paste0(matrix_genes, "_", cell_type)

    rownames(PB$expr) <- matrix_append

    all_matrix[[cell_type]] <- PB$expr
  }

  all_columns <- unique(unlist(lapply(all_matrix, colnames)))

  for (cell_type in names(all_matrix)) {
    missing_cols <- setdiff(all_columns, colnames(all_matrix[[cell_type]]))
    if (length(missing_cols) > 0) {
      zero_matrix <- Matrix(0, nrow = nrow(all_matrix[[cell_type]]), ncol = length(missing_cols), sparse = TRUE)
      colnames(zero_matrix) <- missing_cols
      all_matrix[[cell_type]] <- cbind(all_matrix[[cell_type]], zero_matrix)
    }

    all_matrix[[cell_type]] <- all_matrix[[cell_type]][, all_columns]
  }

  final_matrix <- do.call(rbind, all_matrix)
  matrixmtx <- as.matrix(final_matrix)
  final_dgCmatrix <- as(matrixmtx, "dgCMatrix")
  saveRDS(final_dgCmatrix, paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.rds"))
  gc()
}

studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq")

for (study in studies){
  matrix_path <- paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.rds")
  combined_matrix <- readRDS(matrix_path)

  control_cols <- colnames(combined_matrix)[str_detect(colnames(combined_matrix), "disorderno")]
  scz_cols <- colnames(combined_matrix)[str_detect(colnames(combined_matrix), "disorderyes")]

  control_matrix <- combined_matrix[, control_cols, drop = FALSE]
  scz_matrix <- combined_matrix[, scz_cols, drop = FALSE]

  saveRDS(control_matrix, paste0("data/data_processed/", study, "/Pseudobulk/", study, "_control_matrix.rds"))
  saveRDS(scz_matrix, paste0("data/data_processed/", study, "/Pseudobulk/", study, "_scz_matrix.rds"))
}


  combined_matrix <- readRDS("data/data_processed/Ling/Pseudobulk/Ling_combined_matrix.rds")

  meta_control <- ling_meta |>
    filter(disorder == "no") |>
    dplyr::select(patientID) |>
    pull()

  meta_scz <- ling_meta |>
    filter(disorder == "yes") |>
    dplyr::select(patientID) |>
    pull()

  control_matrix <- combined_matrix[, colnames(combined_matrix) %in% meta_control]
  scz_matrix <- combined_matrix[, colnames(combined_matrix) %in% meta_scz]

  saveRDS(control_matrix, paste0("data/data_processed/Ling/Pseudobulk/Ling_control_matrix.rds"))
  saveRDS(scz_matrix, paste0("data/data_processed/Ling/Pseudobulk/Ling_scz_matrix.rds"))

### matrices for peer
  studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")
  for (study in studies){
    mtx_path <- paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.rds")
    comb_mtx <- readRDS(mtx_path)

    cpm_mtx <- do_cpm_log(comb_mtx, log = TRUE)
    cpm_mtx <- as.data.frame(t(cpm_mtx))
    write_csv(cpm_mtx, paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.csv"))
  }

### covs
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")

for (study in studies){
  meta_path <- paste0("data/data_processed/", study, "/", study, "-patient.rds")
  meta <- readRDS(meta_path) |>
    dplyr::select(patientID, sex, age, disorder) |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))

  meta$disorder <- as.numeric(meta$disorder == "yes")

  matrix_path <- paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.rds")
  combined_matrix <- readRDS(matrix_path)

  correct_order <- colnames(combined_matrix) |>
    str_remove(":(disorderyes|disorderno)") |>
    str_remove("^patientID")

  meta <- meta |> filter(patientID %in% correct_order)

  meta <- meta |> arrange(factor(patientID, levels = correct_order))

  meta <- meta |>
    dplyr::select(-patientID)

  covs_path <- paste0("data/data_processed/", study, "/", study, "-peercovs.csv")
  write_csv(meta, covs_path)
}

for (study in studies){
  covs_path <- paste0("data/data_processed/", study, "/", study, "-peercovs.csv")
  covs <- fread(covs_path)

  disorder_only <- covs |> dplyr::select(disorder)
  write_csv(disorder_only, paste0("data/data_processed/", study, "/", study, "-disorder.csv"))

  age_only <- covs |> dplyr::select(age)
  write_csv(age_only, paste0("data/data_processed/", study, "/", study, "-age.csv"))

  agedisorder <- covs |>dplyr::select(disorder, age)
  write_csv(agedisorder, paste0("data/data_processed/", study, "/", study, "-agedisorder.csv"))
}

### pca
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")
for (study in studies){
  mtx_path <- paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_matrix.rds")
  comb_mtx <- readRDS(mtx_path)

  cpm_mtx <- do_cpm_log(comb_mtx, log = TRUE)
  cpm_mtx <- as.data.frame(t(cpm_mtx))

  cpm_mtx <- as.matrix(cpm_mtx)
  cpm_mtx <- as(cpm_mtx, "dgCMatrix")

  saveRDS(cpm_mtx, paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))

  PCAmtx <- prcomp(cpm_mtx)
  saveRDS(PCAmtx, paste0("results/Ling/", study, "/", study, "_PCA.rds"))
}

### peer vs pca
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")
peers <- c("nocovs", "peercovs", "justage", "agedisorder", "justdisorder")

for (study in studies){

  all_list <- list()
  only4_list <- list()

  for (peer in peers){
    w_matrix <- fread(paste0("results/Ling/", study, "/", peer, "/W.csv")) |> t()
    pca_matrix <- readRDS(paste0("results/Ling/", study, "/", study, "_PCA.rds"))
    pca_scores <- pca_matrix$rotation
    pca_scores <- pca_scores[-1, 1:10]

  if (study == "MultiomeBrain"){
    if (peer == "nocovs"){
      pca_scores <- pca_scores[, 1:9]
    }

    if (peer == "peercovs"){
      pca_scores <- pca_scores[, 1:6]
    }
    if (peer == "justage" | peer == "justdisorder"){
      pca_scores <- pca_scores[, 1:8]
    }
    if (peer == "agedisorder"){
      pca_scores <- pca_scores[, 1:7]
    }
  }
      if (peer == "peercovs"){
        w_matrix <- w_matrix[, -c(1, 2, 3)]
      }
      if (peer == "justage" | peer == "justdisorder"){
        w_matrix <- w_matrix[, -1]
      }
      if (peer == "agedisorder"){
        w_matrix <- w_matrix[, -c(1, 2)]
      }

      mat1_vec <- as.vector(scale(w_matrix))
      mat2_vec <- as.vector(scale(pca_scores))

      differences <- abs(mat1_vec - mat2_vec)
      diff_sd <- sd(differences, na.rm = TRUE)
      threshold <- diff_sd

      df <- data.frame(Mat1 = mat1_vec, Mat2 = mat2_vec)

      df$distance <- abs(df$Mat1 - df$Mat2)
      df$category <- ifelse(df$distance <= threshold, "Close to Line", "Far from Line")
      percentage_close <- sum(df$category == "Close to Line") / nrow(df) * 100

      all_plot <- ggplot(df, aes(x = Mat1, y = Mat2, color = category)) +
        geom_point(alpha = 0.3) +
        labs(title = peer,
             x = "PEER Values",
             y = "PCA Values",
             color = "Category") +
        theme_classic()+
        geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")+
        annotate("text", x = max(df$Mat1, na.rm = TRUE),
                 y = min(df$Mat2, na.rm = TRUE),
                 label = paste0(round(percentage_close, 2), "%"),
                 hjust = 1, vjust = 0, size = 5, color = "red") +
        theme(legend.position="none")+
        scale_color_manual(values = c(
          "Close to Line" = "red",
          "Far from Line" = "black"
        ))

      all_list[[peer]] <- all_plot

      peer4 <- as.vector(scale(w_matrix[,4 ]))
      pca4 <- as.vector(scale(pca_scores[,4 ]))

      differences <- abs(peer4 - pca4)
      diff_sd <- sd(differences, na.rm = TRUE)
      threshold <- diff_sd

      df4 <- data.frame(Mat1 = peer4, Mat2 = pca4)

      df4$distance <- abs(df4$Mat1 - df4$Mat2)
      df4$category <- ifelse(df4$distance <= threshold, "Close to Line", "Far from Line")
      percentage_close <- sum(df4$category == "Close to Line") / nrow(df4) * 100

      only4 <- ggplot(df4, aes(x = Mat1, y = Mat2, color = category)) +
        geom_point(alpha = 0.3) +
        labs(title = peer,
             x = "PEER4 Values",
             y = "PCA4 Values") +
        theme_classic()+
        geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")+
        annotate("text", x = max(df4$Mat1, na.rm = TRUE),
                 y = min(df4$Mat2, na.rm = TRUE),
                 label = paste0(round(percentage_close, 2), "%"),
                 hjust = 1, vjust = 0, size = 5, color = "red") +
        theme(legend.position="none")+
        scale_color_manual(values = c(
          "Close to Line" = "red",
          "Far from Line" = "black"
        ))

      only4_list[[peer]] <- only4
      gc()
  }
  all_plots <- plot_grid(
    plotlist = all_list,
    ncol = 3,
    align = "hv")

  ggsave(plot = all_plots, filename = paste0("results/Ling/", study, "/peervspcaALL.png"))

  only4_plots <- plot_grid(
    plotlist = only4_list,
    ncol = 3,
    align = "hv")

  ggsave(plot = only4_plots, filename = paste0("results/Ling/", study, "/peervspcaONLY4.png"))
}

### ling vs peer
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")
peers <- c("nocovs", "peercovs", "justage", "agedisorder", "justdisorder")

for (study in studies){
  all_plots <- list()
  only4_plots <- list()
  for (peer in peers){
    ling_loadings <- fread("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt") |>
      separate(id, into=c("gene", "cell_type"), sep = "_") |>
      mutate(cell_type = ifelse(cell_type == "glutamatergic", "Exc", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "astrocyte", "Ast", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "gabaergic", "Inh", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "oligodendrocyte", "Oli", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "microglia", "Mic", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "polydendrocyte", "Opc", cell_type)) |>
      filter(cell_type != "endothelia")

    ling_loadings$gene_names <- paste0(ling_loadings$gene, "_", ling_loadings$cell_type)
    ling_loadings <- ling_loadings |>
      dplyr::select(-gene, -cell_type, -V12, -V13, -V14, -V15)

    correct_genes <- ling_loadings$gene_names

    cpm_pb <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
    gene_names <- colnames(cpm_pb)

    if (study != "Ling"){
    gene_names <- sub("^[^.]+\\.", "", gene_names)}
    gene_names <- gene_names[-1]

    w_matrix <- fread(paste0("results/Ling/", study, "/", peer, "/W.csv")) |> t() |> data.frame()
    rownames(w_matrix) <- gene_names

    common_genes <- intersect(correct_genes, gene_names)
    ling_loadings <- ling_loadings[ling_loadings$gene_names %in% common_genes, ]
    w_matrix <- w_matrix[rownames(w_matrix) %in% common_genes,]
    rownames(ling_loadings) <- ling_loadings$gene_names
    ling_loadings <- ling_loadings[, !(names(ling_loadings) %in% "gene_names")]
    ling_loadings <- ling_loadings[rownames(w_matrix), , drop = FALSE]

    if (study == "MultiomeBrain"){
      if (peer == "nocovs"){
        ling_loadings <- ling_loadings[, 1:9]
      }

      if (peer == "peercovs"){
        ling_loadings <- ling_loadings[, 1:6]
      }
      if (peer == "justage" | peer == "justdisorder"){
        ling_loadings <- ling_loadings[, 1:8]
      }
      if (peer == "agedisorder"){
        ling_loadings <- ling_loadings[, 1:7]
      }
    }
    if (peer == "peercovs"){
      w_matrix <- w_matrix[, -c(1, 2, 3)]
    }
    if (peer == "justage" | peer == "justdisorder"){
      w_matrix <- w_matrix[, -1]
    }
    if (peer == "agedisorder"){
      w_matrix <- w_matrix[, -c(1, 2)]
    }

    mat1_vec <- as.vector(scale((as.matrix(w_matrix))))
    mat2_vec <- as.vector(scale(as.matrix(ling_loadings)))

    differences <- abs(mat1_vec - mat2_vec)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    df <- data.frame(Mat1 = mat1_vec, Mat2 = mat2_vec)

    df$distance <- abs(df$Mat1 - df$Mat2)
    df$category <- ifelse(df$distance <= threshold, "Close to Line", "Far from Line")
    percentage_close <- sum(df$category == "Close to Line") / nrow(df) * 100

    all_plot <- ggplot(df, aes(x = Mat1, y = Mat2, color = category)) +
      geom_point(alpha = 0.5) +
      labs(title = peer,
           x = "Rui PEER",
           y = "Ling PEER",
           color = "Category") +  # Legend label
      theme_classic() +
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
      scale_color_manual(values = c("Close to Line" = "red", "Far from Line" = "black")) +
      annotate("text", x = max(df$Mat1, na.rm = TRUE),
               y = min(df$Mat2, na.rm = TRUE),
               label = paste0(round(percentage_close, 2), "%"),
               hjust = 1, vjust = 0, size = 5, color = "red") +
      theme(legend.position="none")

    all_plots[[peer]] <- all_plot

    peer4 <- as.vector(scale(w_matrix[,4 ]))
    ling4 <- as.vector(scale(ling_loadings[,4 ]))

    differences <- abs(peer4 - ling4)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    df4 <- data.frame(Mat1 = peer4, Mat2 = ling4)

    df4$distance <- abs(df4$Mat1 - df4$Mat2)
    df4$category <- ifelse(df4$distance <= threshold, "Close to Line", "Far from Line")
    percentage_close <- sum(df4$category == "Close to Line") / nrow(df4) * 100

    only4 <- ggplot(df4, aes(x = Mat1, y = Mat2, color = category)) +
      geom_point(alpha = 0.3) +
      labs(title = peer,
           x = "Rui PEER4",
           y = "Ling PEER4") +
      theme_classic()+
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")+
      annotate("text", x = max(df4$Mat1, na.rm = TRUE),
               y = min(df4$Mat2, na.rm = TRUE),
               label = paste0(round(percentage_close, 2), "%"),
               hjust = 1, vjust = 0, size = 5, color = "red") +
      theme(legend.position="none")+
      scale_color_manual(values = c("Close to Line" = "red", "Far from Line" = "black"))

    only4_plots[[peer]] <- only4
    gc()
  }
  allplots <- plot_grid(
    plotlist = all_plots,
    ncol = 3,
    align = "hv")

  ggsave(plot = allplots, filename = paste0("results/Ling/", study, "/peervslingALL.png"))

  only4plots <- plot_grid(
    plotlist = only4_plots,
    ncol = 3,
    align = "hv")

  ggsave(plot = only4plots, filename = paste0("results/Ling/", study, "/peervslingONLY4.png"))
}
### ling vs pca
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")

all_list <- list()
only4_list <- list()
for (study in studies){
    ling_loadings <- fread("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt") |>
      separate(id, into=c("gene", "cell_type"), sep = "_") |>
      mutate(cell_type = ifelse(cell_type == "glutamatergic", "Exc", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "astrocyte", "Ast", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "gabaergic", "Inh", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "oligodendrocyte", "Oli", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "microglia", "Mic", cell_type)) |>
      mutate(cell_type = ifelse(cell_type == "polydendrocyte", "Opc", cell_type)) |>
      filter(cell_type != "endothelia")

    ling_loadings$gene_names <- paste0(ling_loadings$gene, "_", ling_loadings$cell_type)
    ling_loadings <- ling_loadings |>
      dplyr::select(-gene, -cell_type, -V12, -V13, -V14, -V15)
    correct_genes <- ling_loadings$gene_names

    pca_matrix <- readRDS(paste0("results/Ling/", study, "/", study, "_PCA.rds"))
    pca_scores <- pca_matrix$rotation
    pca_scores <- pca_scores[, 1:10]
    if (study != "Ling"){
      gene_names <- rownames(pca_scores)
      gene_names <- sub("^[^.]+\\.", "", gene_names)
      rownames(pca_scores) <- gene_names
      }

    common_genes <- intersect(correct_genes, gene_names)
    ling_loadings <- ling_loadings[ling_loadings$gene_names %in% common_genes, ]
    pca_scores <- pca_scores[rownames(pca_scores) %in% common_genes,]
    rownames(ling_loadings) <- ling_loadings$gene_names
    ling_loadings <- ling_loadings[, !(names(ling_loadings) %in% "gene_names")]
    ling_loadings <- ling_loadings[rownames(pca_scores), , drop = FALSE]

    mat1_vec <- as.vector(scale(as.matrix(pca_scores)))
    mat2_vec <- as.vector(scale(as.matrix(ling_loadings)))

    differences <- abs(mat1_vec - mat2_vec)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    df <- data.frame(Mat1 = mat1_vec, Mat2 = mat2_vec)

    df$distance <- abs(df$Mat1 - df$Mat2)
    df$category <- ifelse(df$distance <= threshold, "Close to Line", "Far from Line")
    percentage_close <- sum(df$category == "Close to Line") / nrow(df) * 100

    all_plot <- ggplot(df, aes(x = Mat1, y = Mat2, color = category)) +
      geom_point(alpha = 0.5) +
      labs(title = study,
           x = "PCA",
           y = "Ling PEER",
           color = "Category") +  # Legend label
      theme_classic() +
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
      scale_color_manual(values = c("Close to Line" = "red", "Far from Line" = "black")) +
      annotate("text", x = max(df$Mat1, na.rm = TRUE),
               y = min(df$Mat2, na.rm = TRUE),
               label = paste0(round(percentage_close, 2), "%"),
               hjust = 1, vjust = 0, size = 5, color = "red") +
      theme(legend.position="none")

    all_list[[study]] <- all_plot

    pca4 <- as.vector(scale(pca_scores[,4 ]))
    ling4 <- as.vector(scale(ling_loadings[,4 ]))

    differences <- abs(pca4 - ling4)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    df4 <- data.frame(Mat1 = pca4, Mat2 = ling4)

    df4$distance <- abs(df4$Mat1 - df4$Mat2)
    df4$category <- ifelse(df4$distance <= threshold, "Close to Line", "Far from Line")
    percentage_close <- sum(df4$category == "Close to Line") / nrow(df4) * 100

    only4 <- ggplot(df4, aes(x = Mat1, y = Mat2, color = category)) +
      geom_point(alpha = 0.3) +
      labs(title = study,
           x = "PCA4",
           y = "Ling PEER4") +
      theme_classic()+
      geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")+
      annotate("text", x = max(df4$Mat1, na.rm = TRUE),
               y = min(df4$Mat2, na.rm = TRUE),
               label = paste0(round(percentage_close, 2), "%"),
               hjust = 1, vjust = 0, size = 5, color = "red") +
      theme(legend.position="none")+
      scale_color_manual(values = c("Close to Line" = "red", "Far from Line" = "black"))

    only4_list[[study]] <- only4
    gc()
}
all_plots <- plot_grid(
  plotlist = all_list,
  ncol = 3,
  align = "hv")

ggsave(plot = all_plots, filename = paste0("results/Ling/", study, "/lingvspcaALL.png"))

only4_plots <- plot_grid(
  plotlist = only4_list,
  ncol = 3,
  align = "hv")

ggsave(plot = only4_plots, filename = paste0("results/Ling/", study, "/lingvspcaONLY4.png"))

### get genes that appear across all studies in close to line
studies <- c("Batiuk", "CMC", "MultiomeBrain", "SZBDMulti-Seq", "Ling")
peers <- c("nocovs", "peercovs", "justage", "agedisorder", "justdisorder")

all_plots <- list()
only4_plots <- list()

ling_loadings_all <- fread("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt") %>%
  separate(id, into=c("gene", "cell_type"), sep = "_") %>%
  mutate(cell_type = recode(cell_type,
                            "glutamatergic" = "Exc",
                            "astrocyte" = "Ast",
                            "gabaergic" = "Inh",
                            "oligodendrocyte" = "Oli",
                            "microglia" = "Mic",
                            "polydendrocyte" = "Opc")) %>%
  filter(cell_type != "endothelia") %>%
  mutate(gene_names = paste0(gene, "_", cell_type)) %>%
  dplyr::select(-gene, -cell_type, -V12, -V13, -V14, -V15)

for (study in studies) {
  all_plots[[study]] <- list()
  only4_plots[[study]] <- list()

  for (peer in peers) {
    ling_loadings <- ling_loadings_all
    correct_genes <- ling_loadings$gene_names

    cpm_pb <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
    gene_names <- colnames(cpm_pb)
    if (study != "Ling") {
      gene_names <- sub("^[^.]+\\.", "", gene_names)
    }
    gene_names <- gene_names[-1]

    w_matrix <- fread(paste0("results/Ling/", study, "/", peer, "/W.csv")) %>% t() %>% data.frame()
    rownames(w_matrix) <- gene_names

    common_genes <- intersect(correct_genes, gene_names)
    ling_loadings <- ling_loadings[ling_loadings$gene_names %in% common_genes, ]
    w_matrix <- w_matrix[rownames(w_matrix) %in% common_genes, ]

    rownames(ling_loadings) <- ling_loadings$gene_names
    ling_loadings <- ling_loadings[, !(names(ling_loadings) %in% "gene_names")]
    ling_loadings <- ling_loadings[rownames(w_matrix), , drop = FALSE]

    if (study == "MultiomeBrain") {
      ling_loadings <- switch(peer,
                              nocovs = ling_loadings[, 1:9],
                              peercovs = ling_loadings[, 1:6],
                              justage = ling_loadings[, 1:8],
                              justdisorder = ling_loadings[, 1:8],
                              agedisorder = ling_loadings[, 1:7])
    }

    w_matrix <- switch(peer,
                       peercovs = w_matrix[, -c(1, 2, 3)],
                       justage = w_matrix[, -1],
                       justdisorder = w_matrix[, -1],
                       agedisorder = w_matrix[, -c(1, 2)],
                       w_matrix)

    # All pairs
    mat1_vec <- as.vector(scale(as.matrix(w_matrix)))
    mat2_vec <- as.vector(scale(as.matrix(ling_loadings)))
    differences <- abs(mat1_vec - mat2_vec)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    gene_factor_names <- expand.grid(
      gene = rownames(w_matrix),
      factor = colnames(w_matrix)
    )
    rownames(gene_factor_names) <- paste0(gene_factor_names$gene, "_", gene_factor_names$factor)

    df <- data.frame(
      Mat1 = mat1_vec,
      Mat2 = mat2_vec,
      distance = differences,
      category = ifelse(differences <= threshold, "Close to Line", "Far from Line"),
      row.names = rownames(gene_factor_names)
    )

    all_plots[[study]][[peer]] <- df

    # Only 4th factor
    peer4 <- as.vector(scale(w_matrix[, 4]))
    ling4 <- as.vector(scale(ling_loadings[, 4]))
    differences <- abs(peer4 - ling4)
    diff_sd <- sd(differences, na.rm = TRUE)
    threshold <- diff_sd

    df4 <- data.frame(
      Mat1 = peer4,
      Mat2 = ling4,
      distance = differences,
      category = ifelse(differences <= threshold, "Close to Line", "Far from Line"),
      row.names = rownames(w_matrix)
    )

    only4_plots[[study]][[peer]] <- df4

    gc()
  }
}

common_genes_per_peer <- list()

for (peer in peers) {
  genes_per_study <- list()

  for (study in studies) {
    df <- all_plots[[study]][[peer]]

    close_genes <- rownames(df)[df$category == "Close to Line"]

    genes_per_study[[study]] <- close_genes
  }

  # Get intersection across studies
  common_genes_per_peer[[peer]] <- Reduce(intersect, genes_per_study)
}

common_genes_per_peer4 <- list()

for (peer in peers) {
  genes_per_study <- list()

  for (study in studies) {
    df <- only4_plots[[study]][[peer]]

    close_genes <- rownames(df)[df$category == "Close to Line"]

    genes_per_study[[study]] <- close_genes
  }

  # Get intersection across studies
  common_genes_per_peer4[[peer]] <- Reduce(intersect, genes_per_study)
}

saveRDS(common_genes_per_peer, "results/Ling/common_genes_per_peer.rds")
saveRDS(common_genes_per_peer4, "results/Ling/common_genes_per_peer4.rds")

#### get same patients props
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Opc", "Oli")
szbd_patients <- c("SZ1", "SZ6", "SZ9", "SZ13", "SZ14", "SZ18", "SZ19", "SZ20", "SZ21", "SZ23", "SZ24")
ling_patients <- c("S13794", "S13450", "S02619", "S03190", "S18298", "S08116", "S07568", "S18557", "S07914", "S18978", "S14741")

patient_map <- setNames(ling_patients, szbd_patients)

ling_cell_meta_prop <- fread("data/data_processed/Ling/BA46.jointMetadata.txt") |>
  mutate(predClass = recode(predClass,
                            "glutamatergic" = "Exc",
                            "astrocyte" = "Ast",
                            "gabaergic" = "Inh",
                            "oligodendrocyte" = "Oli",
                            "microglia" = "Mic",
                            "polydendrocyte" = "Opc")) |>
  filter(DONOR %in% ling_patients) |>
  group_by(DONOR, predClass) |>
  count(name = "cell_count") |>
  group_by(DONOR) |>
  mutate(proportion = cell_count / sum(cell_count)) |>
  filter(predClass != "endothelia")

szbd_cell_meta_prop <- list()

for (cell_type in cell_types) {
  path <- paste0("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/SZBDMulti-Seq/0.Raw/", cell_type, "_SZBDMulti-Seq_SZ.rds")
  lists <- readRDS(path)
  meta <- lists$meta

  meta_filt <- meta |>
    filter(patientID %in% szbd_patients) |>
    group_by(patientID) |>
    summarise(cell_count = n(), .groups = "drop") |>
    mutate(predClass = cell_type)

  szbd_cell_meta_prop[[cell_type]] <- meta_filt
}

combined_meta <- bind_rows(szbd_cell_meta_prop)
szbd_cell_meta_prop_df <- combined_meta |>
  group_by(patientID) |>
  mutate(proportion = cell_count / sum(cell_count)) |>
  dplyr::select(patientID, predClass, cell_count, proportion)

szbd_renamed <- szbd_cell_meta_prop_df |>
  mutate(patientID = patient_map[patientID])

ling_clean <- ling_cell_meta_prop |>
  rename(patientID = DONOR) |>
  mutate(cohort = "Ling")

szbd_clean <- szbd_renamed |>
  mutate(cohort = "SZBD")

combined_prop <- bind_rows(ling_clean, szbd_clean)

per_prop <- ggplot(combined_prop, aes(x = predClass, y = proportion, fill = cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ patientID, ncol = 4) +
  labs(title = "Cell Type Proportions per Patient",
       x = "Cell Type",
       y = "Proportion",
       fill = "Cohort") +
  scale_fill_manual(values = c("navy","palevioletred"))+
  theme_bw()

ggsave(plot = per_prop, filename = "results/Ling/celltype_prop.png")

per_number <- ggplot(combined_prop, aes(x = predClass, y = cell_count, fill = cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ patientID, ncol = 4) +
  labs(title = "Cell Type Counts per Patient",
       x = "Cell Type",
       y = "Counts",
       fill = "Cohort") +
  scale_fill_manual(values = c("navy","palevioletred"))+
  theme_bw()
ggsave(plot = per_number, filename = "results/Ling/celltype_count.png")

cell_type_colors <- c(
  "Ast" = "turquoise",
  "Opc" = "coral",
  "Exc" = "pink2",
  "Inh" = "olivedrab3",
  "Mic" = "navy",
  "Oli" = "yellow2"
)

stack_prop <- ggplot(combined_prop, aes(x = cohort, y = proportion, fill = predClass)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_wrap(. ~ patientID, switch = "x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions per Patient in Each Study",
       y = "Proportion", x = NULL, fill = "Cell type") +
  theme_bw()+
  scale_fill_manual(values = cell_type_colors)
ggsave(plot = stack_prop, filename = "results/Ling/Samepatients/celltype_prop_stack.png", width = 8, height =8)

stack_count <- ggplot(combined_prop, aes(x = cohort, y = cell_count, fill = predClass)) +
  geom_bar(stat = "identity", position = position_stack()) +
  facet_wrap(. ~ patientID, switch = "x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Counts per Patient in Each Study",
       y = "Counts", x = NULL, fill = "Cell type") +
  theme_bw()+
  scale_fill_manual(values = cell_type_colors)
ggsave(plot = stack_count, filename = "results/Ling/celltype_count_stack.png", width = 8, height = 8)

continuous_plot <- ggplot(combined_prop, aes(x = predClass, y = cell_count, color = cohort)) +
  geom_point() +
  geom_line(aes(group = cohort)) +
  facet_wrap(. ~ patientID) +
  theme_classic() +
  labs(title = "Cell Type Counts per Patient in Each Study",
       y = "Counts", x = NULL, color = "Study")+
  scale_color_manual(values = c("navy","palevioletred"))
ggsave(plot = continuous_plot, filename = "results/Ling/Samepatients/celltype_count_stack_linear.png", width = 8, height = 8)

### correlation matrices
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Opc", "Oli")
szbd_patients <- c("SZ1", "SZ6", "SZ9", "SZ13", "SZ14", "SZ18", "SZ19", "SZ20", "SZ21", "SZ23", "SZ24")
ling_patients <- c("S13794", "S13450", "S02619", "S03190", "S18298", "S08116", "S07568", "S18557", "S07914", "S18978", "S14741")
patient_map <- setNames(ling_patients, szbd_patients)

first_file <- paste0("data/data_raw/SZBDMulti-Seq/", szbd_patients[1], "-annotated_matrix.txt")
first_data <- fread(first_file, header = TRUE, sep = "\t", data.table = FALSE)
gene_order <- first_data[[1]]
pb_list <- list()
for (szbd in szbd_patients){
  sample_file <- paste0("data/data_raw/SZBDMulti-Seq/", szbd, "-annotated_matrix.txt")
  sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
  rownames(sample_data) <- sample_data[, 1]
  sample_data <- sample_data[, -1]
  sample_data <- sample_data[gene_order, , drop = FALSE]
  sample_summed <- data.frame(rowSums(sample_data))
  colnames(sample_summed) <- patient_map[szbd]  # name the column after the sample ID

  pb_list[[szbd]] <- sample_summed  # add to the list

  gc()
}
combined_szbd <- do.call(cbind, pb_list)
cpm_szbd <- do_cpm_log(combined_szbd)
saveRDS(cpm_szbd, "data/data_processed/SZBDMulti-Seq/Pseudobulk/PseudobulkCPM/full_pseudobulk.rds")


pb_list <- list()
first_file <- fread("data/data_raw/Ling/1203113_S07914/features.tsv.gz") |> dplyr::select(V2)
gene_order <- first_file$V2

for (ling in ling_patients){
  top.dirs <- list.dirs("data/data_raw/Ling", full.names = TRUE, recursive = FALSE)
  matching_dir <- top.dirs[grepl(ling, top.dirs)]

  barcodes <- fread(paste0(matching_dir, "/barcodes.tsv.gz"), header = FALSE) |> pull()
  genes <- fread(paste0(matching_dir, "/features.tsv.gz")) |> dplyr::select(V2) |> pull()
  matrix_mex <- readMM(gzfile(paste0(matching_dir, "/matrix.mtx.gz")))

  rownames(matrix_mex) <- genes
  colnames(matrix_mex) <- barcodes
  matrix_mex <- matrix_mex[gene_order, , drop = FALSE]

  sample_summed <- data.frame(rowSums(matrix_mex))
  colnames(sample_summed) <- ling

  pb_list[[ling]] <- sample_summed
  gc()
}
combined_ling <- do.call(cbind, pb_list)
cpm_ling <- do_cpm_log(combined_ling)
saveRDS(cpm_ling, "data/data_processed/Ling/Pseudobulk/PseudobulkCPM/full_pseudobulk.rds")

### correlation patient vs patient
szbd_pb <- readRDS("data/data_processed/SZBDMulti-Seq/Pseudobulk/PseudobulkCPM/full_pseudobulk.rds")
ling_pb <- readRDS("data/data_processed/Ling/Pseudobulk/PseudobulkCPM/full_pseudobulk.rds")

common_genes <- intersect(rownames(szbd_pb), rownames(ling_pb))
szbd_pb <- szbd_pb[rownames(szbd_pb) %in% common_genes, ]
ling_pb <- ling_pb[rownames(ling_pb) %in% common_genes, ]
ling_pb <- ling_pb[rownames(szbd_pb), , drop = FALSE]
colnames(szbd_pb) <-szbd_patients

cor_mat <- cor(ling_pb, szbd_pb, method = "pearson")
heatmap_all <- pheatmap(cor_mat,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    display_numbers = TRUE,
                    main = "Correlation between Ling and SZBD profiles")
ggsave(heatmap_all, filename = "results/Ling/Samepatients/heatmap_both.png", width = 5, height = 5)

### correlation see pairs
cor_mat_ling <- cor(ling_pb, ling_pb, method = "pearson")
heatmap_ling <- pheatmap(cor_mat_ling,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        display_numbers = TRUE,
                        main = "Ling vs Ling")

ggsave(heatmap_ling, filename = "results/Ling/Samepatients/heatmap_ling.png", width = 5, height = 5)

cor_mat_szbd <- cor(szbd_pb, szbd_pb, method = "pearson")
heatmap_szbd <- pheatmap(cor_mat_szbd,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         display_numbers = TRUE,
                         main = "SZBD vs SZBD")

ggsave(heatmap_szbd, filename = "results/Ling/Samepatients/heatmap_szbd.png", width = 5, height = 5)

checked_patients <- list()
for (ling_patient in ling_patients){
  ling_correlations <- cor_mat_ling[ling_patient, ]
  ling_correlations <- ling_correlations[names(ling_correlations) != ling_patient]

  ling_top_match <- names(which.max(ling_correlations))

  szbd_patient <- names(patient_map)[patient_map == ling_patient]
  szbd_top_match_target <- names(patient_map)[patient_map == ling_top_match]

  if (length(szbd_patient) == 0 || length(szbd_top_match_target) == 0) {
    warning(paste("Missing mapping for", ling_patient, "or", ling_top_match))
    next
  }

 szbd_correlations <- cor_mat_szbd[szbd_patient, ]
  szbd_correlations <- szbd_correlations[names(szbd_correlations) != szbd_patient]
  szbd_top_match <- names(which.max(szbd_correlations))

  is_same_match <- (szbd_top_match == szbd_top_match_target)

  checked_patients[[ling_patient]] <- list(
    ling_match = ling_top_match,
    szbd_patient = szbd_patient,
    szbd_match = szbd_top_match,
    expected_szbd_match = szbd_top_match_target,
    match_same = is_same_match
  )
}

checked_df <- do.call(rbind, lapply(names(checked_patients), function(p) {
  cbind(patient = p, as.data.frame(checked_patients[[p]], stringsAsFactors = FALSE)
  )}))

saveRDS(checked_df, "results/Ling/Samepatients/check_matches.rds")

### correlation add random patient, see correlation
for (ling_patient in ling_patients){
  szbd_patient <- names(patient_map)[patient_map == ling_patient]

  col_szbd <- szbd_pb[, colnames(szbd_pb) == szbd_patient, drop = FALSE]
  new_pb <- cbind(ling_pb, col_szbd)

  cor_mat <- cor(new_pb, new_pb, method = "pearson")
  heatmap <- pheatmap(cor_mat,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           display_numbers = TRUE,
                           main = paste0(ling_patient, "-", szbd_patient))

  ggsave(heatmap, filename = paste0("results/Ling/Samepatients/ling_cors/", ling_patient, "_add.png"), width = 5, height = 5)
}

for (szbd_patient in szbd_patients){
  ling_patient <- patient_map[szbd_patient]

  col_ling <- ling_pb[, colnames(ling_pb) == ling_patient, drop = FALSE]
  new_pb <- cbind(szbd_pb, col_ling)

  cor_mat <- cor(new_pb, new_pb, method = "pearson")
  heatmap <- pheatmap(cor_mat,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      display_numbers = TRUE,
                      main = paste0(szbd_patient, "-", ling_patient))

  ggsave(heatmap, filename = paste0("results/Ling/Samepatients/szbd_cors/", szbd_patient, "_add.png"), width = 5, height = 5)
}

### plot expr against expr
matches <- readRDS("results/Ling/Samepatients/check_matches.rds")

for (ling_patient in ling_patients){
  match <- matches |>
    filter(patient == ling_patient)

  match_ling <- match |> dplyr::select(ling_match) |> pull()

  match_szbd <- match |> dplyr::select(szbd_patient) |> pull()

  patient_col <- ling_pb[, colnames(ling_pb) == ling_patient, drop = FALSE]
  match_col <- ling_pb[, colnames(ling_pb) == match_ling, drop = FALSE]
  match_szbd_col <- szbd_pb[, colnames(szbd_pb) == match_szbd, drop = FALSE]

  all_3 <- data.frame(
    ling_og = as.vector(patient_col),
    ling_match = as.vector(match_col),
    szbd_match = as.vector(match_szbd_col)
  )

  ling_ling <- ggplot(all_3, aes(x = ling_og, y = ling_match)) +
    geom_point(alpha = 0.3) +
    labs(x = ling_patient,
         y = match_ling) +
    theme_classic()+
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")

  ling_szbd <- ggplot(all_3, aes(x = ling_og, y = szbd_match)) +
    geom_point(alpha = 0.3) +
    labs(x = ling_patient,
         y = match_szbd) +
    theme_classic()+
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")

  both_plots <- plot_grid(plotlist = list(ling_ling, ling_szbd))

  ggsave(both_plots, filename = paste0("results/Ling/Samepatients/ling_cors/", ling_patient, "_scatter.png"), width = 10, height = 5)
}

### CORRELATIONS FOR EXC NEURONS
# wrangle pbs
szbd_cpm <- readRDS("data/data_processed/SZBDMulti-Seq/Pseudobulk/PseudobulkCPM/Exc_SZBDMulti-Seq_SZ.rds")
szbd_expr <- szbd_cpm$expr
clean_names <- sub(".*(SZ\\d+):.*", "\\1", colnames(szbd_expr))
colnames(szbd_expr) <- clean_names
szbd_filt <- szbd_expr[, colnames(szbd_expr) %in% szbd_patients, drop = FALSE]
szbd_filt <- szbd_filt[, szbd_patients, drop = FALSE]

ling_pb <- readRDS("data/data_processed/Ling/Pseudobulk/PseudobulkRaw/Exc_Ling_SZ.rds")
ling_expr <- ling_pb$expr
ling_filt <- ling_expr[, colnames(ling_expr) %in% ling_patients, drop = FALSE]
ling_filt <- ling_filt[, ling_patients, drop =FALSE]
ling_cpm <- do_cpm_log(ling_filt)

szbd_genes <- rownames(szbd_filt)
ling_genes <- rownames(ling_cpm)
common_genes <- intersect(szbd_genes, ling_genes)

szbd_gene_common <- szbd_filt[rownames(szbd_filt) %in% common_genes, , drop = FALSE]
szbd_gene_ordered <- szbd_gene_common[common_genes, , drop = FALSE]

ling_gene_common <- ling_cpm[rownames(ling_cpm) %in% common_genes, , drop = FALSE]
ling_gene_ordered <- ling_gene_common[common_genes, , drop = FALSE]
ling_gene_ordered <- as.matrix(ling_gene_ordered)

saveRDS(szbd_gene_ordered, "data/data_processed/SZBDMulti-Seq/Pseudobulk/PseudobulkCPM/exc_filt_pseudobulk.rds")
saveRDS(ling_gene_ordered, "data/data_processed/Ling/Pseudobulk/PseudobulkCPM/exc_filt_pseudobulk.rds")

# all
szbd_pb <- readRDS("data/data_processed/SZBDMulti-Seq/Pseudobulk/PseudobulkCPM/exc_filt_pseudobulk.rds")
ling_pb <- readRDS("data/data_processed/Ling/Pseudobulk/PseudobulkCPM/exc_filt_pseudobulk.rds")

cor_mat <- cor(ling_pb, szbd_pb, method = "pearson")
heatmap_all <- pheatmap(cor_mat,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        display_numbers = TRUE,
                        main = "Ling vs SZBD (Exc)")
ggsave(heatmap_all, filename = "results/Ling/Samepatients/heatmap_both_exc.png", width = 5, height = 5)

# correlation see pairs
cor_mat_ling <- cor(ling_pb, ling_pb, method = "pearson")
heatmap_ling <- pheatmap(cor_mat_ling,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         display_numbers = TRUE,
                         main = "Ling vs Ling (Exc)")

ggsave(heatmap_ling, filename = "results/Ling/Samepatients/heatmap_ling_exc.png", width = 5, height = 5)

cor_mat_szbd <- cor(szbd_pb, szbd_pb, method = "pearson")
heatmap_szbd <- pheatmap(cor_mat_szbd,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         display_numbers = TRUE,
                         main = "SZBD vs SZBD (Exc)")

ggsave(heatmap_szbd, filename = "results/Ling/Samepatients/heatmap_szbd_exc.png", width = 5, height = 5)

checked_patients <- list()
for (ling_patient in ling_patients){
  ling_correlations <- cor_mat_ling[ling_patient, ]
  ling_correlations <- ling_correlations[names(ling_correlations) != ling_patient]

  ling_top_match <- names(which.max(ling_correlations))

  szbd_patient <- names(patient_map)[patient_map == ling_patient]
  szbd_top_match_target <- names(patient_map)[patient_map == ling_top_match]

  if (length(szbd_patient) == 0 || length(szbd_top_match_target) == 0) {
    warning(paste("Missing mapping for", ling_patient, "or", ling_top_match))
    next
  }

  szbd_correlations <- cor_mat_szbd[szbd_patient, ]
  szbd_correlations <- szbd_correlations[names(szbd_correlations) != szbd_patient]
  szbd_top_match <- names(which.max(szbd_correlations))

  is_same_match <- (szbd_top_match == szbd_top_match_target)

  checked_patients[[ling_patient]] <- list(
    ling_match = ling_top_match,
    szbd_patient = szbd_patient,
    szbd_match = szbd_top_match,
    expected_szbd_match = szbd_top_match_target,
    match_same = is_same_match
  )
}

checked_df <- do.call(rbind, lapply(names(checked_patients), function(p) {
  cbind(patient = p, as.data.frame(checked_patients[[p]], stringsAsFactors = FALSE)
  )}))

saveRDS(checked_df, "results/Ling/Samepatients/check_matches_exc.rds")

# correlation add random patient, see correlation
for (ling_patient in ling_patients){
  szbd_patient <- names(patient_map)[patient_map == ling_patient]

  col_szbd <- szbd_pb[, colnames(szbd_pb) == szbd_patient, drop = FALSE]
  new_pb <- cbind(ling_pb, col_szbd)

  cor_mat <- cor(new_pb, new_pb, method = "pearson")
  heatmap <- pheatmap(cor_mat,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      display_numbers = TRUE,
                      main = paste0(ling_patient, "-", szbd_patient, " (Exc)"))

  ggsave(heatmap, filename = paste0("results/Ling/Samepatients/ling_cors/", ling_patient, "_add_exc.png"), width = 5, height = 5)
}

for (szbd_patient in szbd_patients){
  ling_patient <- patient_map[szbd_patient]

  col_ling <- ling_pb[, colnames(ling_pb) == ling_patient, drop = FALSE]
  new_pb <- cbind(szbd_pb, col_ling)

  cor_mat <- cor(new_pb, new_pb, method = "pearson")
  heatmap <- pheatmap(cor_mat,
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      display_numbers = TRUE,
                      main = paste0(szbd_patient, "-", ling_patient, " (Exc)"))

  ggsave(heatmap, filename = paste0("results/Ling/Samepatients/szbd_cors/", szbd_patient, "_add_exc.png"), width = 5, height = 5)
}

# scatter plots
matches <- readRDS("results/Ling/Samepatients/check_matches_exc.rds")

for (ling_patient in ling_patients){
  match <- matches |>
    filter(patient == ling_patient)

  match_ling <- match |> dplyr::select(ling_match) |> pull()

  match_szbd <- match |> dplyr::select(szbd_patient) |> pull()

  patient_col <- ling_pb[, colnames(ling_pb) == ling_patient, drop = FALSE]
  match_col <- ling_pb[, colnames(ling_pb) == match_ling, drop = FALSE]
  match_szbd_col <- szbd_pb[, colnames(szbd_pb) == match_szbd, drop = FALSE]

  all_3 <- data.frame(
    ling_og = as.vector(patient_col),
    ling_match = as.vector(match_col),
    szbd_match = as.vector(match_szbd_col)
  )

  ling_ling <- ggplot(all_3, aes(x = ling_og, y = ling_match)) +
    geom_point(alpha = 0.3) +
    labs(x = ling_patient,
         y = match_ling) +
    theme_classic()+
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")

  ling_szbd <- ggplot(all_3, aes(x = ling_og, y = szbd_match)) +
    geom_point(alpha = 0.3) +
    labs(x = ling_patient,
         y = match_szbd) +
    theme_classic()+
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")

  both_plots <- plot_grid(plotlist = list(ling_ling, ling_szbd))

  ggsave(both_plots, filename = paste0("results/Ling/Samepatients/ling_cors/", ling_patient, "_scatter_exc.png"), width = 10, height = 5)
}

### lf4 against lf4
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Opc", "Oli")
szbd_patients <- c("SZ1", "SZ6", "SZ9", "SZ13", "SZ14", "SZ18", "SZ19", "SZ20", "SZ21", "SZ23", "SZ24")
ling_patients <- c("S13794", "S13450", "S02619", "S03190", "S18298", "S08116", "S07568", "S18557", "S07914", "S18978", "S14741")
patient_map <- setNames(ling_patients, szbd_patients)
studies <- c("Ling", "SZBDMulti-Seq")
merged_df <- list()
for (study in studies){
  peer_data <- fread(paste0("results/Ling/", study, "/nocovs/X.csv"))
  matrix <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
  colnames(peer_data) <- rownames(matrix)
  meta <- readRDS(paste0("data/data_processed/", study, "/", study, "-patient.rds"))
  peer_data <- t(peer_data)
  peer_data <- as.data.frame(peer_data)
  peer_factors <- paste0("PEER", 1:10)
  colnames(peer_data) <- peer_factors

  if (study != "Ling"){
    clean_peer_ids <- sub("^patientID", "", rownames(peer_data))
    clean_peer_ids <- sub(":.*", "", clean_peer_ids)
  rownames(peer_data) <- clean_peer_ids
  }
  peer_data$patientID <- rownames(peer_data)
  merged <- merge(peer_data, meta, by = "patientID")
  merged_df[[study]] <- merged
}

SZBD_merged <- merged_df[["SZBDMulti-Seq"]]
SZBD_patients <- SZBD_merged[SZBD_merged$patientID %in% szbd_patients,]
Ling_merged <- merged_df[["Ling"]]
Ling_patients <- Ling_merged[Ling_merged$patientID %in% ling_patients, ]
SZBD_patients$Ling_match <- patient_map[SZBD_patients$patientID]

full_merged <- merge(SZBD_patients, Ling_patients, by.x = "Ling_match", by.y = "patientID")
peer_factors <- paste0("PEER", 1:10)
scatter_list <- list()
for (PEER in peer_factors){
  scatterplot <- ggplot(full_merged) +
    geom_point(aes(x = .data[[paste0(PEER, ".x")]],
                   y = .data[[paste0(PEER, ".y")]])) +
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
    theme_classic() +
    labs(x = "SZBD", y = "Ling", title = PEER)
scatter_list[[PEER]] <- scatterplot
}
all_scatters <- plot_grid(plotlist = scatter_list,ncol = 5)
ggsave(all_scatters, filename = "results/Ling/Samepatients/scatters_LFs.png", width = 10, height = 5)
### get gemma DEA
# set_gemma_user(username = "xxxxx", password = "xxxxxxxxx")
sc_experiment_list <- c("Ling-2024", "CMC", "SZBDMulti-Seq")
dea_res <- list()
for (experiment in sc_experiment_list){
  sc_DE <- get_differential_expression_values(experiment)
  dea_res[[experiment]] <- sc_DE
}

dea_meta <- list()
for (experiment in sc_experiment_list){
  sc_meta <- get_dataset_differential_expression_analyses(experiment)
  dea_meta[[experiment]] <- sc_meta
}

saveRDS(dea_res, "results/DEA/gemma_dea.rds")
saveRDS(dea_meta, "results/DEA/gemma_dea_meta.rds")

### rename things
age_list <- list()
for (experiment in sc_experiment_list){
  age_meta <- dea_meta[[experiment]] |>
    filter(factor.category == "age") |>
    dplyr::select(result.ID, subsetFactor) |>
    mutate(cell_type = map_chr(subsetFactor, ~ .x$value)) |>
    dplyr::select(-subsetFactor) |>
    unique()

  ids <- age_meta$result.ID

  for (i in seq_len(nrow(age_meta))) {
    id <- age_meta$result.ID[i]
    id <- as.character(id)
    cell_type <- age_meta$cell_type[i]

    dea_subset <- dea_res[[experiment]][[id]]

    age_list[[experiment]][[cell_type]] <- dea_subset
  }

}
saveRDS(age_list, "results/DEA/age_gemma.rds")

disease_list <- list()
for (experiment in sc_experiment_list){
  disease_meta <- dea_meta[[experiment]] |>
    filter(factor.category == "disease") |>
    dplyr::select(result.ID, subsetFactor) |>
    mutate(cell_type = map_chr(subsetFactor, ~ .x$value)) |>
    dplyr::select(-subsetFactor) |>
    unique()

  ids <- disease_meta$result.ID

  for (i in seq_len(nrow(disease_meta))) {
    id <- disease_meta$result.ID[i]
    id <- as.character(id)
    cell_type <- disease_meta$cell_type[i]

    dea_subset <- dea_res[[experiment]][[id]]

    disease_list[[experiment]][[cell_type]] <- dea_subset
  }
}
saveRDS(disease_list, "results/DEA/disease_gemma.rds")

### get down and up regulated genes for each study
age_list <- readRDS("results/DEA/age_gemma.rds")
disease_list <- readRDS("results/DEA/disease_gemma.rds")

genes_list <- list()

for (experiment in sc_experiment_list){
  age_deas <- age_list[[experiment]]
  disease_deas <- disease_list[[experiment]]

  ast_age <- age_deas[["astrocyte"]]
  ast_age_genes <- ast_age |>
    filter(contrast_pvalue <= 0.05 & contrast_log2fc >= 0.05)
  ast_age_up <- ast_age_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  ast_age_down <- ast_age_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  astrocytes <- list(upregulated = ast_age_up, downregulated = ast_age_down)
  genes_list[[experiment]][["age"]][["astrocytes"]] <- astrocytes

  exc_age_list <- age_deas[grep("glutamatergic", names(age_deas))]
  exc_age <- do.call(rbind, exc_age_list)
  exc_age_genes <- exc_age |>
    filter(contrast_pvalue <= 0.05& contrast_log2fc >= 0.05)
  exc_age_up <- exc_age_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  exc_age_down <- exc_age_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  excitatory <- list(upregulated = exc_age_up, downregulated = exc_age_down)
  genes_list[[experiment]][["age"]][["excitatory"]] <- excitatory

  inh_age_list <- age_deas[grep("GABAergic", names(age_deas))]
  inh_age <- do.call(rbind, inh_age_list)
  inh_age_genes <- inh_age |>
    filter(contrast_pvalue <= 0.05& contrast_log2fc >= 0.05)
  inh_age_up <- inh_age_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  inh_age_down <- inh_age_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  inhibitory <- list(upregulated = inh_age_up, downregulated = inh_age_down)
  genes_list[[experiment]][["age"]][["inhibitory"]] <- inhibitory

  if (experiment == "Ling-2024"){
    corrected_pvalue <- "contrast_301802_pvalue"
    contrast_log2fc <- "contrast_301802_log2fc"
  }

  if (experiment == "CMC"){
    corrected_pvalue <- "contrast_301166_pvalue"
    contrast_log2fc <- "contrast_301166_log2fc"
  }
  if (experiment == "SZBDMulti-Seq"){
    corrected_pvalue <- "contrast_301372_pvalue"
    contrast_log2fc <- "contrast_301372_log2fc"
  }

  ast_disease <- disease_deas[["astrocyte"]]
  ast_disease_genes <- ast_disease |>
    filter(corrected_pvalue <= 0.05& contrast_log2fc >= 0.05)
  ast_disease_up <- ast_disease_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  ast_disease_down <- ast_disease_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  astrocytes <- list(upregulated = ast_disease_up, downregulated = ast_disease_down)
  genes_list[[experiment]][["disease"]][["astrocytes"]] <- astrocytes

  exc_disease_list <- disease_deas[grep("glutamatergic", names(disease_deas))]
  exc_disease <- do.call(rbind, exc_disease_list)
  exc_disease_genes <- exc_disease |>
    filter(corrected_pvalue <= 0.05& contrast_log2fc >= 0.05)
  exc_disease_up <- exc_disease_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  exc_disease_down <- exc_disease_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  excitatory <- list(upregulated = exc_disease_up, downregulated = exc_disease_down)
  genes_list[[experiment]][["disease"]][["excitatory"]] <- excitatory

  inh_disease_list <- disease_deas[grep("GABAergic", names(disease_deas))]
  inh_disease <- do.call(rbind, inh_disease_list)
  inh_disease_genes <- inh_disease |>
    filter(corrected_pvalue <= 0.05& contrast_log2fc >= 0.05)
  inh_disease_up <- inh_disease_genes |>
    filter(contrast_log2fc > 0) |>
    pull(GeneSymbol)
  inh_disease_down <- inh_disease_genes |>
    filter(contrast_log2fc < 0) |>
    pull(GeneSymbol)
  inhibitory <- list(upregulated = inh_disease_up, downregulated = inh_disease_down)
  genes_list[[experiment]][["disease"]][["inhibitory"]] <- inhibitory
}
saveRDS(genes_list, "results/DEA/gemma_DEGs.rds")

### rug plots
studies <- c("CMC", "Ling-2024", "SZBDMulti-Seq")

ll<-read.table("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt", header=T, sep='\t')

ll.ex<-ll[ grepl(x=ll$id,pattern=".+_glutamatergic") ,]
ll.ex[,1]<-gsub(ll.ex[,1], pattern="_glutamatergic", r="")
ll.ex.lf4<-ll.ex[,c(1,5)] # isolate factor 4
row.names(ll.ex.lf4)<-ll.ex.lf4[,1]
ll.ex.lf4<-ll.ex.lf4[,2, drop=F]

ll.in<-ll[ grepl(x=ll$id,pattern=".+_gabaergic") ,]
ll.in[,1]<-gsub(ll.in[,1], pattern="_gabaergic", r="")
ll.in.lf4<-ll.in[,c(1,5)] # isolate factor 4
row.names(ll.in.lf4)<-ll.in.lf4[,1]
ll.in.lf4<-ll.in.lf4[,2, drop=F]

ll.as<-ll[grepl(x=ll$id,pattern=".+_astrocyte"),]
ll.as[,1]<-gsub(ll.as[,1], pattern="_astrocyte", r="")
ll.as.lf4<-ll.as[,c(1,5)]
row.names(ll.as.lf4)<-ll.as.lf4[,1]
ll.as.lf4<-ll.as.lf4[,2, drop=F]

cell_types <- c("excitatory", "inhibitory", "astrocytes")
contrasts <- c("age", "disease")

for (study in studies){
  for (contrast in contrasts){
    all_plots <- list()
  for (cell_type in cell_types){

    deg_up <- genes_list[[study]][[contrast]][[cell_type]][["upregulated"]]
    deg_down <- genes_list[[study]][[contrast]][[cell_type]][["downregulated"]]

    if (cell_type == "astrocytes"){
      p_up <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Astrocytes",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.as.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.as.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Astrocytes",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.as.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }

    if (cell_type == "excitatory"){
      p_up <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Ex. Neurons",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.ex.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.ex.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="Ex. Neurons",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.ex.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }

    if (cell_type == "inhibitory"){
      p_up <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="In. Neurons",
             subtitle = "Rug: top up-regulated genes")

      p_up <-p_up + geom_rug(data=ll.in.lf4[deg_up,,drop=F], aes(x=PEER_4), sides='b')

      p_down <-ggplot(ll.in.lf4, aes(x=PEER_4)) +
        geom_histogram(bins=100) +
        theme_classic() +
        labs(title="In. Neurons",
             subtitle = "Rug: top down-regulated genes")

      p_down <-p_down + geom_rug(data=ll.in.lf4[deg_down,,drop=F], aes(x=PEER_4), sides='b')
    }
    all_plots[[paste0(cell_type, "_up")]] <- p_up
    all_plots[[paste0(cell_type, "_down")]] <- p_down

  }
    plot_list <- plot_grid(plotlist = all_plots, ncol = 2, align = "hv")

    title_plot <- ggdraw() +
      draw_label(paste0(study, ": ", contrast, " contrast"), x = 0.5, y = 0.5, size = 12, hjust = 0.5)
    all_plots_with_title <- plot_grid(title_plot, plot_list, ncol = 1, rel_heights = c(0.1, 1))

    final_path <- paste0("results/Ling/", study, "/rug_", contrast, ".png")
    ggsave(plot = all_plots_with_title, filename = final_path, height = 12, width = 10)
  }}

### LFs plots
studies <- c("Batiuk", "CMC", "MultiomeBrain", "Ling", "SZBDMulti-Seq")
peers <- c("nocovs")

for (study in studies){
  for (peer in peers){
    peer_data <- fread(paste0("results/Ling/", study, "/", peer, "/X.csv"))
    matrix <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
    colnames(peer_data) <- rownames(matrix)
    meta <- readRDS(paste0("data/data_processed/", study, "/", study, "-patient.rds"))
    peer_data <- t(peer_data)
    peer_data <- as.data.frame(peer_data)

    plotlist <- list()
    for (i in seq_along(colnames(peer_data))){
      if (study != "Ling"){
        clean_peer_ids <- sub("^patientID", "", rownames(peer_data))
        clean_peer_ids <- sub(":.*", "", clean_peer_ids)
        rownames(peer_data) <- clean_peer_ids
      }
      df <- data.frame(
        peer = peer_data[,i],
        patientID = rownames(peer_data)
      )

      merged <- merge(df, meta, by = "patientID")
      colnames(merged) <- make.unique(colnames(merged))
      merged <- merged |>
        mutate(age = as.numeric(str_replace(age, "\\+", "")))

      scatterplot <- ggplot(merged, aes(x = age, y = peer, color = disorder)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.3, aes(fill=disorder)) +
        scale_x_continuous(breaks = seq(0, 100, by = 10)) +
        theme_classic()+
        scale_fill_manual(values = c("black","red"))+
        scale_color_manual(values = c("black","red"))+
        labs(x = "Age (years)",
             y = paste0("PEER", i),
             color = "Schizophrenia",
             fill = "Schizophrenia")+
        theme(legend.position = "none")

      boxplot <- ggplot(merged, aes(x = disorder, y = peer)) +
        geom_violin(trim = FALSE, aes(color = disorder)) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "royalblue4", fill = NA) +
        geom_jitter(width = 0.15, alpha = 0.2, size = 1, aes(color = disorder)) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = paste0("PEER", i), x = "") +
        theme(
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
        scale_fill_manual(values = c("black","red"))+
        scale_color_manual(values = c("black","red"))

      pair_plot <- cowplot::plot_grid(scatterplot, boxplot, ncol = 1, align = "v")
      plotlist[[length(plotlist) + 1]] <- pair_plot
      gc()
    }
    all_plots <- cowplot::plot_grid(plotlist = plotlist, ncol = 5, align = "hv")
    cowplot::save_plot(all_plots, filename = paste0("results/Ling/", study, "/", peer, "/LFexpression.png"), base_height = 9, base_width = 12)
    gc()
  }
}

### schizophrenia association with factors
studies <- c("Batiuk", "CMC", "Ling", "SZBDMulti-Seq", "MultiomeBrain")

plot_list <- list()
for (study in studies){
  peer_data <- fread(paste0("results/Ling/", study, "/nocovs/X.csv"))
  matrix <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
  colnames(peer_data) <- rownames(matrix)
  meta <- readRDS(paste0("data/data_processed/", study, "/", study, "-patient.rds"))
  peer_data <- t(peer_data)
  peer_data <- as.data.frame(peer_data)
  peer_factors <- paste0("PEER", 1:10)
  if (study == "MultiomeBrain"){
    peer_factors <- paste0("PEER", 1:9)
  }
  colnames(peer_data) <- peer_factors

  if (study != "Ling"){
    clean_peer_ids <- sub("^patientID", "", rownames(peer_data))
    clean_peer_ids <- sub(":.*", "", clean_peer_ids)
    rownames(peer_data) <- clean_peer_ids
  }
  peer_data$patientID <- rownames(peer_data)
  merged <- merge(peer_data, meta, by = "patientID")

  merged <- merged |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  merged$disorder <- factor(merged$disorder)
  merged$sex <- factor(merged$sex)
  merged$disorder <- relevel(merged$disorder, ref = "no")

  pvals <- numeric(length(peer_factors))

  for (i in seq_along(peer_factors)) {
    formula <- as.formula(paste(peer_factors[i], "~ disorder + age"))
    model <- lm(formula, data = merged)
    pvals[i] <- summary(model)$coefficients["disorderyes", "Pr(>|t|)"]
  }

  observed <- -log10(sort(pvals))
  expected <- -log10(ppoints(length(pvals)))

  qq_df <- data.frame(
    expected = expected,
    observed = observed,
    factor = seq_along(peer_factors)[order(pvals)]
  )

  qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_text(aes(label = factor), vjust = -0.5, size = 4) +
    theme_classic() +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = study
    )

  plot_list[[study]] <- qq_plot
}

ling_loadings <- read_xlsx("data/data_raw/41586_2024_7109_MOESM6_ESM.xlsx", sheet = 1)
ling_meta <- read_xlsx("data/data_processed/Ling/41586_2024_7109_MOESM5_ESM.xlsx")
merged <- merge(ling_loadings, ling_meta, by = "Donor")

merged <- merged |>
  mutate(age = as.numeric(str_replace(Age, "\\+", "")))
merged$Schizophrenia <- factor(merged$Schizophrenia)
merged$Sex <- factor(merged$Sex)
merged$Schizophrenia <- relevel(merged$Schizophrenia, ref = "Unaffected")

peer_factors <- paste0("PEER_", 1:10)
pvals <- numeric(length(peer_factors))

for (i in seq_along(peer_factors)) {
  formula <- as.formula(paste(peer_factors[i], "~ Schizophrenia+Age"))
  model <- lm(formula, data = merged)
  pvals[i] <- summary(model)$coefficients["SchizophreniaAffected", "Pr(>|t|)"]
}

observed <- -log10(sort(pvals))
expected <- -log10(ppoints(length(pvals)))

qq_df <- data.frame(
  expected = expected,
  observed = observed,
  factor = seq_along(peer_factors)[order(pvals)]
)

qq_plot <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_text(aes(label = factor), vjust = -0.5, size = 4) +
  theme_classic() +
  labs(
    x = "Expected -log10(p)",
    y = "Observed -log10(p)",
    title = "Ling by Ling")

plot_list[["Ling by Ling"]] <- qq_plot

all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv", greedy = FALSE)
ggsave(all_plots, filename = "results/Ling/qq_plots.png", width = 12, height = 10)

### get which cell types have the biggest influence
studies <- c("Batiuk", "CMC", "Ling", "SZBDMulti-Seq", "MultiomeBrain")

plot_list <- list()
for (study in studies){
  cpm_pb <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
  gene_names <- colnames(cpm_pb)

  if (study != "Ling"){
    gene_names <- sub("^[^.]+\\.", "", gene_names)}
  gene_names <- gene_names[-1]

  w_matrix <- fread(paste0("results/Ling/", study, "/nocovs/W.csv")) |> t() |> data.frame()
  rownames(w_matrix) <- gene_names

  w_matrix4 <- w_matrix[,4]
  abs_loadings <- abs(w_matrix4)
  top_indices <- order(abs_loadings, decreasing = TRUE)[1:1000]
  top_traits <- rownames(w_matrix)[top_indices]
  top_df <- data.frame(
    trait = rownames(w_matrix)[top_indices],
    LF4_loading = w_matrix4[top_indices]
  )
  cell_types <- sapply(strsplit(top_df$trait, "_"), function(x) x[2])
  cell_type_counts <- table(cell_types)
  cell_type_counts <- as.data.frame(cell_type_counts)

  cell_types_present <- unique(cell_types)
  cell_type_colors <- c(
    "Ast" = "turquoise",
    "Opc" = "coral",
    "Exc" = "pink2",
    "Inh" = "olivedrab3",
    "Mic" = "navy",
    "Oli" = "yellow2",
    "Gli" = "slategray2"
  )
  cell_type_colors_subset <- cell_type_colors[cell_types_present]
  cell_type_counts$cell_types <- factor(cell_type_counts$cell_types, levels = cell_type_counts$cell_types[order(-cell_type_counts$Freq)])
  hist <- ggplot(cell_type_counts)+
    geom_bar(aes(x = cell_types, y = Freq, fill = cell_types), stat = "identity")+
    labs(x = NULL, y = "Count", title = study) +
    scale_fill_manual(values = cell_type_colors_subset) +
    theme_classic()+
    theme(legend.position = "none")+
    ylim(c(0,1000))

  plot_list[[study]] <- hist
}

ling_loadings <- read.table("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt", header=T, sep='\t') |>
  dplyr::select(id, PEER_4)
ling_loadings$PEER_4 <- abs(ling_loadings$PEER_4)
ling_loadings <- ling_loadings[order(ling_loadings$PEER_4, decreasing = TRUE), ]
ling_loadings <- ling_loadings[1:1000, ]
cell_types <- sapply(strsplit(ling_loadings$id, "_"), function(x) x[2])
cell_type_counts <- table(cell_types)
cell_type_counts <- as.data.frame(cell_type_counts)
cell_type_counts <- cell_type_counts |>
  mutate(cell_types = as.character(cell_types)) |>
  mutate(cell_types = ifelse(cell_types == "gabaergic", "Inh", cell_types)) |>
  mutate(cell_types = ifelse(cell_types == "glutamatergic", "Exc", cell_types)) |>
  mutate(cell_types = ifelse(cell_types == "oligodendrocyte", "Oli", cell_types)) |>
  mutate(cell_types = ifelse(cell_types == "polydendrocyte", "Opc", cell_types)) |>
  mutate(cell_types = ifelse(cell_types == "astrocyte", "Ast", cell_types))

cell_types_present <- unique(cell_type_counts$cell_types)
cell_type_colors <- c(
  "Ast" = "turquoise",
  "Opc" = "coral",
  "Exc" = "pink2",
  "Inh" = "olivedrab3",
  "Mic" = "navy",
  "Oli" = "yellow2",
  "Gli" = "slategray2"
)
cell_type_colors_subset <- cell_type_colors[cell_types_present]

cell_type_counts$cell_types <- factor(cell_type_counts$cell_types, levels = cell_type_counts$cell_types[order(-cell_type_counts$Freq)])
hist <- ggplot(cell_type_counts)+
  geom_bar(aes(x = cell_types, y = Freq, fill = cell_types), stat = "identity")+
  labs(x = NULL, y = "Count", title = "Ling (by Ling)") +
  scale_fill_manual(values = cell_type_colors_subset) +
  theme_classic()+
  theme(legend.position = "none")+
  ylim(c(0,1000))

plot_list[["lingbyling"]] <- hist

all_plots <- plot_grid(plotlist = plot_list, ncol = 3)
ggsave(all_plots, filename = "results/Ling/cell_typesLF4.png")

### get cell type heatmap fig1g
studies <- c("CMC", "Ling", "SZBDMulti-Seq", "MultiomeBrain")
for (study in studies){
  cpm_pb <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
  gene_names <- colnames(cpm_pb)

  if (study != "Ling"){
    gene_names <- sub("^[^.]+\\.", "", gene_names)}
  gene_names <- gene_names[-1]

  w_matrix <- fread(paste0("results/Ling/", study, "/nocovs/W.csv")) |> t() |> data.frame()
  rownames(w_matrix) <- gene_names

  peer_list <- list()
  required_cell_types <- c("Exc", "Inh", "Ast", "Mic", "Opc", "Oli")

  for (i in seq_along(colnames(w_matrix))){
    w_matrix4 <- w_matrix[, i]
    abs_loadings <- abs(w_matrix4)
    top_indices <- order(abs_loadings, decreasing = TRUE)[1:1000]

    top_df <- data.frame(
      trait = rownames(w_matrix)[top_indices],
      LF4_loading = w_matrix4[top_indices]
    )

    cell_types <- sapply(strsplit(top_df$trait, "_"), function(x) x[2])

    cell_type_counts <- table(cell_types)
    cell_type_counts <- as.data.frame(cell_type_counts)

    missing_cell_types <- setdiff(required_cell_types, cell_type_counts$cell_types)
    if (length(missing_cell_types) > 0) {
      missing_rows <- data.frame(cell_types = missing_cell_types, Freq = 0)
      cell_type_counts <- rbind(cell_type_counts, missing_rows)
    }
    cell_type_counts$cell_types <- factor(cell_type_counts$cell_types, levels = required_cell_types)
    cell_type_counts <- cell_type_counts[order(cell_type_counts$cell_types), ]
    colnames(cell_type_counts)[2] <- paste0("PEER", i)

    peer_list[[i]] <- cell_type_counts
  }
  all_peers <- Reduce(function(x, y) merge(x, y, by = "cell_types", all = TRUE), peer_list)

  rownames(all_peers) <- all_peers$cell_types
  all_peers <- all_peers[, -1]
  all_peers <- all_peers[required_cell_types, ]
  heatmap_all <- pheatmap(all_peers,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          display_numbers = TRUE,
                          main = study, legend = TRUE,
                          fontsize_row = 16,
                          fontsize_col = 16)

 ggsave(heatmap_all, filename = paste0("results/Ling/", study, "/heatmaptop1000.png"),width = 5, height = 5)
}

### look at overlap of 1000 genes
studies <- c("CMC", "Ling", "SZBDMulti-Seq", "MultiomeBrain", "Batiuk")

ling_loadings <- read.table("data/data_raw/41586_2024_7109_MOESM6_ESM.Ling.GeneLoadings.txt", header=T, sep='\t') |>
  dplyr::select(id, PEER_4)
ling_loadings$PEER_4 <- abs(ling_loadings$PEER_4)
ling_loadings <- ling_loadings[order(ling_loadings$PEER_4, decreasing = TRUE), ]
ling_loadings <- ling_loadings[1:1000, ]
ling_genes <- sub("_.*", "", ling_loadings$id)

overlaps <- list()
for (study in studies){
  cpm_pb <- readRDS(paste0("data/data_processed/", study, "/Pseudobulk/", study, "_combined_cpm_matrix.rds"))
  gene_names <- colnames(cpm_pb)

  if (study != "Ling"){
    gene_names <- sub("^[^.]+\\.", "", gene_names)}
  gene_names <- gene_names[-1]

  w_matrix <- fread(paste0("results/Ling/", study, "/nocovs/W.csv")) |> t() |> data.frame()
  rownames(w_matrix) <- gene_names

  w_matrix4 <- w_matrix[,4]
  abs_loadings <- abs(w_matrix4)
  top_indices <- order(abs_loadings, decreasing = TRUE)[1:1000]
  top_traits <- rownames(w_matrix)[top_indices]
  top_df <- data.frame(
    trait = rownames(w_matrix)[top_indices],
    LF4_loading = w_matrix4[top_indices]
  )

  genes <- sub("_.*", "", top_df$trait)
  common_genes <-intersect(genes, ling_genes)
  overlaps[[study]] <- common_genes
}
}
### helpers
do_cpm_log <- function(mtx, log = FALSE) {
  colsums <- Matrix::colSums(mtx)
  cpm_result <- Matrix::t(Matrix::t(mtx) / colsums * 1e6)

  if (log) {
    cpm_result <- base::log1p(cpm_result)
  }

  return(cpm_result)
}

main()
