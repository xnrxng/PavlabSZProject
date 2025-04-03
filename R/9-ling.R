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
  for (peer in peers){
    w_matrix <- fread(paste0("results/Ling/", study, "/", peer, "/W.csv")) |> t()
    pca_matrix <- readRDS(paste0("results/Ling/", study, "/", study, "_PCA.rds"))
    pca_scores <- pca_matrix$rotation
    pca_scores <- pca_scores[-1, ]

    peer_list <- list()

  if (study == "MultiomeBrain"){




  }




  }


}

mat1_vec <- as.vector(w_matrix)
mat2_vec <- as.vector(pca_scores)

df <- data.frame(Mat1 = mat1_vec, Mat2 = mat2_vec)

# Plot
test <- ggplot(df, aes(x = Mat1, y = Mat2)) +
  geom_point(alpha = 0.3) +
  labs(title = "Scatter Plot of Matrix Values",
       x = "PEER Values",
       y = "PCA Values") +
  theme_classic()+
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed")


print(test)





}
### helpers``
do_cpm_log <- function(mtx, log = FALSE) {
  colsums <- Matrix::colSums(mtx)
  cpm_result <- Matrix::t(Matrix::t(mtx) / colsums * 1e6)

  if (log) {
    cpm_result <- base::log1p(cpm_result)
  }

  return(cpm_result)
}

main()
