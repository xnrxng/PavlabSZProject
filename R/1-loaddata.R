# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script downloads the data and wrangles it.
# Usage: R/1-loaddata.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)
library(readxl)
library(biomaRt)

main <- function() {
  raw_metadata <- data.table::fread("https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt")
  saveRDS(raw_metadata, file.path("data/data_raw/raw_metadata.rds"))

  clean_metadata <- raw_metadata |>
    filter(Cohort != "ROSMAP") |>
    mutate(Disorder = ifelse(tolower(Disorder) == "control", "Control", Disorder))

  saveRDS(clean_metadata, file.path("data/data_processed/clean_metadata.rds"))

  ### metadata RDS files for each cohort
  CMC_meta <- clean_metadata |>
    filter(Cohort == "CMC") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(CMC_meta, "data/data_processed/CMC/CMC-patient.rds")

  SZBD_meta <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(SZBD_meta, "data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patient.rds")

  MB_meta <- clean_metadata |>
    filter(Cohort == "MultiomeBrain") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(MB_meta, "data/data_processed/MultiomeBrain/MultiomeBrain-patient.rds")

  ### creation of lists of IDs
  CMC_list <- clean_metadata |>
    filter(Cohort == "CMC") |>
    pull(Individual_ID)

  SZBD_list <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq" & Disorder != "Bipolar Disorder") |>
    pull(Individual_ID)

  MB_list <- clean_metadata |>
    filter(Cohort == "MultiomeBrain" & Disorder != "Bipolar Disorder") |>
    pull(Individual_ID)

  ### generate RDS cell type files
  common_genes_CMC <- generate_genes("CMC", CMC_list)
  common_genes_SZBD <- generate_genes("SZBDMulti-Seq", SZBD_list)
  common_genes_MB <- generate_genes("MultiomeBrain", MB_list)

  generate_bycelltype("Micro", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("Micro", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("Micro", common_genes_CMC, "CMC", CMC_list)

  generate_bycelltype("Astro", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("Astro", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("Astro", common_genes_CMC, "CMC", CMC_list)

  generate_bycelltype("Oligo", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("Oligo", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("Oligo", common_genes_CMC, "CMC", CMC_list)

  generate_bycelltype("OPC", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("OPC", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("OPC", common_genes_CMC, "CMC", CMC_list)

  generate_bycelltype("excitatory", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("excitatory", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("excitatory", common_genes_CMC, "CMC", CMC_list)

  generate_bycelltype("inhibitory", common_genes_MB, "MultiomeBrain", MB_list)
  generate_bycelltype("inhibitory", common_genes_SZBD, "SZBDMulti-Seq", SZBD_list)
  generate_bycelltype("inhibitory", common_genes_CMC, "CMC", CMC_list)


  ### filter each matrix as well as cpm log them
  savefiltered("CMC")
  savefiltered("SZBDMulti-Seq")
  savefiltered("MultiomeBrain")

  ### wrangle batiuk
  batiuk_annotations <- readRDS("data/data_raw/Batiuk/annotations_final.RDS")
  batiuk_matrices <- readRDS("data/data_raw/Batiuk/snRNA-seq_raw_countmatrices.RDS")
  batiuk_meta <- read_excel("data/data_raw/Batiuk/sciadv.abn8367_tables_s1_to_s6.xlsx", sheet = 1)

  colnames(batiuk_meta) <- batiuk_meta[1, ]
  batiuk_meta <- batiuk_meta[-1, ]
  batiuk_meta_clean <- batiuk_meta |>
    dplyr::select(Identifier, Diagnosis, Age, Gender) |>
    rename(patientID = Identifier, disorder = Diagnosis, age = Age, sex = Gender) |>
    mutate(brainRegion = "PFC", study = "Batiuk",
           disorder = ifelse(disorder == "Scz", "yes", "no"),
           sex = ifelse(sex == "F", "female", "male")) |>
    filter(patientID != "MB8  (technical replicate of MB8-2)")

  saveRDS(batiuk_meta_clean, "data/data_processed/Batiuk/Batiuk-patient.rds")


  exc_neurons <- c("L2_3_CUX2_FREM3", "L2_CUX2_LAMP5", "L5_6_THEMIS", "L5_6_FEZF2_TLE4", "L3_CUX2_PRSS12",
                   "L4_RORB_SCHLAP1", "L4_5_FEZF2_LRRK1", "L5_FEZF2_ADRA1A")

  inh_neurons <- c("ID2_LAMP5", "VIP", "ID2_PAX6", "ID2_NCKAP5", "PVALB", "SST")

  batiuk_df <- data.frame(
    cell_ID = names(batiuk_annotations$med),
    cell_type = unname(batiuk_annotations$med)
  )

  batiuk_annotations_clean <- batiuk_df |>
  mutate(cell_type = ifelse(cell_type %in% exc_neurons, "Exc", cell_type),
         cell_type = ifelse(cell_type %in% inh_neurons, "Inh", cell_type))

  batiuk_matrices$MB8 <- NULL

  glia_rds_file <- batiuk_cell_types(celltype = "Glia", sample_list = batiuk_matrices, metadata = batiuk_meta_clean, cell_IDs = batiuk_annotations_clean)
  exc_rds_file <- batiuk_cell_types(celltype = "Exc", sample_list = batiuk_matrices, metadata = batiuk_meta_clean, cell_IDs = batiuk_annotations_clean)
  inh_rds_file <- batiuk_cell_types(celltype = "Inh", sample_list = batiuk_matrices, metadata = batiuk_meta_clean, cell_IDs = batiuk_annotations_clean)

  saveRDS(glia_rds_file, "data/data_processed/Batiuk/0.Raw/Gli_Batiuk_SZ.rds")
  saveRDS(exc_rds_file, "data/data_processed/Batiuk/0.Raw/Exc_Batiuk_SZ.rds")
  saveRDS(inh_rds_file, "data/data_processed/Batiuk/0.Raw/Inh_Batiuk_SZ.rds")

  ### create filtered V1
  ast_filt <- create_filteredV1("Ast")
  saveRDS(ast_filt$cmc, "data/data_processed/CMC/FilteredV1/Ast_CMC_SZ.rds")
  saveRDS(ast_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Ast_SZBDMulti-Seq_SZ.rds")
  saveRDS(ast_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Ast_MultiomeBrain_SZ.rds")

  exc_filt <- create_filteredV1("Exc")
  saveRDS(exc_filt$cmc, "data/data_processed/CMC/FilteredV1/Exc_CMC_SZ.rds")
  saveRDS(exc_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Exc_SZBDMulti-Seq_SZ.rds")
  saveRDS(exc_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Exc_MultiomeBrain_SZ.rds")
  saveRDS(exc_filt$batiuk, "data/data_processed/Batiuk/FilteredV1/Exc_Batiuk_SZ.rds")

  inh_filt <- create_filteredV1("Inh")
  saveRDS(inh_filt$cmc, "data/data_processed/CMC/FilteredV1/Inh_CMC_SZ.rds")
  saveRDS(inh_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Inh_SZBDMulti-Seq_SZ.rds")
  saveRDS(inh_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Inh_MultiomeBrain_SZ.rds")
  saveRDS(inh_filt$batiuk, "data/data_processed/Batiuk/FilteredV1/Inh_Batiuk_SZ.rds")

  mic_filt <- create_filteredV1("Mic")
  saveRDS(mic_filt$cmc, "data/data_processed/CMC/FilteredV1/Mic_CMC_SZ.rds")
  saveRDS(mic_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Mic_SZBDMulti-Seq_SZ.rds")
  saveRDS(mic_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Mic_MultiomeBrain_SZ.rds")

  opc_filt <- create_filteredV1("Opc")
  saveRDS(opc_filt$cmc, "data/data_processed/CMC/FilteredV1/Opc_CMC_SZ.rds")
  saveRDS(opc_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Opc_SZBDMulti-Seq_SZ.rds")
  saveRDS(opc_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Opc_MultiomeBrain_SZ.rds")

  oli_filt <- create_filteredV1("Oli")
  saveRDS(oli_filt$cmc, "data/data_processed/CMC/FilteredV1/Oli_CMC_SZ.rds")
  saveRDS(oli_filt$szbd, "data/data_processed/SZBDMulti-Seq/FilteredV1/Oli_SZBDMulti-Seq_SZ.rds")
  saveRDS(oli_filt$mb, "data/data_processed/MultiomeBrain/FilteredV1/Oli_MultiomeBrain_SZ.rds")

  gli_raw <- readRDS("data/data_processed/Batiuk/0.Raw/Gli_Batiuk_SZ.rds")
  gli_raw$expr <- cleanCtmat(gli_raw$expr)
  gli_raw$meta <- gli_raw$meta[rownames(gli_raw$meta) %in% colnames(gli_raw$expr), ]
  saveRDS(gli_raw, "data/data_processed/Batiuk/FilteredV1/Gli_Batiuk_SZ.rds")

  ### save batiuk subtypes
  for (excitatory in exc_neurons){
    exc_matrix <- readRDS("data/data_processed/Batiuk/FilteredV1/Exc_Batiuk_SZ.rds")

    target_cells <- batiuk_df$cell_ID[batiuk_df$cell_type == excitatory]
    filtered_matrix <- exc_matrix$expr[, colnames(exc_matrix$expr) %in% target_cells]
    filtered_meta <- exc_matrix$meta[rownames(exc_matrix$meta) %in% colnames(filtered_matrix), ]
    final_list <- list(meta = filtered_meta, expr = filtered_matrix)

    final_path <- paste0("data/data_processed/Batiuk/Excitatory/FilteredV1/", excitatory, "_Batiuk_SZ.rds")
    saveRDS(final_list, final_path)
  }

  for (inhibitory in inh_neurons){
    inh_matrix <- readRDS("data/data_processed/Batiuk/FilteredV1/Inh_Batiuk_SZ.rds")

    target_cells <- batiuk_df$cell_ID[batiuk_df$cell_type == inhibitory]
    filtered_matrix <- inh_matrix$expr[, colnames(inh_matrix$expr) %in% target_cells]
    filtered_meta <- inh_matrix$meta[rownames(inh_matrix$meta) %in% colnames(filtered_matrix), ]
    final_list <- list(meta = filtered_meta, expr = filtered_matrix)

    final_path <- paste0("data/data_processed/Batiuk/Inhibitory/FilteredV1/", inhibitory, "_Batiuk_SZ.rds")
    saveRDS(final_list, final_path)
  }
}

### helper functions
generate_bycelltype <- function(cell_type, common_genes, cohort, sample_list) {
  final_list <- list()
  metadata <- list()

  excitatory_cell_types <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 CT",
                               "L6 IT Car3", "L5 ET", "L5/6 NP", "L6b")

  inhibitory_cell_types <- c("Sst", "Sst Chodl", "Pvalb", "Chandelier", "Pax6",
                               "Lamp5 Lhx6", "Lamp5", "Sncg", "Vip")

  for (sample in sample_list) {
    sample_file <- paste0("data/data_raw/", cohort, "/", sample, "-annotated_matrix.txt")
    sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
    rownames(sample_data) <- sample_data[, 1]
    sample_data <- sample_data[, -1]
    sample_data <- sample_data[rownames(sample_data) %in% common_genes, , drop = FALSE]
    sample_data[is.na(sample_data)] <- 0
    colnames(sample_data) <- gsub("\\..*", "", colnames(sample_data))

    individual_metadata <- clean_metadata[clean_metadata$Individual_ID == sample, ]

    cell_counter <- 1

    for (col_idx in seq_along(colnames(sample_data))){
      column <- colnames(sample_data)[col_idx]

      if (cell_type == "excitatory" && column %in% excitatory_cell_types) {
        new_column_name <- paste0(cohort, "_" , sample, "_", cell_type, cell_counter)
      } else if (cell_type == "inhibitory" && column %in% inhibitory_cell_types) {
        new_column_name <- paste0(cohort, "_" , sample, "_", cell_type, cell_counter)
      } else if (column == cell_type) {
        new_column_name <- paste0(cohort, "_" , sample, "_", cell_type, cell_counter)
      } else {
        next
      }
      cell_counter <- cell_counter + 1
      final_list[[new_column_name]] <- sample_data[, col_idx]

      metadata[[length(metadata) + 1]] <- data.frame(
        cell_ID = new_column_name,
        cohort = individual_metadata$Cohort,
        patientID = individual_metadata$Individual_ID,
        sex = individual_metadata$Biological_Sex,
        age = individual_metadata$Age_death,
        disorder = individual_metadata$Disorder,
        stringsAsFactors = FALSE
      )
    }
  }
  metadata <- do.call(rbind, metadata)
  metadata <- metadata |>
    mutate(disorder = ifelse(disorder == "Schizophrenia", "yes", "no"))
  rownames(metadata) <- metadata[, 1]
  metadata <- metadata[, -1]

  final_dataframe <- as.data.frame(final_list, check.names = FALSE)
  rownames(final_dataframe) <- common_genes
  final_matrix <- as.matrix(final_dataframe)
  final_dgCmatrix <- as(final_matrix, "dgCMatrix")
  finalRDSfile <- list(meta = metadata, expr = final_dgCmatrix)

  capitalize_first_three <- function(cell_type) {
    short_name <- substr(cell_type, 1, 3)
    short_name <- paste0(toupper(substr(short_name, 1, 1)), tolower(substr(short_name, 2, 3)))
    return(short_name)
  }

  output_path <- paste0("data/data_processed/", cohort, "/0.Raw/", capitalize_first_three(cell_type), "_", cohort, "_SZ.rds")
  saveRDS(finalRDSfile, output_path)
  return(finalRDSfile)
}

generate_genes <- function(cohort, sample_list) {
  gene_list <- future_lapply(sample_list, function(sample) {
    sample_file <- paste0("data/data_raw/", cohort, "/", sample, "-annotated_matrix.txt")
    sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
    rownames(sample_data) <- sample_data[, 1]
    sample_data <- sample_data[, -1]
    gene_names <- rownames(sample_data)
    return(gene_names)
    })
  all_genes <- unlist(gene_list)
  common_genes <- unique(all_genes)
  return(common_genes)
}

cleanCtmat <- function(ctmat, geneThr = 0.02, sampleThr = 0.05) {
  nCellsTot <- ncol(ctmat)
  nCellsThr <- geneThr * nCellsTot

  genesNN <- data.frame(gene = rownames(ctmat), nCells = rowSums(ctmat > 0), check.names = FALSE)  |>
    filter(nCells > nCellsThr)
  ctmat <- ctmat[genesNN$gene, ]

  cellsNN <- data.frame(cell = colnames(ctmat), nGenes = colSums(ctmat > 0), check.names = FALSE) |>
    mutate(perc = percent_rank(nGenes)) |>
    filter(perc > sampleThr)
  ctmat <- ctmat[, cellsNN$cell]

  return(ctmat)
}

do_cpm_log <- function(mtx, log) {
  colsums <- colSums(mtx)
  cpm_result <- t(t(mtx) / colsums * 1e6)

  if (log) {
    cpm_result <- log1p(cpm_result)
  }

  return(cpm_result)
}

savefiltered <- function(cohort) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")

  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)

  for (cell in cell_list) {
    full_path <- paste0("data/data_processed/", cohort, "/0.Raw/", cell)
    unfiltered <- readRDS(full_path)
    expr_matrix <- unfiltered$expr
    metadata <- unfiltered$meta

    less_than_20_cells <- metadata |>
      group_by(patientID) |>
      tally() |>
      filter(n <20) |>
      pull(patientID)

    columns_to_remove <- sapply(colnames(expr_matrix), function(col) {
      any(sapply(less_than_20_cells, function(str) grepl(str, col)))
    })

    matrix_patientsplus20 <- expr_matrix[, !columns_to_remove]

    matrix_columns <- colnames(matrix_patientsplus20)
    meta_rows <- rownames(metadata)
    common_rows <- meta_rows[meta_rows %in% matrix_columns]
    filtered_meta <- metadata[common_rows, , drop = FALSE]

    filtered <- cleanCtmat(matrix_patientsplus20)
    filteredfile <- list(meta = filtered_meta, expr = filtered)
    filtered_path <- paste0("data/data_processed/", cohort, "/Scrapped/Filtered/", cell)
    saveRDS(filteredfile, filtered_path)

    cpmlogged <- do_cpm_log(filtered, TRUE)
    cpmfile <- list(meta = filtered_meta, expr = cpmlogged)
    cpm_path <- paste0("data/data_processed/", cohort, "/Scrapped/FilteredandCPMLog/", cell)
    saveRDS(cpmfile, cpm_path)
  }
  return(NULL)
}

batiuk_cell_types <- function(celltype, sample_list, metadata, cell_IDs) {
  final_list <- list()
  meta_list <- list()

  all_genes <- rownames(sample_list$MB7)

  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_annotations <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                            filters = "external_gene_name", values = rownames(sample_list$MB7), mart = ensembl)
  unique_genes <- unique(gene_annotations$external_gene_name)

  keep_genes <- rownames(sample_list$MB7[unique_genes, ])

  patient_IDs <- metadata |> dplyr::select(patientID) |> pull()
  valid_cell_IDs <- cell_IDs[cell_IDs$cell_type == celltype, "cell_ID"]

  for (patient in patient_IDs) {
    sample <- sample_list[[patient]]
    columns_to_keep <- intersect(colnames(sample), valid_cell_IDs)
    sample_filtered <- sample[keep_genes, columns_to_keep]

    individual_metadata <- metadata[metadata$patientID == patient, ]

    for (column in colnames(sample_filtered)){
      final_list[[column]] <- sample_filtered[, column]

      meta_list[[length(meta_list) + 1]] <- data.frame(
        cell_ID = column,
        cohort = "Batiuk",
        patientID = individual_metadata$patientID,
        sex = individual_metadata$sex,
        age = individual_metadata$age,
        disorder = individual_metadata$disorder,
        stringsAsFactors = FALSE)
    }
  }

  cell_metadata <- do.call(rbind, meta_list)
  rownames(cell_metadata) <- cell_metadata[, 1]
  cell_metadata <- cell_metadata[, -1]

  final_dataframe <- as.data.frame(final_list, check.names = FALSE)
  rownames(final_dataframe) <- keep_genes
  final_matrix <- as.matrix(final_dataframe)
  final_dgCmatrix <- as(final_matrix, "dgCMatrix")
  finalRDSfile <- list(meta = cell_metadata, expr = final_dgCmatrix)
  return(finalRDSfile)
}

create_filteredV1 <- function(cell_type) {
  batiuk_path <- paste0("data/data_processed/Batiuk/0.Raw/", cell_type, "_Batiuk_SZ.rds")
  cmc_path <- paste0("data/data_processed/CMC/0.Raw/", cell_type, "_CMC_SZ.rds")
  szbd_path <- paste0("data/data_processed/SZBDMulti-Seq/0.Raw/", cell_type, "_SZBDMulti-Seq_SZ.rds")
  mb_path <- paste0("data/data_processed/MultiomeBrain/0.Raw/", cell_type, "_MultiomeBrain_SZ.rds")

  all_data <- list()

  if (file.exists(batiuk_path)) {
    batiuk <- readRDS(batiuk_path)
    batiuk$expr <- cleanCtmat(batiuk$expr)
    all_data$batiuk <- batiuk
  } else {
    message("File not found: ", batiuk_path)
  }

  if (file.exists(cmc_path)) {
    cmc <- readRDS(cmc_path)
    cmc$expr <- cleanCtmat(cmc$expr)
    all_data$cmc <- cmc
  } else {
    message("File not found: ", cmc_path)
  }

  if (file.exists(szbd_path)) {
    szbd <- readRDS(szbd_path)
    szbd$expr <- cleanCtmat(szbd$expr)
    all_data$szbd <- szbd
  } else {
    message("File not found: ", szbd_path)
  }

  if (file.exists(mb_path)) {
    mb <- readRDS(mb_path)
    mb$expr <- cleanCtmat(mb$expr)
    all_data$mb <- mb
  } else {
    message("File not found: ", mb_path)
  }

  all_genes <- c()
  for (file in all_data){
    genes <- rownames(file$expr)
    all_genes <- c(all_genes, genes)
  }

  string_counts <- table(all_genes)
  filtered_genes <- all_genes[all_genes %in% names(string_counts[string_counts >= 3])]
  filtered_genes <- unique(filtered_genes)

  for (name in names(all_data)){
    file <- all_data[[name]]
    missing_genes <- setdiff(filtered_genes, rownames(file$expr))
    zero_matrix <- Matrix(0, nrow = length(missing_genes), ncol = ncol(file$expr), dimnames = list(missing_genes, colnames(file$expr)))
    filt <- rbind(file$expr, zero_matrix)
    file$expr <- filt[rownames(filt) %in% filtered_genes,]
    file$meta <- file$meta[rownames(file$meta) %in% colnames(file$expr), ]

    all_data[[name]] <- file
  }

  return(all_data)
}

main()
