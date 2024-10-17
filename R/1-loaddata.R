# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script downloads the data.
# Usage: R/loaddata.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)

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
  saveRDS(CMC_meta, "data/data_processed/CMC/CMC-patients.rds")
  
  SZBD_meta <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(SZBD_meta, "data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patients.rds")
  
  MB_meta <- clean_metadata |>
    filter(Cohort == "MultiomeBrain") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(MB_meta, "data/data_processed/MultiomeBrain/MultiomeBrain-patients.rds")
  
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
  
  CMC_control <- clean_metadata |>
    filter(Cohort == "CMC" & Disorder == "Control") |>
    pull(Individual_ID)
  
  CMC_sz <- clean_metadata |>
    filter(Cohort == "CMC" & Disorder == "Schizophrenia") |>
    pull(Individual_ID)
  
  SZBD_control <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq" & Disorder == "Control") |>
    pull(Individual_ID)
  
  SZBD_sz <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq" & Disorder == "Schizophrenia") |>
    pull(Individual_ID)
  
  MB_control <- clean_metadata |>
    filter(Cohort == "MultiomeBrain" & Disorder == "Control") |>
    pull(Individual_ID)
  
  MB_sz <- clean_metadata |>
    filter(Cohort == "MultiomeBrain" & Disorder == "Schizophrenia") |>
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
}

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
  
  output_path <- paste0("data/data_processed/", cohort, "/", capitalize_first_three(cell_type), "_", cohort, "_SZ.rds")
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

cleanCtmat <- function(ctmat, geneThr = 0.05, sampleThr = 0.05) {
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
    full_path <- paste0("data/data_processed/", cohort, "/", cell)
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
    filtered_path <- paste0("data/data_processed/", cohort, "/Filtered/Filt_", cell)
    saveRDS(filteredfile, filtered_path)
    
    cpmlogged <- do_cpm_log(filtered, TRUE)
    cpmfile <- list(meta = filtered_meta, expr = cpmlogged)
    cpm_path <- paste0("data/data_processed/", cohort, "/CPMLog/CPM_", cell)
    saveRDS(cpmfile, cpm_path)
  }
  return(NULL)
}

main()
