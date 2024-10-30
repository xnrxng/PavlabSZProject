# Author: Rui Xiang Yu
# Date: 2024 October 29th
# This script downloads the ROSMAP Alzheimers data and wrangles it.
# Usage: R/5-alzheimerswrangle.R

library(tidyverse)
library(data.table)
library(Matrix)
library(future.apply)

main <- function() {
  
  ### create patient metadata
  raw_clinical_data <- fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/meta/ROSMAP_clinical.csv") |>
    select(individualID, msex, age_death) |>
    rename(patientID = individualID, sex = msex, age = age_death) |>
    mutate(sex = ifelse(sex == 0, "female", "male"))
  
  raw_metadata <- fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/meta/snRNAseqAD_TREM2_biospecimen_metadata.csv") |>
    select(individualID, specimenID, tissue, BrodmannArea) |>
    rename(patientID = individualID) |>
    mutate(brainRegion = paste0(tissue, ", ", BrodmannArea)) |>
    select(patientID, specimenID, brainRegion)

  ad_specimenIDs <- c("AD1", "AD2", "AD3", "AD5", "AD7", "AD8", "AD9", "AD10", "AD11", "AD12", "AD13")
  ctrl_specimenIDs <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C11", "C12")
  trem_specimenIDs <- c("P2", "P3", "P5", "P6", "P7", "P9", "P10", "P11", "P12", "P13")
  all_specimenIDs <- c(ad_specimenIDs, ctrl_specimenIDs, trem_specimenIDs)
  
  filtered_meta <- merge(raw_clinical_data, raw_metadata, by = "patientID", all = TRUE)
  filtered_meta <- filtered_meta[filtered_meta$specimenID %in% all_specimenIDs, ]
  filtered_meta <- filtered_meta |>
    mutate(study = "ROSMAP-TREM2-2020")
  filtered_meta$ADdiag2types <- ifelse(filtered_meta$specimenID %in% c(ad_specimenIDs, trem_specimenIDs), "yes", "no")
  filtered_meta$ADdiag3types <- ifelse(filtered_meta$specimenID %in% ctrl_specimenIDs, "Ctrl",
                            ifelse(filtered_meta$specimenID %in% ad_specimenIDs, "WT",
                                   ifelse(filtered_meta$specimenID %in% trem_specimenIDs, "R62H", "no")))
  
  ### wrangle cell ID  
  cell_IDs <- readxl::read_excel("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/meta/clusters_cellID.xlsx", sheet = 2) |>
    select(Barcodes, Sample, Label)
  cell_IDs$Barcodes <- gsub("^TWCC-", "", cell_IDs$Barcodes)
  cell_IDs$Sample <- gsub("^TWCC-", "", cell_IDs$Sample)
  cell_IDs$Label <- gsub("\\d+$", "",  cell_IDs$Label)
  
  ### generate final RDS file
  create_combined_matrix(all_specimenIDs)
  common_genes <- generate_genes(all_specimenIDs)
  
  mic_file <- generate_bycelltype("Micro", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(mic_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Mic_ROSMAP-TREM2-2020.rds")
  
  exc_file <- generate_bycelltype("Ex", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(exc_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Exc_ROSMAP-TREM2-2020.rds")
  
  oli_file <- generate_bycelltype("Oli", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(oli_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Oli_ROSMAP-TREM2-2020.rds")
  
  ast_file <- generate_bycelltype("Astro", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(ast_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Ast_ROSMAP-TREM2-2020.rds")
  
  opc_file <- generate_bycelltype("OPC", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(opc_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Opc_ROSMAP-TREM2-2020.rds")
  
  inh_file <- generate_bycelltype("In", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(inh_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/Inh_ROSMAP-TREM2-2020.rds")
  
  end_file <- generate_bycelltype("Endo", common_genes, all_specimenIDs, filtered_meta, cell_IDs)
  saveRDS(end_file, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/End_ROSMAP-TREM2-2020.rds")
  
  ### filter and cpmlog
  savefiltered()
  
  ### 
  filtered_meta <- select(filtered_meta, -specimenID)
  saveRDS(filtered_meta, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/ROSMAP-TREM2-2020-patient.rds")
}

### helper functions
create_combined_matrix <- function(specimenIDs) {
  for (patient in specimenIDs){
    barcode_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/expr/", patient, "_barcodes.tsv")
    barcodes <- fread(barcode_file, header = FALSE) |>
      mutate(V1 = paste0(patient, "_", sub("-1$", "", V1))) |>
      pull()
    
    gene_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/expr/", patient, "_features.tsv")
    genes <- fread(gene_file, header = FALSE) |>
      select(V2) |>
      pull()
    
    matrix_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/0-source/expr/", patient, "_matrix.mtx")
    matrix_expr <- readMM(matrix_file)
    rownames(matrix_expr) <- genes
    colnames(matrix_expr) <- barcodes
    
    final_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/matrices/", patient, "combined_matrix.rds")
    saveRDS(matrix_expr, final_file)
  }}

generate_genes <- function(sample_list) {
  gene_list <- future_lapply(sample_list, function(sample) {
    sample_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/matrices/", sample, "combined_matrix.rds")
    sample_data <- readRDS(sample_file)
    gene_names <- rownames(sample_data)
    return(gene_names)
  })
  all_genes <- unlist(gene_list)
  common_genes <- unique(all_genes)
  return(common_genes)
}

generate_bycelltype <- function(cell_type, common_genes, sample_list, metadata, cell_IDs) {
  final_list <- list()
  meta_list <- list() 
  
  valid_cell_IDs <- cell_IDs[cell_IDs$Label == cell_type, "Barcodes"] |>
    na.omit() |>
    pull()
  
  for (sample in sample_list) {
    sample_file <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/matrices/", sample, "combined_matrix.rds")
    sample_data <- readRDS(sample_file)
    
    standardized_data <- data.frame(matrix(NA, nrow = length(common_genes), ncol = ncol(sample_data)))
    rownames(standardized_data) <- common_genes
    colnames(standardized_data) <- colnames(sample_data)
    present_genes <- rownames(sample_data)[rownames(sample_data) %in% common_genes]
    standardized_data[present_genes, ] <- sample_data[present_genes, , drop = FALSE]
    standardized_data <- standardized_data[, colnames(standardized_data) %in% valid_cell_IDs, drop = FALSE]
    
    individual_metadata <- metadata[metadata$specimenID == sample, ]
    
    for (column in colnames(standardized_data)){
          final_list[[column]] <- standardized_data[, column]
          
          meta_list[[length(meta_list) + 1]] <- data.frame(
            cell_ID = column,
            patientID = individual_metadata$patientID,
            sex = individual_metadata$sex,
            age = individual_metadata$age,
            brainRegion = individual_metadata$brainRegion,
            ADdiag2types = individual_metadata$ADdiag2types,
            ADdiag3types = individual_metadata$ADdiag3types,
            stringsAsFactors = FALSE)
        }
      }
  
  cell_metadata <- do.call(rbind, meta_list)
  rownames(cell_metadata) <- cell_metadata[, 1]
  cell_metadata <- cell_metadata[, -1]
  
  final_dataframe <- as.data.frame(final_list, check.names = FALSE)
  rownames(final_dataframe) <- common_genes
  final_matrix <- as.matrix(final_dataframe)
  final_dgCmatrix <- as(final_matrix, "dgCMatrix") 
  finalRDSfile <- list(meta = cell_metadata, expr = final_dgCmatrix)
  return(finalRDSfile)
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

savefiltered <- function() {
  cell_list <- c("Exc", "Inh", "Opc", "End", "Mic", "Ast", "Oli")
  
  for (cell in cell_list) {
    full_path <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/raw/", cell, "_ROSMAP-TREM2-2020.rds")
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
    filtered_path <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/filtered/Filt_", cell, "_ROSMAP-TREM2-2020.rds")
    saveRDS(filteredfile, filtered_path)
    
    cpmlogged <- do_cpm_log(filtered, TRUE)
    cpmfile <- list(meta = filtered_meta, expr = cpmlogged)
    cpm_path <- paste0("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/filteredandcpmlog/CPM_", cell, "_ROSMAP-TREM2-2020.rds")
    saveRDS(cpmfile, cpm_path)
  }
  return(NULL)
}

main()