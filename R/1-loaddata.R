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
  saveRDS(CMC_meta, "data/data_processed/CMC-patients.rds")
  
  SZBD_meta <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(SZBD_meta, "data/data_processed/SZBDMulti-Seq-patients.rds")
  
  MB_meta <- clean_metadata |>
    filter(Cohort == "MultiomeBrain") |>
    rename(patientID = Individual_ID, sex = Biological_Sex, disorder = Disorder, study = Cohort, age = Age_death) |>
    mutate(brainRegion = "PFC",
           disorder = ifelse(disorder == "Schizophrenia", "yes", "no")) |>
    select(patientID, sex, age, disorder, brainRegion, study)
  saveRDS(MB_meta, "data/data_processed/MultiomeBrain-patients.rds")
  
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
  common_genes <- generate_genes()
  
  generate_bycelltype("excitatory", common_genes)
  generate_bycelltype("inhibitory", common_genes)
  generate_bycelltype("Astro", common_genes)
  generate_bycelltype("Oligo", common_genes)
  generate_bycelltype("OPC", common_genes)
  generate_bycelltype("Micro", common_genes)
  
}

generate_bycelltype <- function(cell_type, common_genes) {
  
  final_list <- list()
  
  metadata <- list() 
  
  generate_individual <- function(cohort, sample_list, final_list, metadata) {
    excitatory_cell_types <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 CT", 
                               "L6 IT Car3", "L5 ET", "LE/6 NP", "L6b IT")
    
    inhibitory_cell_types <- c("Sst", "Sst Chodl", "Pvalb", "Chandelier", "Pax6",
                               "Lamp5 Lhx6", "Lamp5", "Sncg", "Vip")
    
    for (sample in sample_list) {
      sample_file <- paste0("data/data_raw/", cohort, "/", sample, "-annotated_matrix.txt")
      sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
      rownames(sample_data) <- sample_data[, 1]
      sample_data <- sample_data[, -1]
      sample_data <- sample_data[rownames(sample_data) %in% common_genes, , drop = FALSE]
      sample_data[is.na(sample_data)] <- 0
      
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
    return(list(metadata = metadata, final_list = final_list))
  }
  
  result1 <- generate_individual("CMC", CMC_list, final_list, metadata)
  result2 <- generate_individual("SZBDMulti-Seq", SZBD_list, final_list, metadata)
  result3 <- generate_individual("MultiomeBrain", MB_list, final_list, metadata)
  
  metadata <- do.call(rbind, c(result1$metadata, result2$metadata, result3$metadata))
  final_list <- c(result1$final_list, result2$final_list, result3$final_list)
  
  metadata <- metadata |> mutate(disorder = ifelse(disorder == "Schizophrenia", "yes", "no"))
  rownames(metadata) <- metadata[, 1]
  metadata <- metadata[, -1]
  
  final_dataframe <- as.data.frame(final_list)
  rownames(final_dataframe) <- common_genes
  final_matrix <- as.matrix(final_dataframe)
  final_dgCmatrix <- as(final_matrix, "dgCMatrix") 
  
  finalRDSfile <- list(meta = metadata, expr = final_dgCmatrix)
  
  output_path <- paste0("data/data_processed/", cell_type, ".rds")
  saveRDS(finalRDSfile, output_path)
}

generate_genes <- function() {
  all_genes_list <- list()
  
  get_all_genes <- function(cohort, sample_list) {
    gene_list <- future_lapply(sample_list, function(sample) {
      sample_file <- paste0("data/data_raw/", cohort, "/", sample, "-annotated_matrix.txt")
      sample_data <- fread(sample_file, header = TRUE, sep = "\t", data.table = FALSE)
      rownames(sample_data) <- sample_data[, 1]
      sample_data <- sample_data[, -1]
      return(rownames(sample_data))
    })
    return(gene_list)
  }
  
  all_genes_list$CMC <- get_all_genes("CMC", CMC_list)
  all_genes_list$SZBDMultiSeq <- get_all_genes("SZBDMulti-Seq", SZBD_list)
  all_genes_list$MultiomeBrain <- get_all_genes("MultiomeBrain", MB_list)
  common_genes <- Reduce(dplyr::intersect, unlist(all_genes_list, recursive = FALSE))
  
  return(common_genes)
}

main()