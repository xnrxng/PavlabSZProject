# Author: Rui Xiang Yu
# Date: 2024 September 19th
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/countcellsandgenes.R

library(tidyverse)
library(data.table)

main <- function() {
  clean_metadata <- read_csv("data/data_processed/clean_metadata.csv")
  
  CMC_list <- clean_metadata |>
    filter(Cohort == "CMC") |>
    pull(Individual_ID)
  
  SZBD_list <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq") |>
    pull(Individual_ID)
  
  MB_list <- clean_metadata |>
    filter(Cohort == "MultiomeBrain") |>
    pull(Individual_ID)
  
  results <- data.frame(
    Cohort = character(),
    Genes = integer(),
    Total_Cells = integer(),
    stringsAsFactors = FALSE
  )

  results <- process_cohort(CMC_list, "CMC", results)
  results <- process_cohort(SZBD_list, "SZBDMulti-Seq", results)
  results <- process_cohort(MB_list, "MultiomeBrain", results)

  write.csv(results, file.path("results/5-cells_genes.csv"), row.names = FALSE)
  
  }
  
process_cohort <- function(sample_list, cohort_name, results_df) {
  total_cells <- 0
  genes <- 0
  for (i in seq_along(sample_list)) {
    sample <- sample_list[i]
    sample_file <- paste0("data/data_raw/", cohort_name, "/", sample, "-annotated_matrix.txt")
    sample_data <- data.table::fread(sample_file)
    total_cells <- total_cells + ncol(sample_data)

    if (i == 1) {
      genes <- nrow(sample_data)
    }
  }

  results_df <- rbind(results_df, data.frame(
    Cohort = cohort_name,
    Genes = genes,
    Total_Cells = total_cells,
    stringsAsFactors = FALSE
  ))
  
  return(results_df)
}

main()