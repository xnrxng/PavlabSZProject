# Author: Rui Xiang Yu
# Date: 2024 September 19th
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/countcellsandgenes.R

library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)

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
  
  summary_df <- data.frame(
    Individual_ID = character(),
    Cohort = character(),
    Cell_Type = character(),
    Cell_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  summary_df <- process_sample(CMC_list, "CMC", summary_df)
  summary_df <- process_sample(SZBD_list, "SZBDMulti-Seq", summary_df)
  summary_df <- process_sample(MB_list, "MultiomeBrain", summary_df)
  
  celltype_perpatient <- summary_df |>
    filter(Cell_Type != 'featurekey')

  celltypeplot <- ggplot(celltype_perpatient, aes(x = Cell_Count, y = Cell_Type, color = Cohort)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 3, alpha = 0.8) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type"
    ) +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, 'lines')        
    ) +
    scale_color_brewer(palette = "Set2")
  
  ggsave(file.path("results/6-celltypedistribution.png"), celltypeplot)
  
  depth_df <- data.frame(
    Individual_ID = character(),
    Cohort = character(),
    Cell_Type = character(),
    Seq_Depth = integer(),
    stringsAsFactors = FALSE
  )

  depth_df <- sum_columns(CMC_list, "CMC", depth_df)
  depth_df <- sum_columns(SZBD_list, "SZBDMulti-Seq", depth_df)
  depth_df <- sum_columns(MB_list, "MultiomeBrain", depth_df)
  
  depthplot <- ggplot(depth_df, aes(x = Seq_Depth, y = Cell_Type, color = Cohort)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 3, alpha = 0.8) +
    labs(x = "Sequencing Depth (number  of reads)",
         y = "Cell type") +
    theme_minimal() + 
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 10),
      legend.key.size = unit(1.2, 'lines')        
    ) +
    scale_color_brewer(palette = "Set2")
    
  ggsave(file.path("results/7-depthdistribution.png"), depthplot)
  
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

process_sample <- function(sample_list, cohort_name, summary_df) {
  for (sample in sample_list) {
    sample_file <- paste0("data/data_raw/", cohort_name, "/", sample, "-annotated_matrix.txt")
    sample_data <- data.table::fread(sample_file)
    cell_types <- colnames(sample_data)
    cell_type_counts <- as.data.frame(table(cell_types))
    cell_type_counts <- cell_type_counts |>
      rename(Cell_Type = cell_types, Cell_Count = Freq) |>
      mutate(Individual_ID = sample, Cohort = cohort_name)
    summary_df <- bind_rows(summary_df, cell_type_counts)
  }
  return(summary_df)
}

sum_columns <- function(sample_list, cohort_name, depth_df) {
  for (sample in sample_list) {
    sample_file <- paste0("data/data_raw/", cohort_name, "/", sample, "-annotated_matrix.txt")
    sample_data <- data.table::fread(sample_file)
    numeric_data <- sample_data[, sapply(sample_data, is.numeric), with = FALSE]
    seq_depth <- colSums(numeric_data)
    
    seq_depth_df <- data.frame(
      Cell_Type = colnames(numeric_data),
      Seq_Depth = as.integer(seq_depth),
      Cohort = cohort_name,
      Individual_ID = sample,
      stringsAsFactors = FALSE
    )

    depth_df <- bind_rows(depth_df, seq_depth_df)
  }
  return(depth_df)
}

main()