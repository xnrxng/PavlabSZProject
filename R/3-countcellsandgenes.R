# Author: Rui Xiang Yu
# Date: 2024 September 19th
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/countcellsandgenes.R

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggforce)

main <- function() {
  clean_metadata <- readRDS("data/data_processed/clean_metadata.rds")
  
  CMC_list <- clean_metadata |>
    filter(Cohort == "CMC") |>
    pull(Individual_ID)
  
  SZBD_list <- clean_metadata |>
    filter(Cohort == "SZBDMulti-Seq" & Disorder != "Bipolar Disorder") |>
    pull(Individual_ID)
  
  MB_list <- clean_metadata |>
    filter(Cohort == "MultiomeBrain" & Disorder != "Bipolar Disorder") |>
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

  saveRDS(results, file.path("results/5-cells_genes.rds"))

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
  
  ### save filtered matrixes
  for (CMCsample in CMC_list) {
    sample_file <- paste0("data/data_raw/CMC/", CMCsample, "-annotated_matrix.txt")
    sample_data <- read.table(sample_file, check.names = FALSE, header = TRUE, row.names = 1, sep ="\t")
    filtered <- cleanCtmat(as.data.frame(sample_data))
    if (ncol(filtered) < 20) {
      next}
    output_path <- paste0("data/data_processed/CMC/", CMCsample, "-annotated_matrix.rds")
    write_rds(filtered, output_path)
    cpmlognormalized <- do_cpm_log(filtered, log=TRUE)
    CPM_path <- paste0("data/data_processed/CMC_CPM/", CMCsample, "-CPM_annotated_matrix.rds")
    write_rds(cpmlognormalized, CPM_path)
  }
  
  for (SZBDsample in SZBD_list) {
    sample_file <- paste0("data/data_raw/SZBDMulti-Seq/", SZBDsample, "-annotated_matrix.txt")
    sample_data <- read.table(sample_file, check.names = FALSE, header = TRUE, row.names = 1, sep ="\t")
    filtered <- cleanCtmat(as.data.frame(sample_data))
    if (ncol(filtered) < 20) {
      next}
    output_path <- paste0("data/data_processed/SZBDMulti-Seq/", SZBDsample, "-annotated_matrix.rds")
    write_rds(filtered, output_path)
    cpmlognormalized <- do_cpm_log(filtered, log=TRUE)
    CPM_path <- paste0("data/data_processed/SZBDMulti-Seq_CPM/", SZBDsample, "-CPM_annotated_matrix.rds")
    write_rds(cpmlognormalized, CPM_path)
  }
  
  for (MBsample in MB_list) {
    sample_file <- paste0("data/data_raw/MultiomeBrain/", MBsample, "-annotated_matrix.txt")
    sample_data <- read.table(sample_file, check.names = FALSE, header = TRUE, row.names = 1, sep ="\t")
    filtered <- cleanCtmat(as.data.frame(sample_data))
    if (ncol(filtered) < 20) {
      next}
    output_path <- paste0("data/data_processed/MultiomeBrain/", MBsample, "-annotated_matrix.rds")
    write_rds(filtered, output_path)
    cpmlognormalized <- do_cpm_log(filtered, log=TRUE)
    CPM_path <- paste0("data/data_processed/MultiomeBrain_CPM/", MBsample, "-CPM_annotated_matrix.rds")
    write_rds(cpmlognormalized, CPM_path)
  }
  
  ### results data frame but with filtered matrices
  results_filtered <- data.frame(
    Cohort = character(),
    Genes = integer(),
    Total_Cells = integer(),
    stringsAsFactors = FALSE
  )
  
  results_filtered <- process_cohort_filtered(CMC_list, "CMC", results_filtered)
  results_filtered <- process_cohort_filtered(SZBD_list, "SZBDMulti-Seq", results_filtered)
  results_filtered <- process_cohort_filtered(MB_list, "MultiomeBrain", results_filtered)
  
  write_rds(results_filtered, file.path("results/5.2-cells_genes.rds"))
  
  ### cell type plot
  summary_df <- data.frame(
    Individual_ID = character(),
    Cohort = character(),
    Cell_Type = character(),
    Cell_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  summary_df <- process_sample(CMC_control, CMC_sz, "CMC", summary_df)
  summary_df <- process_sample(SZBD_control, SZBD_sz, "SZBDMulti-Seq", summary_df)
  summary_df <- process_sample(MB_control, MB_sz, "MultiomeBrain", summary_df)
  
  celltype_perpatient <- summary_df |>
    filter(Cell_Type != 'featurekey')

  CMCcelltypeplot <- celltype_perpatient |>
    filter(Cohort == "CMC") |>
    ggplot(aes(x = Cell_Count, y = Cell_Type, color = Disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.5) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "CMC") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
    xlim(c(0, 5000)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      strip.background = element_blank(),
      panel.spacing = unit(1, "lines")) +
    scale_color_brewer(palette = "Set2")
  
  SZBDcelltypeplot <- celltype_perpatient |>
    filter(Cohort == "SZBDMulti-Seq") |>
    ggplot(aes(x = Cell_Count, y = Cell_Type, color = Disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.5) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "SZBDMulti-Seq") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
    xlim(c(0, 5000)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = "none",
      strip.background = element_blank(),
      panel.spacing = unit(1, "lines")) +
    scale_color_brewer(palette = "Set2")
  
  MBcelltypeplot <- celltype_perpatient |>
    filter(Cohort == "MultiomeBrain") |>
    ggplot(aes(x = Cell_Count, y = Cell_Type, color = Disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.5) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "MultiomeBrain") +
    scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000)) +
    xlim(c(0, 5000)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_blank(),
      panel.spacing = unit(1, "lines")) +
    scale_color_brewer(palette = "Set2")
  
  celltypecombined_plots <- plot_grid(
    plotlist = list(CMCcelltypeplot, SZBDcelltypeplot, MBcelltypeplot+ theme(legend.position = "none")), ncol = 1, align = "hv")
  
  disorderlegend <- get_legend(MBcelltypeplot+ theme(legend.box.margin = margin(0, 0, 0, 3),
                                             legend.key.size = unit(1, 'cm'),
                                             legend.title = element_text(size = 15)))
  
  x_label_plot <- ggdraw() + 
    draw_label("Number of Cells (Frequency)", x = 0.5, y = 0.5, size = 14)
  
  celltypefinal_plot <- plot_grid(combined_plots, x_label_plot, ncol = 1, rel_heights = c(1, 0.1))
  
  celltypefinal_plot_with_legend <- plot_grid(
    final_plot,
    legend,
    ncol = 2,  
    rel_widths = c(0.85, 0.15)
  )
  
  ggsave(file.path("results/6-celltypedistribution.png"), celltypefinal_plot_with_legend, width = 10, height = 15)
  
  ###sequencing depth plot
  depth_df <- data.frame(
    Individual_ID = character(),
    Cohort = character(),
    Cell_Type = character(),
    Seq_Depth = integer(),
    Disorder = character(),
    stringsAsFactors = FALSE
  )

  depth_df <- sum_columns(CMC_control, CMC_sz, "CMC", depth_df)
  depth_df <- sum_columns(SZBD_control, SZBD_sz, "SZBDMulti-Seq", depth_df)
  depth_df <- sum_columns(MB_control, MB_sz, "MultiomeBrain", depth_df)
  
  depth_df_cpm <- data.frame(
    Individual_ID = character(),
    Cohort = character(),
    Cell_Type = character(),
    Seq_Depth = integer(),
    Disorder = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_cpm <- sum_columns_cpm(CMC_control, CMC_sz, "CMC", depth_df_cpm)
  depth_df_cpm <- sum_columns_cpm(SZBD_control, SZBD_sz, "SZBDMulti-Seq", depth_df_cpm)
  depth_df_cpm <- sum_columns_cpm(MB_control, MB_sz, "MultiomeBrain", depth_df_cpm)
  
  depthplot <- ggplot(depth_df, aes(x = Seq_Depth, y = Cell_Type, color = Disorder)) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = Disorder), size=0.4, alpha=0.5) +
    xlim(0, 300000) +
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
      strip.placement = "outside") +
    scale_color_brewer(palette = "Set2") +
    facet_wrap(~ Cohort, ncol =1, scales = "free_x")
      
  ggsave(file.path("results/7-depthdistribution.png"), depthplot, width = 8, height = 15)

  ### cell type abundance per condition plot
  abundance_df <- data.frame(Cohort = character(),
                             Disorder = character(),
                             Cell_Type = character(),
                             Cell_Count = integer(),
                             stringsAsFactors = FALSE)

  abundance_df <- abundance_conditions(CMC_control, CMC_sz, "CMC", abundance_df)
  abundance_df <- abundance_conditions(SZBD_control, SZBD_sz, "SZBDMulti-Seq", abundance_df)
  abundance_df <- abundance_conditions(MB_control, MB_sz, "MultiomeBrain", abundance_df)

  abundance_plot <- ggplot(abundance_df, aes(x = Cell_Type, y = Cell_Count, fill = Disorder)) +
    geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Cohort, ncol=1) +
  labs(x = "Cell Type",
       y = "Cell Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")) +
    scale_fill_brewer(palette = "Pastel1")
  
  ggsave(file.path("results/8-abundancecondition.png"), abundance_plot, width = 8, height = 12)
    
  }


### helper functions  
process_cohort <- function(sample_list, cohort_name, results_df) {
  total_cells <- 0
  genes <- 0
  for (i in seq_along(sample_list)) {
    sample <- sample_list[i]
    sample_file <- paste0("data/data_raw/", cohort_name, "/", sample, "-annotated_matrix.txt")
    sample_data <- read.table(sample_file, check.names = FALSE, header = TRUE, row.names = 1, sep ="\t")
    total_cells <- total_cells + ncol(sample_data)

    if (i == 1) {
      genes <- nrow(sample_data)
    }
  }

  results_df <<- rbind(results_df, data.frame(
    Cohort = cohort_name,
    Genes = genes,
    Total_Cells = total_cells,
    stringsAsFactors = FALSE
  ))
  
  return(results_df)
}

process_cohort_filtered <- function(sample_list, cohort_name, results_df) {
  total_cells <- 0
  genes <- 0
  for (i in seq_along(sample_list)) {
    sample <- sample_list[i]
    sample_file <- paste0("data/data_processed/", cohort_name, "/", sample, "-annotated_matrix.rds")
    sample_data <- read_rds(sample_file)
    total_cells <- total_cells + ncol(sample_data)
    
    if (i == 1) {
      genes <- nrow(sample_data)
    }
  }
  
  results_df <<- rbind(results_df, data.frame(
    Cohort = cohort_name,
    Genes = genes,
    Total_Cells = total_cells,
    stringsAsFactors = FALSE
  ))
  
  return(results_df)
}

process_sample <- function(control_list, disease_list, cohort_name, summary_df) {
  
  process_individual <- function(sample_list, disorder_label) {
    for (sample in sample_list) {
      sample_file <- paste0("data/data_processed/", cohort_name, "/", sample, "-annotated_matrix.rds")
      sample_data <- read_rds(sample_file)
      cell_types <- colnames(sample_data)
      cell_type_counts <- as.data.frame(table(cell_types))
      cell_type_counts <- cell_type_counts |>
      rename(Cell_Type = cell_types, Cell_Count = Freq) |>
      mutate(Individual_ID = sample, Cohort = cohort_name, Disorder = disorder_label)
      summary_df <<- bind_rows(summary_df, cell_type_counts)}}

  process_individual(control_list, "Control")
  process_individual(disease_list, "Schizophrenia")
  
  return(summary_df)
}

sum_columns <- function(control_list, disease_list, cohort_name, depth_df) {
  
  sum_columns_individual <- function(sample_list, disorder_label) {
    for (sample in sample_list) {
      sample_file <- paste0("data/data_processed/", cohort_name, "/", sample, "-annotated_matrix.rds")
      sample_data <- read_rds(sample_file)
      numeric_data <- sample_data[, sapply(sample_data, is.numeric)]
      seq_depth <- colSums(numeric_data)
      cleaned_cell_type <- sub("\\.\\d+$", "", colnames(numeric_data))
      
      seq_depth_df <- data.frame(
      Cell_Type = cleaned_cell_type,
      Seq_Depth = as.integer(seq_depth),
      Cohort = cohort_name,
      Individual_ID = sample,
      Disorder = disorder_label,
      stringsAsFactors = FALSE)
      depth_df <<- bind_rows(depth_df, seq_depth_df)
      }}
  
  sum_columns_individual(control_list, "Control")
  sum_columns_individual(disease_list, "Schizophrenia")
  return(depth_df)
}

sum_columns_cpm <- function(control_list, disease_list, cohort_name, depth_df) {
  
  sum_columns_individual <- function(sample_list, disorder_label) {
    for (sample in sample_list) {
      sample_file <- paste0("data/data_processed/", cohort_name, "_CPM/", sample, "-CPM_annotated_matrix.rds")
      sample_data <- read_rds(sample_file)
      #numeric_data <- sample_data[, sapply(sample_data, is.numeric)]
      seq_depth <- colSums(sample_data)
      cleaned_cell_type <- sub("\\.\\d+$", "", colnames(sample_data))
      
      seq_depth_df <- data.frame(
        Cell_Type = cleaned_cell_type,
        Seq_Depth = as.integer(seq_depth),
        Cohort = cohort_name,
        Individual_ID = sample,
        Disorder = disorder_label,
        stringsAsFactors = FALSE)
      depth_df <<- bind_rows(depth_df, seq_depth_df)
    }}
  
  sum_columns_individual(control_list, "Control")
  sum_columns_individual(disease_list, "Schizophrenia")
  return(depth_df)
}

abundance_conditions <- function(control_list, disease_list, cohort_name, summary_df) {
  control_counts <- data.frame(Cell_Type = character(), Cell_Count = integer(), Disorder = character(), Cohort = character(), stringsAsFactors = FALSE)
  disease_counts <- data.frame(Cell_Type = character(), Cell_Count = integer(), Disorder = character(), Cohort = character(), stringsAsFactors = FALSE)
  
  for (control in control_list) {
    sample_file <- paste0("data/data_processed/", cohort_name, "/", control, "-annotated_matrix.rds")
    sample_data <- as.data.frame(read_rds(sample_file))
    cell_types <- colnames(sample_data)
    cell_type_counts <- as.data.frame(table(cell_types)) |>
      rename(Cell_Type = cell_types, Cell_Count = Freq)
    control_counts <<- bind_rows(control_counts, cell_type_counts)
  }
  
  control_counts <- control_counts %>%
    group_by(Cell_Type) %>%
    summarise(Cell_Count = sum(Cell_Count)) %>%
    mutate(Disorder = "Control", Cohort = cohort_name)
  
  for (sample in disease_list) {
    sample_file <- paste0("data/data_processed/", cohort_name, "/", sample, "-annotated_matrix.rds")
    sample_data <- as.data.frame(read_rds(sample_file))
    cell_types <- colnames(sample_data)
    cell_type_counts <- as.data.frame(table(cell_types)) |>
      rename(Cell_Type = cell_types, Cell_Count = Freq)
    disease_counts <- bind_rows(disease_counts, cell_type_counts)
  }
  disease_counts <- disease_counts %>%
    group_by(Cell_Type) %>%
    summarise(Cell_Count = sum(Cell_Count)) %>%
    mutate(Disorder = "Schizophrenia", Cohort = cohort_name)

  summary_df <- bind_rows(summary_df, control_counts, disease_counts)
  
  return(summary_df)
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
  
  colnames(ctmat) <- gsub("\\..*", "", colnames(ctmat))
  
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


main()
