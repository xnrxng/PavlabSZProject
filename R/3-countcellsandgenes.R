# Author: Rui Xiang Yu
# Date: 2024 September 19th
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/3-countcellsandgenes.R

library(tidyverse)
library(ggplot2)
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
  
  ### filtered cohort numbers
  # results_filtered <- data.frame(
  #   Cohort = character(),
  #   Genes = integer(),
  #   Total_Cells = integer(),
  #   stringsAsFactors = FALSE
  # )
  # 
  # results_filtered <- process_cohort_filtered("CMC", results_filtered)
  # results_filtered <- process_cohort_filtered("SZBDMulti-Seq", results_filtered)
  # results_filtered <- process_cohort_filtered("MultiomeBrain", results_filtered)
  # 
  # saveRDS(results_filtered, file.path("results/5.2-cells_genes.rds"))

  ### cell type plot
  celltype_summary_CMC <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )
  
  celltype_summary_SZBD <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )
  
  celltype_summary_MB <- data.frame(
    patientID = character(),
    cells = integer(),
    disorder = character(),
    cell_type = character(),
    stringsAsFactors = FALSE
  )
  
  celltype_summary_CMC <- countcells_pertype("CMC", celltype_summary_CMC)
  celltype_summary_SZBD <- countcells_pertype("SZBDMulti-Seq", celltype_summary_SZBD)
  celltype_summary_MB <- countcells_pertype("MultiomeBrain", celltype_summary_MB)

  CMCcelltypeplot <- celltype_summary_CMC |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = cells, y = cell_type, color = disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.8) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "CMC") +
    xlim(c(0, 12000)) +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2"))
  
  SZBDcelltypeplot <- celltype_summary_SZBD |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = cells, y = cell_type, color = disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.8) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "SZBDMulti-Seq") +
    xlim(c(0, 12000)) +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2"))
  
  MBcelltypeplot <- celltype_summary_MB |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = cells, y = cell_type, color = disorder)) +
    geom_point(position = position_jitter(width = 0.2, height = 0.1), size = 1, alpha = 0.8) +
    labs(
      x = "Number of Cells (Frequency)", 
      y = "Cell Type",
      title = "MultiomeBrain") +
    xlim(c(0, 12000)) +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2"))
  
  celltypecombined_plots <- plot_grid(
    plotlist = list(CMCcelltypeplot, SZBDcelltypeplot, MBcelltypeplot+ theme(legend.position = "none")), ncol = 1, align = "hv")
  
  disorderlegend <- get_legend(MBcelltypeplot+ labs(color = "Disorder") + theme(legend.box.margin = margin(0, 0, 0, 3),
                                             legend.key.size = unit(1, 'cm'),
                                             legend.title = element_text(size = 15)))
  
  x_label_plot <- ggdraw() + 
    draw_label("Number of Cells (Frequency)", x = 0.5, y = 0.5, size = 14)
  
  celltypefinal_plot <- plot_grid(celltypecombined_plots, x_label_plot, ncol = 1, rel_heights = c(1, 0.1))
  
  celltypefinal_plot_with_legend <- plot_grid(
    celltypefinal_plot,
    disorderlegend,
    ncol = 2,  
    rel_widths = c(0.85, 0.15)
  )
  
  ggsave(file.path("results/6-celltypedistribution.png"), celltypefinal_plot_with_legend, width = 10, height = 12)
  
  ###sequencing depth plot
  depth_df_CMC <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_SZBD <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )

  depth_df_MB <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_CMC <- sum_columns("CMC", depth_df_CMC)
  depth_df_SZBD <- sum_columns("SZBDMulti-Seq", depth_df_SZBD)
  depth_df_MB <- sum_columns("MultiomeBrain", depth_df_MB)
  combined_depthdf <- rbind(depth_df_CMC, depth_df_SZBD, depth_df_MB)
  
  depthplot <- combined_depthdf |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = seq_depth, y = cell_type, color = disorder)) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = disorder), size=0.4, alpha=0.5) +
    xlim(0, 300000) +
    labs(x = "Sequencing Depth (number  of reads)",
         y = "Cell type",
         color = "Disorder") +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    facet_wrap(~ cohort, ncol =1, scales = "free_x")
      
  ggsave(file.path("results/7-depthdistribution.png"), depthplot, width = 8, height = 15)
  
  depthplot_unscaled <- combined_depthdf |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = seq_depth, y = cell_type, color = disorder)) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = disorder), size=0.4, alpha=0.5) +
    labs(x = "Sequencing Depth (number  of reads)",
         y = "Cell type",
         color = "Disorder") +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    facet_wrap(~ cohort, ncol =1, scales = "free_x")
  
  ggsave(file.path("results/7.2-depthdistribution_unscaled.png"), depthplot_unscaled, width = 8, height = 15)
  
  ### sequencing depth plot for CPM samples
  depth_df_CMC_CPM <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_SZBD_CPM <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_MB_CPM <- data.frame(
    patientID = character(),
    cell_type = character(),
    seq_depth = integer(),
    disorder = character(),
    cohort = character(),
    stringsAsFactors = FALSE
  )
  
  depth_df_CMC_CPM <- sum_columns_cpm("CMC", depth_df_CMC_CPM)
  depth_df_SZBD_CPM <- sum_columns_cpm("SZBDMulti-Seq", depth_df_SZBD_CPM)
  depth_df_MB_CPM <- sum_columns_cpm("MultiomeBrain", depth_df_MB_CPM)
  combined_depthdf_cpm <- rbind(depth_df_CMC_CPM, depth_df_SZBD_CPM, depth_df_MB_CPM)
  
  
  depthplot_cpm <- combined_depthdf |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = seq_depth, y = cell_type, color = disorder)) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = disorder), size=0.4, alpha=0.5) +
    xlim(0, 300000) +
    labs(x = "Sequencing Depth (number  of reads)",
         y = "Cell type",
         color = "Disorder") +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    facet_wrap(~ cohort, ncol =1, scales = "free_x")
  
  ggsave(file.path("results/7.3-depthdistribution_cpm.png"), depthplot_cpm, width = 8, height = 15)
  
  depthplot_cpm_unscaled <- combined_depthdf |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = seq_depth, y = cell_type, color = disorder)) +
    geom_boxplot(alpha = 0.5, color = "black", outlier.shape = NA) +
    geom_jitter(aes(color = disorder), size=0.4, alpha=0.5) +
    labs(x = "Sequencing Depth (number  of reads)",
         y = "Cell type",
         color = "Disorder") +
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
    scale_color_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2")) +
    facet_wrap(~ cohort, ncol =1, scales = "free_x")
  
  ggsave(file.path("results/7.4-depthdistribution_cpm_unscaled.png"), depthplot_cpm_unscaled, width = 8, height = 15)
  
  ### cell type abundance per condition plot
  abundance_df <- data.frame(cohort = character(),
                             disorder = character(),
                             cell_type = character(),
                             cells = integer(),
                             stringsAsFactors = FALSE)

  abundance_df <- abundance_conditions("CMC", abundance_df)
  abundance_df <- abundance_conditions("SZBDMulti-Seq", abundance_df)
  abundance_df <- abundance_conditions("MultiomeBrain", abundance_df)

  abundance_plot <- abundance_df |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = cell_type, y = cells, fill = disorder)) +
    geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ cohort, ncol=1, scales = "free_x") +
  labs(x = "Cell Type",
       y = "Cell Count",
       fill = "Disorder") +
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
    scale_fill_manual(values = c("Control" = "grey21", "Schizophrenia" = "firebrick2"))
  
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

# process_cohort_filtered <- function(cohort, results_df) {
#   astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
#   excitatory <- paste0("Exc_", cohort, "_SZ.rds")
#   inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
#   microglia <- paste0("Mic_", cohort, "_SZ.rds")
#   oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
#   opc <- paste0("Opc_", cohort, "_SZ.rds")
#   
#   cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
#   
#   total_cells <- 0
#   genes <- 0
#   for (i in seq_along(cell_list)) {
#     sample <- cell_list[i]
#     sample_file <- paste0("data/data_processed/", cohort, "/FilteredV1/", sample)
#     sample_data <- readRDS(sample_file)
#     total_cells <- total_cells + ncol(sample_data$expr)
#     
#     if (i == 1) {
#       genes <- nrow(sample_data$expr)
#     }
#   }
#   
#   results_df <- rbind(results_df, data.frame(
#     Cohort = cohort,
#     Genes = genes,
#     Total_Cells = total_cells,
#     stringsAsFactors = FALSE
#   ))
#   
#   return(results_df)
# }

countcells_pertype <- function(cohort, results_df) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")
  
  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)

  for (cell_type_file in cell_list) {
    cell_type <- strsplit(cell_type_file, "_")[[1]][1]
    
    sample_file <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell_type_file)
    sample_data <- readRDS(sample_file)
    
    n_cells <- sample_data$meta |>
      group_by(patientID, disorder) |>
      summarize(cells = n()) |>
      mutate(cell_type = cell_type) |>
      select(patientID, cells, disorder, cell_type)
    results_df <- bind_rows(results_df, n_cells)
  }
  
  return(results_df)
}

sum_columns <- function(cohort, results_df) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")
  
  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
  
  for (cell_type_file in cell_list) {
    cell_type <- strsplit(cell_type_file, "_")[[1]][1]
    
    sample_file <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell_type_file)
    sample_data <- readRDS(sample_file)
    
    expr_matrix <- sample_data$expr
    metadata <- sample_data$meta
  
    seq_depth <- colSums(expr_matrix)
    
    seq_depth_df <- data.frame(
      cell_type = cell_type,
      seq_depth = as.integer(seq_depth),
      patientID = metadata$patientID[match(colnames(expr_matrix), rownames(metadata))],
      disorder = metadata$disorder[match(colnames(expr_matrix), rownames(metadata))],
      cohort = cohort,
      stringsAsFactors = FALSE)
    results_df <- bind_rows(results_df, seq_depth_df)
  }
  return(results_df)
}

sum_columns_cpm <- function(cohort, results_df) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")
  
  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
  
  for (cell_type_file in cell_list) {
    cell_type <- strsplit(cell_type_file, "_")[[1]][1]
    
    sample_file <- paste0("data/data_processed/", cohort, "/SingleCell/CPM/", cell_type_file)
    sample_data <- readRDS(sample_file)
    
    expr_matrix <- sample_data$expr
    metadata <- sample_data$meta
    
    seq_depth <- colSums(expr_matrix)
    
    seq_depth_df <- data.frame(
      cell_type = cell_type,
      seq_depth = as.integer(seq_depth),
      patientID = metadata$patientID[match(colnames(expr_matrix), rownames(metadata))],
      disorder = metadata$disorder[match(colnames(expr_matrix), rownames(metadata))],
      cohort = cohort,
      stringsAsFactors = FALSE)
    results_df <- bind_rows(results_df, seq_depth_df)
  }
  return(results_df)
}

abundance_conditions <- function(cohort, summary_df) {
  astrocyte <- paste0("Ast_", cohort, "_SZ.rds")
  excitatory <- paste0("Exc_", cohort, "_SZ.rds")
  inhibitory <- paste0("Inh_", cohort, "_SZ.rds")
  microglia <- paste0("Mic_", cohort, "_SZ.rds")
  oligodendrocyte <- paste0("Oli_", cohort, "_SZ.rds")
  opc <- paste0("Opc_", cohort, "_SZ.rds")
  
  cell_list <- c(astrocyte, excitatory, inhibitory, microglia, oligodendrocyte, opc)
  
  for (cell_type_file in cell_list) {
    cell_type <- strsplit(cell_type_file, "_")[[1]][1]
    
    sample_file <- paste0("data/data_processed/", cohort, "/FilteredV1/", cell_type_file)
    sample_data <- readRDS(sample_file)
    
    metadata <- sample_data$meta
    
    cells_per_condition <- metadata |>
      group_by(disorder) |>
      summarize(cells = n()) |>
      mutate(cell_type = cell_type, cohort = cohort) |>
      select(cells, disorder, cell_type, cohort)
    summary_df <- bind_rows(summary_df, cells_per_condition)
  }
  
  return(summary_df)
}

main()
