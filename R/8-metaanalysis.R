# Author: Rui Xiang Yu
# Date: 2024 November 27th
# This script performs a metaanalysis.
# Usage: R/8-metaanalysis.R

library(tidyverse)
library(magrittr)
library(dplyr)
library(Matrix)
library(Seurat)
library(stringr)
library(metafor)

main <- function(){

  ### create a master list
  cohort_list <- c("CMC", "SZBDMulti-Seq")

  DGE = list()

  cell_types <- c("Ast", "Oli", "Opc", "Exc", "Inh", "Mic")
  for (cell_type in cell_types) {
    DGE[[cell_type]] = list() # Initialize only once per cell type
    for (cohort in cohort_list) {
      file_path <- paste0("results/DEA/", cohort, "/LimmaCPMLogAgeSex/DEAresults_",
                          cell_type, "_", cohort, "_SZ.rds")
      dge_data <- readRDS(file_path)
      DGE[[cell_type]][[cohort]] <- dge_data
    }
  }

  saveRDS(DGE, "results/DEA/RuzickaMetaanalysis/DGE_meta.rds")

  ### analysis ======================================================================

  combined.analysis.tables = list()
  for (cell_type in names(DGE)) {
    df_list = DGE[[cell_type]]

    ### 1. filter the genes
    common_genes <-  Reduce(intersect, lapply(df_list, function(df) rownames(df)))
    df_list_filtered <- lapply(df_list, filter_dataframe, common_genes)

    ### 2. calculate se
    tbls = lapply( df_list_filtered, function(tab){
      tab$se = tab$groupdisorderyes / tab$F
      tab
    })

    ### 3. calculate gene tables
    gene.tbls = lapply(1:nrow(tbls[[1]]), function(i) {
      dfs = lapply(1:length(tbls), function(k) tbls[[k]][i, ])
      df = do.call("rbind", dfs)
    })
    names(gene.tbls) = tbls[[1]]$gene


    ### 4 combined analysis table
    combined.analysis.tbl = do.call(rbind, lapply(names(gene.tbls), function(gene){
      x = suppressWarnings(metafor::rma(yi=groupdisorderyes, sei=se, data = gene.tbls[[gene]], method="FE"))
      combined.tbl = data.frame( gene = gene,
                                 logFC     = x$beta,
                                 se        = x$se,
                                 tstat = x$zval,
                                 P.Value   = x$pval)
      return(combined.tbl)
    }))
    rownames(combined.analysis.tbl) = names(gene.tbls)
    combined.analysis.tbl = combined.analysis.tbl[order(combined.analysis.tbl$P.Value), ]

    combined.analysis.tables[[cell_type]] = combined.analysis.tbl
  }



  names(combined.analysis.tables) == names(DGE)

  ### local FDR adjustment
  combined.analysis.tables.v2 = lapply(combined.analysis.tables, function(DF) {
    DF$adj.P.Val = p.adjust(DF$P.Value, method = "fdr")
    #DF$adj.P.Val.BH = p.adjust(DF$P.Value, method = "BH")
    return(DF)
  })

  ### global FDR adjustment
  DF.total = do.call(rbind, combined.analysis.tables.v2)
  DF.total$adj.P.Val.global = p.adjust(DF.total$P.Value, method = "fdr")

  ### split data back to cell types
  ff = factor(unlist(lapply(names(combined.analysis.tables.v2), function(celltype) rep(celltype, nrow(combined.analysis.tables.v2[[celltype]])))), names(DGE))
  combined.analysis.tables.v3 = split(DF.total, ff)

  saveRDS(combined.analysis.tables.v3, "results/DEA/RuzickaMetaanalysis/DEAresults_metaanalysis.rds")

  ### papers thresholds
  print(sort(sapply(combined.analysis.tables.v3, function(DF) sum( (DF$adj.P.Val < 0.05) & (abs(DF$logFC) > 0.1) )), decreasing = T)[1:15])
  print(sort(sapply(combined.analysis.tables.v3, function(DF) sum( (DF$adj.P.Val.global < 0.05) & (abs(DF$logFC) > 0.1))), decreasing = T)[1:15])

  ### my thresholds
  print(sort(sapply(combined.analysis.tables.v3, function(DF) sum( (DF$adj.P.Val < 0.05) & (abs(DF$logFC) > 0.5) )), decreasing = T)[1:15])
  print(sort(sapply(combined.analysis.tables.v3, function(DF) sum( (DF$adj.P.Val.global < 0.05) & (abs(DF$logFC) > 0.5))), decreasing = T)[1:15])


  ### visualization
  plot_list <- list()
  for (cell_type in names(combined.analysis.tables.v3)){
    data <- combined.analysis.tables.v3[[cell_type]]
      res_his <- data |>
        ggplot(aes(x = logFC)) +
        geom_histogram() +
        theme_classic()+
        labs(
          title = cell_type,
          x = "logFC",
          y = NULL
        )

    plot_list[[cell_type]] <- res_his}

  all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

  title_plot <- ggdraw() +
    draw_label("CMC & SZBDMulti-Seq MetaAnalysis", x = 0.5, y = 0.5, size = 12, hjust = 0.5)

  all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

  final_path <- paste0("results/DEA/RuzickaMetaanalysis/logFCplot_ruzicka.png")
  cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

plot_list <- list()
for (cell_type in names(combined.analysis.tables.v3)){
  data <- combined.analysis.tables.v3[[cell_type]]
  res_his <- data |>
    ggplot(aes(x = P.Value)) +
    geom_histogram() +
    theme_classic()+
    labs(
      title = cell_type,
      x = "P-value",
      y = NULL
    )

  plot_list[[cell_type]] <- res_his}

all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

title_plot <- ggdraw() +
  draw_label("CMC & SZBDMulti-Seq MetaAnalysis", x = 0.5, y = 0.5, size = 12, hjust = 0.5)

all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

final_path <- paste0("results/DEA/RuzickaMetaanalysis/pvalplot_ruzicka.png")
cowplot::save_plot(plot = all_plots_with_title, filename = final_path)

plot_list <- list()
for (cell_type in names(combined.analysis.tables.v3)){
  data <- combined.analysis.tables.v3[[cell_type]]
  res_volcano <- EnhancedVolcano(data, lab = data$gene,
                                 x = 'logFC', y = 'P.Value', title = cell_type, subtitle = NULL, FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "none", caption = NULL,
                                 axisLabSize = 10, titleLabSize = 10, labSize = 3, pointSize = 1)

  plot_list[[cell_type]] <- res_volcano}

all_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")

title_plot <- ggdraw() +
  draw_label("CMC & SZBDMulti-Seq MetaAnalysis", x = 0.5, y = 0.5, size = 12, hjust = 0.5)

all_plots_with_title <- plot_grid(title_plot, all_plots, ncol = 1, rel_heights = c(0.1, 1))

final_path <- paste0("results/DEA/RuzickaMetaanalysis/volcanoplot_ruzicka.png")
ggsave(plot = all_plots_with_title, filename = final_path, height = 9, width = 10)

  }

### helper functions
  filter_dataframe <- function(df, common_genes) {
    df |>
      rownames_to_column("gene") |>
      filter(gene %in% common_genes) |>
      mutate(gene = factor(gene, levels = common_genes)) |>
      arrange(gene)
  }

main()
