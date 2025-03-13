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

}

main()
