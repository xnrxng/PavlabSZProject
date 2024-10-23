# Author: Rui Xiang Yu
# Date: 2024 October 21st
# This script obtains the number of single cells and genes in the Schizophrenia cohorts.
# Usage: R/4-szagedistributionandumap.R

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

main <- function() {
  CMC_meta <- readRDS("data/data_processed/CMC/CMC-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  MB_meta <- readRDS("data/data_processed/MultiomeBrain/MultiomeBrain-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  SZBD_meta <- readRDS("data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patients.rds") |>
    mutate(age = as.numeric(str_replace(age, "\\+", "")))
  
  CMC_age_hist <- CMC_meta |>
    mutate(disorder = recode(disorder, "yes" = "Schizophrenia", "no" = "Control")) |>
    ggplot(aes(x = age)) +
    geom_histogram(fill = "blue4") +
    labs(x = "Age at death (years)",
         y = "Number of patients",
         title = "CMC") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))
  
  print(CMC_age_hist)
  
  
  
}

main()
