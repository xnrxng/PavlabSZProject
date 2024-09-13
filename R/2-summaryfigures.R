# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script obtains summarized tables/figures of the data.
# Usage: R/summaryfigures.R

library(tidyverse)
library(ggplot2)
library(stringr)
library(RColorBrewer)

main <- function() {
  raw_metadata <- read_csv("data/data_raw/raw_metadata.csv") |>
    filter(Cohort != "ROSMAP")
  
  total_patients <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort) |>
    summarize(
      Total_Patients = n(),
      Mean_Age_Death = round(mean(Age, na.rm = TRUE)),
      Disorder_studied = paste(unique(Disorder[Disorder != "control" & Disorder != "Control"]), collapse = ", "),
      Disorder_studied = ifelse(Disorder_studied == "", "None", Disorder_studied)
    )
  
  write.csv(total_patients, file.path("results/1-total_patients.csv"))
  
  patientpercondition <- raw_metadata |>
    group_by(Cohort, Disorder) |>
    summarize(Number_of_Patients = n())
  
  write.csv(patientpercondition, file.path("results/2-patients_per_condition.csv"))
  
  elderlypatients <- raw_metadata |>
    filter(Cohort == "CMC" | Cohort == "MultiomeBrain" | Cohort == "SZBDMulti-Seq") |>
    group_by(Cohort, Disorder) |>
    filter(Age_death == "89+") |>
    summarize(Plus89_patients = n())
  
  write.csv(elderlypatients, file.path("results/3-elderlypatients.csv"))
  
  meanage_summary <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort, Disorder) |>
    summarize(Mean_Age = round(mean(Age)))
  
  CMCageplot <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "CMC") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "CMC"), 
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", size = 1) +
    labs(x = "Age at Death (years)",
         y = "Density",
         fill = "Disorder", 
         colour = "Mean of Age") +
    scale_color_manual(labels = c("Control", 
                                  "Schizophrenia"),
                       values = c("darkseagreen", "coral")) +
    theme_bw(base_size = 20, base_line_size  = 1) +
    scale_fill_brewer(palette = "Set2") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
  
  DevBrainageplot <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "DevBrain") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "DevBrain"), 
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", size = 1) +
    labs(x = "Age at Death (years)",
         y = "Density",
         fill = "Disorder", 
         colour = "Mean of Age") +
    scale_color_manual(labels = c("Autism SD", 
                                  "Control", "Williams Syndrome"),
                       values = c("darkseagreen", "coral", "purple")) +
    theme_bw(base_size = 20, base_line_size  = 1) +
    scale_fill_brewer(palette = "Set2") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
  
Girgentiageplot <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "Girgenti-snMultiome") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "Girgenti-snMultiome"), 
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", size = 1) +
    labs(x = "Age at Death (years)",
         y = "Density",
         fill = "Disorder", 
         colour = "Mean of Age") +
    scale_color_manual(labels = c("Control"),
                       values = c("darkseagreen")) +
    theme_bw(base_size = 20, base_line_size  = 1) +
    scale_fill_brewer(palette = "Set2") +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )

ISOhubageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "IsoHuB") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "IsoHuB"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Control"),
                     values = c("darkseagreen")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

LIBDageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "LIBD") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "LIBD"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Control"),
                     values = c("darkseagreen")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )
  
Maageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "Ma_et_al") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "Ma_et_al"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Control"),
                     values = c("darkseagreen")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

MBageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "MultiomeBrain") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "MultiomeBrain"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Bipolar disorder", 
                                "Control", "Schizophrenia"),
                     values = c("darkseagreen", "coral", "purple")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

PTSDageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "PTSDBrainomics") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "PTSDBrainomics"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Control", 
                                "MDD", "PTSD"),
                     values = c("darkseagreen", "coral", "purple")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

SZBDageplot <- raw_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "SZBDMulti-Seq") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "SZBDMulti-seq"), 
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", size = 1) +
  labs(x = "Age at Death (years)",
       y = "Density",
       fill = "Disorder", 
       colour = "Mean of Age") +
  scale_color_manual(labels = c("Bipolar disorder", 
                                "Control", "Schizophrenia"),
                     values = c("darkseagreen", "coral", "purple")) +
  theme_bw(base_size = 20, base_line_size  = 1) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

  print(SZBDageplot)
  
  plot_grid(CMCage_plot, ncol = 2)
  
  ggsave(file.path("results/4-agedistribution"), histogram)
}

main()
