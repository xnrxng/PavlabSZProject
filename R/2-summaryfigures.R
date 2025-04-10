# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script obtains summarized tables/figures of the data.
# Usage: R/2-summaryfigures.R

library(tidyverse)
library(ggplot2)
library(stringr)
library(cowplot)

main <- function() {
  ### get total patients per cohort
  clean_metadata <- readRDS("data/data_processed/clean_metadata.rds")
  batiuk_meta <- readRDS("data/data_processed/Batiuk/Batiuk-patient.rds")

  total_patients <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort) |>
    summarize(
      Patients = n(),
      Mean_Age = round(mean(Age, na.rm = TRUE)),
      Disorder_studied = paste(unique(Disorder[Disorder != "Control"]), collapse = ", "),
      Disorder_studied = ifelse(Disorder_studied == "", "None", Disorder_studied)
    )

  batiuk_total <- data.frame(
    Cohort = "Batiuk",
    Patients = length(batiuk_meta$patientID),
    Mean_Age = mean(as.numeric(batiuk_meta$age)),
    Disorder_studied = "Schizophrenia",
    stringsAsFactors = FALSE
  )

  all_patients <- rbind(total_patients, batiuk_total)

  saveRDS(all_patients, file.path("results/1-total_patients.rds"))

  ### get number of patients per condition per cohort
  patientpercondition <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort, Disorder) |>
    summarize(
      Total = n(),
      Mean_Age = round(mean(Age, na.rm = TRUE)),
      N_Male = sum(Biological_Sex == "male", na.rm = TRUE),
      N_Female = sum(Biological_Sex == "female", na.rm = TRUE)
    )

  batiuk_per_condition <- batiuk_meta |>
    group_by(disorder) |>
    summarize(
      Total = n(),
      Mean_Age = round(mean(as.numeric(age, na.rm = TRUE))),
      N_Male = sum(sex == "male", na.rm = TRUE),
      N_Female = sum(sex == "female", na.rm = TRUE)
    ) |>
    rename(Disorder = disorder) |>
    mutate(Cohort = "Batiuk",
           Disorder = ifelse(Disorder == "yes", "Schizophrenia", "Control"))

  allpatients_per_condition <- rbind(patientpercondition, batiuk_per_condition)

  saveRDS(allpatients_per_condition, file.path("results/2-patients_per_condition.rds"))

  ### get number of elderly patients in schizophrenia cohorts
  elderlypatients <- clean_metadata |>
    filter(Cohort == "CMC" | Cohort == "MultiomeBrain" | Cohort == "SZBDMulti-Seq") |>
    group_by(Cohort, Disorder) |>
    filter(Age_death == "89+") |>
    summarize(Plus89_patients = n())

  saveRDS(elderlypatients, file.path("results/3-elderlypatients.rds"))

  ### make a plot of age distributions
  meanage_summary <- allpatients_per_condition |>
    select(Cohort, Disorder, Mean_Age)

  x_limits <- c(0, 90)
  y_limits <- c(0, 0.2)

color_palette <- c(
    "Control" = "turquoise",
    "Schizophrenia" = "coral",
    "ASD" = "pink2",
    "Bipolar Disorder" = "yellow2",
    "MDD" = "olivedrab3",
    "PTSD" = "navy",
    "Williams Syndrome" = "red3"
  )

  CMCageplot <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "CMC") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "CMC"),
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", linewidth = 1) +
    labs(title = "CMC") +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_bw(base_size = 10, base_line_size  = 1) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

  DevBrainageplot <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "DevBrain") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "DevBrain"),
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", linewidth = 1) +
    labs(title = "DevBrain") +
    scale_color_manual(values = color_palette) +
    scale_fill_manual(values = color_palette) +
    theme_bw(base_size = 10, base_line_size  = 1) +
    xlim(x_limits) +
    ylim(y_limits) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

Girgentiageplot <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    filter(Cohort == "Girgenti-snMultiome") |>
    ggplot(aes(x = Age, fill = Disorder)) +
    geom_density(alpha = 0.3) +
    geom_vline(data = meanage_summary |> filter(Cohort == "Girgenti-snMultiome"),
               aes(xintercept = Mean_Age, colour = Disorder),
               linetype = "dashed", linewidth = 1) +
    labs(title = "Girgenti-snMultiome") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
    theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )

ISOhubageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "IsoHuB") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "IsoHuB"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "IsoHuB") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

LIBDageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "LIBD") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "LIBD"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "LIBD") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

Maageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "Ma_et_al") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "Ma_et_al"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Ma-Sestan") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

MBageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "MultiomeBrain") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "MultiomeBrain"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "MultiomeBrain") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

PTSDageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "PTSDBrainomics") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "PTSDBrainomics"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "PTSDBrainomics") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

SZBDageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "SZBDMulti-Seq") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "SZBDMulti-Seq"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "SZBDMulti-Seq") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

UCLAageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "UCLA-ASD") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "UCLA-ASD"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "UCLA-ASD") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

Velageplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  filter(Cohort == "Velmeshev_et_al") |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "Velmeshev_et_al"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Velmeshev_et_al") +
  scale_color_manual(values = color_palette, guide = "none") +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

Batiukageplot <- batiuk_meta |>
  mutate(disorder = ifelse(disorder == "yes", "Schizophrenia", "Control"),
         age = as.numeric(age)) |>
  ggplot(aes(x = age, fill = disorder)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meanage_summary |> filter(Cohort == "Batiuk"),
             aes(xintercept = Mean_Age, colour = Disorder),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Batiuk") +
  scale_color_manual(values = color_palette, guide = "none") +
  scale_fill_manual(values = color_palette) +
  theme_bw(base_size = 10, base_line_size  = 1) +
  xlim(x_limits) +
  ylim(y_limits) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

plot_list <- list(CMCageplot, DevBrainageplot, Girgentiageplot, Batiukageplot, ISOhubageplot,
                  LIBDageplot, SZBDageplot, Maageplot, PTSDageplot,
                  MBageplot, UCLAageplot, Velageplot)

combined_plots <- plot_grid(
  plotlist = plot_list,
  ncol = 3,
  align = "hv"
)

dummyplot <- clean_metadata |>
  mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
  ggplot(aes(x = Age, fill = Disorder)) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  theme(legend.box.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size = 10))

legend <- get_legend(dummyplot)

x_label_plot <- ggdraw() +
  draw_label("Age at Death (years)", x = 0.5, y = 0.5, size = 14)

y_label_plot <- ggdraw() +
  draw_label("Density", x = 0.5, y = 0.5, angle = 90, size = 14)

final_plot <- plot_grid(
  plot_grid(
    y_label_plot,
    combined_plots,
    ncol = 2,
    rel_widths = c(0.1, 1)
  ),
  x_label_plot,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

final_plot_with_legend <- plot_grid(
  final_plot,
  legend,
  ncol = 2,
  rel_widths = c(0.85, 0.15)
)

ggsave(file.path("results/4-agedistribution.png"), final_plot_with_legend, width = 12,
       height = 8, units = "in", dpi = 300)

###
patientsexdisorder_data <- clean_metadata |>
  filter(Cohort == 'CMC' | Cohort == "MultiomeBrain" | Cohort == "SZBDMulti-Seq") |>
  filter(Disorder == "Control" | Disorder == "Schizophrenia") |>
  group_by(Cohort, Disorder, Biological_Sex) |>
  summarize(count = n())

batiuksexdisorder_data <- batiuk_meta |>
  group_by(disorder, sex) |>
  summarize(count = n()) |>
  rename(Disorder = disorder, Biological_Sex = sex) |>
  mutate(Cohort = "Batiuk",
         Disorder = ifelse(Disorder == "yes", "Schizophrenia", "Control"))

allpatientsexdisorder_data <- rbind(patientsexdisorder_data, batiuksexdisorder_data)

patient_disorderbarplot <- allpatientsexdisorder_data |>
  ggplot(aes(x = Disorder, y = count)) +
  geom_bar(aes(fill = Biological_Sex), stat = "identity") +
  facet_wrap(~ Cohort) +
  scale_fill_manual(values = c("male" = "hotpink1", "female" = "mediumseagreen")) +
  labs(y = "Number of patients", fill = "Biological \n sex") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60))

ggsave(file.path("results/12-patientbarplot.png"), patient_disorderbarplot, width = 8,
       height = 8, units = "in", dpi = 300)
}

main()
