# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script obtains summarized tables/figures of the data.
# Usage: R/summaryfigures.R

library(tidyverse)
library(ggplot2)
library(stringr)
library(cowplot)

main <- function() {
  clean_metadata <- read_csv("data/data_raw/raw_metadata.csv") |>
    filter(Cohort != "ROSMAP") |>
    mutate(Disorder = ifelse(tolower(Disorder) == "control", "Control", Disorder))

write.csv(clean_metadata, file.path("data/data_processed/clean_metadata.csv"), row.names = FALSE)
  
  total_patients <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort) |>
    summarize(
      Patients = n(),
      Mean_Age = round(mean(Age, na.rm = TRUE)),
      Disorder_studied = paste(unique(Disorder[Disorder != "Control"]), collapse = ", "),
      Disorder_studied = ifelse(Disorder_studied == "", "None", Disorder_studied)
    )
  
  write.csv(total_patients, file.path("results/1-total_patients.csv"), row.names = FALSE)
  
  patientpercondition <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort, Disorder) |>
    summarize(
      Number_of_Patients = n(),
      Mean_Age = round(mean(Age, na.rm = TRUE)),
      N_Male = sum(Biological_Sex == "male", na.rm = TRUE),
      N_Female = sum(Biological_Sex == "female", na.rm = TRUE)
    )
  
  
  write.csv(patientpercondition, file.path("results/2-patients_per_condition.csv"), row.names = FALSE)
  
  elderlypatients <- clean_metadata |>
    filter(Cohort == "CMC" | Cohort == "MultiomeBrain" | Cohort == "SZBDMulti-Seq") |>
    group_by(Cohort, Disorder) |>
    filter(Age_death == "89+") |>
    summarize(Plus89_patients = n())
  
  write.csv(elderlypatients, file.path("results/3-elderlypatients.csv"), row.names = FALSE)
  
  meanage_summary <- clean_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort, Disorder) |>
    summarize(Mean_Age = round(mean(Age)))
  
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
               linetype = "dashed", size = 1) +
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
               linetype = "dashed", size = 1) +
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
               linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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
             linetype = "dashed", size = 1) +
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

plot_list <- list(CMCageplot, DevBrainageplot, Girgentiageplot, ISOhubageplot,
                  LIBDageplot, Maageplot, MBageplot, PTSDageplot, SZBDageplot,
                  UCLAageplot, Velageplot)

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
}

main()
