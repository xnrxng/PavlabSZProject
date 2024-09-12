# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script obtains summarized tables/figures of the data.
# Usage: R/summaryfigures.R

library(tidyverse)
library(stringr)

main <- function() {
  raw_metadata <- read_csv("data/data_raw/raw_metadata.csv")
  
  total_patients <- raw_metadata |>
    mutate(Age = as.numeric(str_replace(Age_death, "\\+", ""))) |>
    group_by(Cohort) |>
    summarize(
      Total_Patients = n(),
      Mean_Age_Death = round(mean(Age, na.rm = TRUE)),
      Disease_studied = paste(unique(Disorder[Disorder != "control" & Disorder != "Control"]), collapse = ", "),
      Disease_studied = ifelse(Disease_studied == "", "None", Disease_studied)
    )
  
  patientpercondition <- raw_metadata |>
    group_by(Cohort, Disorder) |>
    summarize(Number_of_Patients = n())
}

main()
