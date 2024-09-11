# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script obtains summarized tables/figures of the data.
# Usage: R/obtaintable.R

library(tidyverse)

main <- function() {
  raw_metadata <- read_csv("data/data_raw/raw_metadata.csv")
  total_patients <- raw_metadata |>
    unique(Individual_ID) |>
    mutate(Total_Patients = count(Cohort))
}

main()
