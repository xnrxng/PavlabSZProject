# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script downloads the metadata information and obtains summarized tables/figures of its information.
# Usage: R/loadmetadata_obtaintable.R

library(tidyverse)
library(data.table)

main <- function() {
  raw_metadata <- data.table::fread("https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt")
  write_csv(raw_metadata, file.path("data/raw_data"))
  
  total_patients <- raw_metadata |>
    mutate(Total_Patients = count(Cohort))
}

main()