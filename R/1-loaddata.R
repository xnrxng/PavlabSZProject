# Author: Rui Xiang Yu
# Date: 2024 September 11th
# This script downloads the data To be expanded.
# Usage: R/loadmetadata_obtaintable.R

library(tidyverse)
library(data.table)

main <- function() {
  raw_metadata <- data.table::fread("https://brainscope.gersteinlab.org/data/sample_metadata/PEC2_sample_metadata.txt")
  write_csv(raw_metadata, file.path("data/data_raw/raw_metadata.csv"))
}

main()