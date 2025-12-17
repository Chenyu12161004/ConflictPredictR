## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

library(tidyverse)

full_data <- read_csv("~/Desktop/capstone project/ConflictPredictR/data-raw/clinvar_conflicting.csv")

set.seed(215)
sample_variants <- full_data %>%
  group_by(CLASS) %>%
  slice_sample(n = 5) %>%
  ungroup() %>%
  select(-CLASS)

usethis::use_data(sample_variants, overwrite = TRUE)
