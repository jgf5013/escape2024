################################################################################
#
# Table 1
#
################################################################################

library(dplyr)
library(ggplot2)
library(escape2024)

fs::dir_create("results")


################################################################################
# Final size based on LSHTM estimates for R_0
# mass action reference model
# table

# read data
if (file.exists("data/pathogen_R0_values.csv")) {
  df_reproduction_numbers <- readr::read_csv("data/pathogen_R0_values.csv"
  )
} else {
  warning("data file not found, dummy dataset is used")
  df_reproduction_numbers <- tibble(
    Pathogen = "test",
    Archetype = 1,
    R0 = 1.5
  )
}

# calculate mass action final size for each R_0 value
df_reproduction_numbers <- df_reproduction_numbers |>
  rowwise() |>
  mutate(
    final_size_theory = theory_reference_final_size(r_0 = `Basic reproduction number`),
    final_size_theory = final_size_theory |> round(2)
  ) |>
  select(Pathogen,
    "Basic reproduction number",
    "Mass action final size" = final_size_theory
  )

# write table to disk
df_reproduction_numbers |>
  readr::write_csv("results/mass_action_final_size.csv")
