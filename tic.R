# tic.R
library(tic)

# Define stages
drat_stage <- stage("drat")
pkgdown_stage <- stage("pkgdown")

# Define a custom stage for WASM build
wasm_stage <- stage("wasm") %>%
  add_step(step_script("./tools/build_wasm.R"))

# Define the overall pipeline
pipeline <- c(drat_stage, pkgdown_stage, wasm_stage)

# Run the pipeline on GitHub Actions
if (ci_on("github_actions")) {
  run_on = c("linux", "macos", "windows") # Adjust as needed
  get_stage("drat") %>%
    add_step(step_install_cran()) %>%
    add_step(step_build_pkg()) %>%
    add_step(step_check()) %>%
    add_step(step_drat(repo = "jgf5013/drat")) # Replace with your drat repo

  get_stage("pkgdown") %>%
    add_step(step_install_cran(c("pkgdown"))) %>%
    add_step(step_build_pkg()) %>%
    add_step(step_pkgdown())

  get_stage("wasm") %>%
    add_step(step_install_cran(c("r-wasm/wasm", "deSolve", "ggplot2", "rootSolve", "tidyr"))) %>%
    add_step(step_script("./tools/build_wasm.R"))

  run_pipeline()
}