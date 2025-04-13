library(tic)

# Define stages (keep other stages you might have)
# drat_stage <- stage("drat") # If you still have it
# pkgdown_stage <- stage("pkgdown") # If you still have it

wasm_stage <- stage("wasm") %>%
  add_step(step_install_cran(c("r-wasm/wasm", "deSolve", "ggplot2", "rootSolve", "tidyr"))) %>%
  add_step(step_script("./tools/build_wasm.R"))

# Define the overall pipeline
pipeline <- c(
  # drat_stage,
  # pkgdown_stage,
  wasm_stage
  # Add other stages in your desired order
)

# Run the pipeline on GitHub Actions
if (ci_on("github_actions")) {
  # Configuration for other stages (if you have them)
  # get_stage("drat") %>% ...
  # get_stage("pkgdown") %>% ...

  get_stage("wasm") %>%
    add_step(step_install_cran(c("r-wasm/wasm", "deSolve", "ggplot2", "rootSolve", "tidyr"))) %>%
    add_step(step_script("./tools/build_wasm.R"))

  run_pipeline()
}