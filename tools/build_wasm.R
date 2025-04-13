#!/usr/bin/env Rscript

# tools/build_wasm.R

# Check if the 'wasm' package is available
if (requireNamespace("wasm", quietly = TRUE)) {
  message("Attempting to build WASM binary...")
  try({
    message("Installing required packages...")
    install.packages(c("deSolve", "ggplot2", "rootSolve", "tidyr"), repos = "https://cloud.r-project.org/")
    message("Packages installed successfully.")
    message("Building WASM binary...")
    # Assuming your package source is at the root of the repository
    # and you want to place the WASM artifacts in 'inst/wasm'
    wasm::build_wasm(
      pkg = "./R",
      dest = "inst/wasm",
      export = "*" # Or specify functions to export if needed
      # You might need to adjust other parameters based on your package
    )
    message("WASM build completed. Files saved in inst/wasm/")
  }, error = function(e) {
    message("Error during WASM build:")
    message(e)
  })
} else {
  message("Warning: 'wasm' package is not available. Skipping WASM build.")
}