#!/usr/bin/env Rscript

# Backward-compatible entry point for the adaptive-overlap model.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(dirname(script_file), mustWork = TRUE)

source(file.path(root, "scripts", "gf_workflow.R"))
run_gf_workflow(
  project_root = root,
  snp_set = "adaptive",
  output_dir = file.path(root, "results", "gradient_forest", "adaptive")
)
