#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("adaptive", "all")) {
  stop("Usage: Rscript scripts/run_gradient_forest.R {adaptive|all}", call. = FALSE)
}
snp_set <- args[1]
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

source(file.path(root, "scripts", "gf_workflow.R"))
run_gf_workflow(
  project_root = root,
  snp_set = snp_set,
  output_dir = file.path(root, "results", "gradient_forest", snp_set)
)
