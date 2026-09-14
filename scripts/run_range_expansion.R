#!/usr/bin/env Rscript

if (!requireNamespace("rangeExpansion", quietly = TRUE)) {
  stop(
    "Missing rangeExpansion. Install the pinned revision with: ",
    "remotes::install_github('BenjaminPeter/rangeExpansion@d485202')"
  )
}

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
genotype_path <- file.path(root, "180_ql.ld_prune.lfmm.csv")
metadata_path <- file.path(root, "data", "sample_metadata.tsv")
output_dir <- file.path(root, "results", "range_expansion")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genotype <- read.csv(
  genotype_path, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.table(
  metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  check.names = FALSE
)
required <- c("sample_id", "population", "longitude", "latitude")
missing_columns <- setdiff(required, names(metadata))
if (length(missing_columns)) {
  stop("Sample metadata is missing: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(rownames(genotype))) stop("Duplicate IDs in genotype matrix")
if (anyDuplicated(metadata$sample_id)) stop("Duplicate IDs in sample metadata")
if (!setequal(rownames(genotype), metadata$sample_id)) {
  stop("Genotype and metadata sample IDs do not match")
}
metadata <- metadata[match(rownames(genotype), metadata$sample_id), ]
if (anyNA(genotype)) stop("The range-expansion genotype matrix contains missing values")
if (!all(as.matrix(genotype) %in% 0:2)) stop("Genotypes must be coded 0, 1, or 2")

# load.data.snapp() expects headerless comma-separated rows beginning with ID.
snapp_path <- tempfile(fileext = ".snapp")
coordinate_path <- tempfile(fileext = ".csv")
on.exit(unlink(c(snapp_path, coordinate_path)), add = TRUE)
write.table(
  cbind(sample_id = rownames(genotype), genotype), snapp_path,
  sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE
)
coordinates <- data.frame(
  id = metadata$sample_id,
  latitude = metadata$latitude,
  longitude = metadata$longitude,
  region = "all",
  outgroup = 0,
  pop = metadata$population
)
write.table(
  coordinates, coordinate_path, sep = ",", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

raw_data <- rangeExpansion::load.data.snapp(
  snapp_path, coordinate_path, ploidy = 2, sep = ","
)

# The pinned package checks raw_data$coord rather than raw_data$coords before
# make.pop(). Supplying this compatibility field preserves named populations.
raw_data$coord <- list(pop = raw_data$coords$pop)
population_data <- rangeExpansion::make.pop(raw_data, ploidy = 2)
psi <- rangeExpansion::get.all.psi(population_data)
rownames(psi) <- population_data$coords$pop
colnames(psi) <- population_data$coords$pop
write.table(
  psi, file.path(output_dir, "pairwise_psi.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

grid_x <- as.integer(Sys.getenv("TDOA_XLEN", "20"))
grid_y <- as.integer(Sys.getenv("TDOA_YLEN", "20"))
if (anyNA(c(grid_x, grid_y)) || grid_x < 2 || grid_y < 2) {
  stop("TDOA_XLEN and TDOA_YLEN must be integers of at least 2")
}
old_directory <- getwd()
setwd(output_dir)
on.exit(setwd(old_directory), add = TRUE)
origin_results <- rangeExpansion::run.regions(
  region = list(NULL), pop = population_data, psi = psi,
  xlen = grid_x, ylen = grid_y
)
origin_summary <- summary(origin_results$tbl[[1]])
write.table(
  origin_summary, "origin_summary.tsv", sep = "\t",
  quote = FALSE, row.names = FALSE
)
saveRDS(origin_results, "origin_results.rds")

pdf("origin_map.pdf", width = 8, height = 6)
plot(origin_results, add.map = TRUE, add.samples = TRUE, add.sample.het = TRUE)
dev.off()
message("Wrote range-expansion outputs to ", output_dir)
