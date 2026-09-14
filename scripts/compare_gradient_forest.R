#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
adaptive_path <- file.path(root, "results", "gradient_forest", "adaptive", "gf_prediction_pca_scores.tsv")
all_path <- file.path(root, "results", "gradient_forest", "all", "gf_prediction_pca_scores.tsv")
for (path in c(adaptive_path, all_path)) {
  if (!file.exists(path)) stop("Missing model output: ", path)
}

adaptive <- read.table(adaptive_path, header = TRUE, sep = "\t")
all_snps <- read.table(all_path, header = TRUE, sep = "\t")
key <- c("longitude", "latitude", "period")
joined <- merge(adaptive, all_snps, by = key, suffixes = c("_adaptive", "_all"))
if (!nrow(joined)) stop("The two prediction tables do not share grid coordinates")

correlations <- data.frame(
  axis = paste0("PC", 1:3),
  correlation = vapply(1:3, function(axis) {
    cor(joined[[paste0("PC", axis, "_adaptive")]], joined[[paste0("PC", axis, "_all")]])
  }, numeric(1))
)
output_dir <- file.path(root, "results", "gradient_forest", "comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  correlations, file.path(output_dir, "pca_correlations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pdf(file.path(output_dir, "pca_comparison.pdf"), width = 8, height = 5)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
for (period in c("current", "lgm")) {
  for (axis in 1:3) {
    keep <- joined$period == period
    x <- joined[[paste0("PC", axis, "_adaptive")]][keep]
    y <- joined[[paste0("PC", axis, "_all")]][keep]
    plot(
      x, y, pch = 16, cex = 0.2, col = rgb(0, 0, 0, 0.15),
      xlab = paste0("Adaptive PC", axis), ylab = paste0("All-SNP PC", axis),
      main = period
    )
    abline(lm(y ~ x), col = "firebrick", lwd = 2)
  }
}
dev.off()
print(correlations)
