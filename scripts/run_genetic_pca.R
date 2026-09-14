#!/usr/bin/env Rscript

required <- c("ggplot2", "vegan")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
genotype_path <- file.path(root, "180_ql.ld_prune.lfmm.csv")
metadata_path <- file.path(root, "data", "sample_metadata.tsv")
output_dir <- file.path(root, "results", "genetic_pca")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genotype <- read.csv(
  genotype_path, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
metadata <- read.table(
  metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  check.names = FALSE
)
if (anyDuplicated(rownames(genotype))) stop("Duplicate IDs in genotype matrix")
if (anyDuplicated(metadata$sample_id)) stop("Duplicate IDs in sample metadata")
if (!setequal(rownames(genotype), metadata$sample_id)) {
  stop("Genotype and metadata sample IDs do not match")
}
if (anyNA(genotype)) stop("Genotype matrix contains missing values")
metadata <- metadata[match(rownames(genotype), metadata$sample_id), ]

pca <- vegan::rda(genotype, scale = FALSE)
eigenvalues <- vegan::eigenvals(pca)
variance <- 100 * eigenvalues / sum(eigenvalues)
n_axes <- min(10, length(eigenvalues))
scores <- as.data.frame(
  vegan::scores(pca, display = "sites", choices = seq_len(n_axes))
)
scores$sample_id <- rownames(scores)
scores <- merge(scores, metadata, by = "sample_id", sort = FALSE)
scores <- scores[match(rownames(genotype), scores$sample_id), ]

write.table(
  scores, file.path(output_dir, "genetic_pca_scores.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(axis = seq_along(eigenvalues), eigenvalue = eigenvalues,
             variance_percent = variance),
  file.path(output_dir, "genetic_pca_eigenvalues.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

variety_colors <- c(
  kuoi = "#4C956C", longinux = "#4C78A8", lativiolaciifolia = "#9C6ADE"
)
plot_pca <- ggplot2::ggplot(
  scores, ggplot2::aes(x = PC1, y = PC2, color = variety)
) +
  ggplot2::geom_point(size = 2, alpha = 0.8) +
  ggplot2::scale_color_manual(values = variety_colors) +
  ggplot2::labs(
    x = sprintf("PC1 (%.2f%%)", variance[1]),
    y = sprintf("PC2 (%.2f%%)", variance[2]),
    color = "Variety"
  ) +
  ggplot2::theme_classic()
ggplot2::ggsave(
  file.path(output_dir, "genetic_pca.pdf"), plot_pca,
  width = 6, height = 5
)
message("Wrote genetic PCA outputs to ", output_dir)
