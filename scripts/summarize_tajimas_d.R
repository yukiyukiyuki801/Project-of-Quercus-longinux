#!/usr/bin/env Rscript

required <- c("ggplot2", "maps")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "."
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
input_dir <- file.path(root, "results", "tajimas_d", "windows")
output_dir <- file.path(root, "results", "tajimas_d", "summary")
metadata_path <- file.path(root, "data", "sample_metadata.tsv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  input_dir, pattern = "_SNP_ONLY_EXPLORATORY\\.Tajima\\.D$", full.names = TRUE
)
if (!length(files)) stop("No Tajima's D window files found in ", input_dir)

read_one <- function(path) {
  values <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  names(values) <- sub("^TajimaD$", "Tajima_D", names(values))
  if (!"Tajima_D" %in% names(values)) stop("No TajimaD column in ", path)
  stem <- sub("_[0-9]+bp_SNP_ONLY_EXPLORATORY\\.Tajima\\.D$", "", basename(path))
  pieces <- strsplit(stem, "_", fixed = TRUE)[[1]]
  values$group_type <- pieces[1]
  values$group <- paste(pieces[-1], collapse = "_")
  values
}

windows <- do.call(rbind, lapply(files, read_one))
windows <- windows[is.finite(windows$Tajima_D), ]
summarize_group <- function(values) {
  data.frame(
    n_windows = length(values), mean = mean(values), median = median(values),
    sd = sd(values), min = min(values), max = max(values)
  )
}
summary_rows <- do.call(
  rbind,
  lapply(split(windows, list(windows$group_type, windows$group), drop = TRUE), function(x) {
    cbind(group_type = x$group_type[1], group = x$group[1], summarize_group(x$Tajima_D))
  })
)
write.table(
  summary_rows, file.path(output_dir, "tajimas_d_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

colors <- c(kuoi = "#4C956C", longinux = "#4C78A8", lativiolaciifolia = "#9C6ADE")
variety_windows <- windows[windows$group_type == "variety", ]
variety_windows$group <- factor(variety_windows$group, levels = names(colors))
p_variety <- ggplot2::ggplot(
  variety_windows, ggplot2::aes(x = group, y = Tajima_D, fill = group)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.2) +
  ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
  ggplot2::labs(
    x = "Variety", y = "Tajima's D",
    subtitle = "Exploratory: filtered LD-pruned variant-only RAD VCF"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none")
ggplot2::ggsave(
  file.path(output_dir, "tajimas_d_by_variety.pdf"), p_variety,
  width = 7, height = 4
)

population_windows <- windows[windows$group_type == "population", ]
p_population <- ggplot2::ggplot(
  population_windows,
  ggplot2::aes(x = reorder(group, Tajima_D, median), y = Tajima_D)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.15) +
  ggplot2::labs(
    x = "Population", y = "Tajima's D",
    subtitle = "Exploratory: filtered LD-pruned variant-only RAD VCF"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
ggplot2::ggsave(
  file.path(output_dir, "tajimas_d_by_population.pdf"), p_population,
  width = 10, height = 5
)

metadata <- read.table(metadata_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
population_metadata <- unique(metadata[, c("population", "variety", "longitude", "latitude")])
population_summary <- summary_rows[summary_rows$group_type == "population", ]
population_summary$mean <- as.numeric(population_summary$mean)
map_values <- merge(
  population_summary, population_metadata,
  by.x = "group", by.y = "population"
)
taiwan <- ggplot2::map_data("world", region = "Taiwan")
p_map <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = taiwan, ggplot2::aes(long, lat, group = group),
    fill = "grey95", color = "grey40", linewidth = 0.25
  ) +
  ggplot2::geom_point(
    data = map_values,
    ggplot2::aes(longitude, latitude, color = mean), size = 3
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "Mean D"
  ) +
  ggplot2::coord_quickmap(xlim = c(120.6, 122.0), ylim = c(21.8, 25.4)) +
  ggplot2::labs(
    x = "Longitude", y = "Latitude",
    subtitle = "Exploratory: filtered LD-pruned variant-only RAD VCF"
  ) +
  ggplot2::theme_classic()
ggplot2::ggsave(
  file.path(output_dir, "tajimas_d_population_map.pdf"), p_map,
  width = 5, height = 7
)
message("Wrote Tajima's D summaries and figures to ", output_dir)
