required_gf_packages <- c("gradientForest", "raster", "rgdal")
gf_predictors <- c("BIO01", "BIO07", "BIO12", "BIO17", "BIO18")
gf_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

check_gf_packages <- function() {
  missing <- required_gf_packages[
    !vapply(required_gf_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))
}

read_gf_inputs <- function(project_root, snp_set) {
  genotype_path <- file.path(project_root, "180_ql.ld_prune.lfmm.csv")
  environment_path <- file.path(project_root, "ql_180_env_variables.csv")
  candidate_path <- file.path(project_root, "ql_rda_fst_overlap_outliers.csv")
  for (path in c(genotype_path, environment_path)) {
    if (!file.exists(path)) stop("Missing required input: ", path)
  }

  genotype <- read.csv(
    genotype_path, row.names = 1, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  environment <- read.csv(
    environment_path, row.names = 1, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(rownames(genotype))) stop("Duplicate sample IDs in genotype matrix")
  if (anyDuplicated(rownames(environment))) stop("Duplicate sample IDs in environment table")
  missing_samples <- setdiff(rownames(genotype), rownames(environment))
  if (length(missing_samples)) {
    stop("Environment table is missing samples: ", paste(missing_samples, collapse = ", "))
  }
  missing_predictors <- setdiff(gf_predictors, colnames(environment))
  if (length(missing_predictors)) {
    stop("Environment table is missing predictors: ", paste(missing_predictors, collapse = ", "))
  }
  environment <- environment[rownames(genotype), gf_predictors, drop = FALSE]
  if (anyNA(environment)) stop("Missing environmental values after matching sample IDs")

  if (snp_set == "adaptive") {
    if (!file.exists(candidate_path)) stop("Missing candidate-locus table: ", candidate_path)
    candidates <- read.csv(candidate_path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"snp" %in% names(candidates)) stop("Candidate table needs a column named 'snp'")
    retained <- intersect(candidates$snp, colnames(genotype))
    removed <- setdiff(candidates$snp, colnames(genotype))
    if (!length(retained)) stop("No candidate loci overlap the LD-pruned matrix")
    message(
      "Adaptive-overlap set: ", length(retained), " of ", nrow(candidates),
      " candidate loci retained; ", length(removed), " absent after LD pruning."
    )
    genotype <- genotype[, retained, drop = FALSE]
  }
  if (anyNA(genotype)) stop("Genotype matrix contains missing values")
  list(genotype = genotype, environment = environment)
}

read_predictor_stack <- function(raster_dir) {
  paths <- file.path(raster_dir, paste0(gf_predictors, ".asc"))
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing raster(s): ", paste(missing, collapse = ", "))
  layers <- lapply(paths, raster::raster)
  stk <- raster::stack(layers)
  names(stk) <- gf_predictors
  raster::crs(stk) <- gf_crs
  stk
}

predict_habitat <- function(model, raster_dir, habitat_path, period) {
  if (!file.exists(habitat_path)) stop("Missing habitat shapefile: ", habitat_path)
  habitat <- rgdal::readOGR(habitat_path, verbose = FALSE)
  raster::crs(habitat) <- gf_crs
  predictors <- read_predictor_stack(raster_dir)
  predictors <- raster::crop(predictors, raster::extent(habitat))
  predictors <- raster::mask(predictors, habitat)
  points <- as.data.frame(raster::rasterToPoints(predictors))
  names(points)[1:2] <- c("longitude", "latitude")
  transformed <- predict(model, points[, gf_predictors, drop = FALSE])
  result <- cbind(points[, c("longitude", "latitude")], transformed)
  result$period <- period
  result
}

scale_zero_one <- function(x) {
  limits <- range(x, na.rm = TRUE)
  if (!all(is.finite(limits)) || diff(limits) == 0) return(rep(0.5, length(x)))
  (x - limits[1]) / diff(limits)
}

run_gf_workflow <- function(project_root, snp_set = c("adaptive", "all"), output_dir) {
  snp_set <- match.arg(snp_set)
  check_gf_packages()
  inputs <- read_gf_inputs(project_root, snp_set)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  tree_count <- as.integer(Sys.getenv("GF_NTREE", "500"))
  if (is.na(tree_count) || tree_count < 1) stop("GF_NTREE must be a positive integer")
  max_level <- log2(0.368 * nrow(inputs$environment) / 2)
  set.seed(1)
  model <- gradientForest::gradientForest(
    cbind(inputs$environment, inputs$genotype),
    predictor.vars = gf_predictors,
    response.vars = colnames(inputs$genotype),
    ntree = tree_count,
    maxLevel = max_level,
    trace = TRUE,
    corr.threshold = 0.50
  )
  saveRDS(model, file.path(output_dir, "gradient_forest_model.rds"))

  importance <- data.frame(
    variable = names(model$overall.imp),
    importance = as.numeric(model$overall.imp)
  )
  write.table(
    importance, file.path(output_dir, "variable_importance.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  current <- predict_habitat(
    model, file.path(project_root, "current_bio19"),
    file.path(project_root, "mask_layer", "current_habitat.shp"), "current"
  )
  lgm <- predict_habitat(
    model, file.path(project_root, "lgm_bio19"),
    file.path(project_root, "mask_layer", "lgm_habitat.shp"), "lgm"
  )
  predictions <- rbind(current, lgm)
  write.table(
    predictions, file.path(output_dir, "gf_predictions.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  pca <- prcomp(predictions[, gf_predictors], center = TRUE, scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:3, drop = FALSE])
  scores$longitude <- predictions$longitude
  scores$latitude <- predictions$latitude
  scores$period <- predictions$period
  write.table(
    scores, file.path(output_dir, "gf_prediction_pca_scores.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  colors <- rgb(
    scale_zero_one(scores$PC1), scale_zero_one(scores$PC2), scale_zero_one(scores$PC3)
  )
  pdf(file.path(output_dir, "gf_prediction_map.pdf"), width = 7, height = 4)
  par(mfrow = c(1, 2), mar = c(3, 3, 2, 1))
  for (period in c("current", "lgm")) {
    keep <- scores$period == period
    plot(
      scores$longitude[keep], scores$latitude[keep], pch = ".", cex = 0.1,
      asp = 1, col = colors[keep], xlab = "Longitude", ylab = "Latitude",
      main = paste(if (snp_set == "adaptive") "Adaptive overlap" else "All SNPs", period)
    )
  }
  dev.off()
  message("Wrote Gradient Forest outputs to ", output_dir)
  invisible(model)
}
