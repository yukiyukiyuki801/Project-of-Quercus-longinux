args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
if (length(script_path) > 0) {
  setwd(dirname(normalizePath(script_path)))
}

old_r_lib <- "/vast/palmer/home.mccleary/ps2327/R/4.2"
if (dir.exists(old_r_lib)) {
  .libPaths(c(old_r_lib, .libPaths()))
}

library(raster)
library(rgdal)
library(vegan)
library(gradientForest)

crs_ll <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
vars <- c("BIO01", "BIO07", "BIO12", "BIO17", "BIO18")

toStk3 <- function(x, varList, rType = "asc", vConvert = FALSE) {
  for (var in varList) {
    if (rType == "asc") {
      inF <- file.path(x, paste0(var, ".asc"))
    } else {
      stop("Only asc rasters are supported by this rerun script")
    }
    r <- raster(inF)
    if (vConvert && names(r) %in% c("MAT", "MWMT", "MCMT", "TD", "EMT", "EXT", "AHM", "SHM", "MAR")) {
      r <- r / 10
    }
    if (var == varList[1]) {
      stk <- r
    } else {
      stk <- stack(stk, r)
    }
  }
  stk
}

fit_gf <- function() {
  gen <- read.csv("180_ql.ld_prune.lfmm.csv", header = TRUE, row.names = 1)
  outlier <- read.csv("ql_rda_fst_overlap_outliers.csv", header = TRUE, row.names = 1)
  gen <- gen[, colnames(gen) %in% rownames(outlier)]

  clim.point <- read.csv("ql_180_env_variables.csv", header = TRUE, row.names = 1)
  clim.point <- clim.point[101:119]
  clim.point <- clim.point[, vars]

  env.gf <- clim.point
  maxLevel <- log2(0.368 * nrow(env.gf) / 2)
  set.seed(1)
  gradientForest(
    cbind(env.gf, gen),
    predictor.vars = colnames(env.gf),
    response.vars = colnames(gen),
    ntree = 500,
    maxLevel = maxLevel,
    trace = TRUE,
    corr.threshold = 0.50
  )
}

make_prediction_df <- function(raster_dir, habitat_shp) {
  range <- readOGR(habitat_shp, verbose = FALSE)
  crs(range) <- crs_ll

  stk <- toStk3(raster_dir, vars, rType = "asc", vConvert = TRUE)
  crs(stk) <- crs_ll
  stk <- crop(stk, extent(range))
  stk <- mask(stk, range)
  crs(stk) <- crs_ll

  pts <- as.data.frame(rasterToPoints(stk))
  names(pts)[1:2] <- c("Longitude", "Latitude")
  pts
}

gf <- fit_gf()

df_current <- make_prediction_df("current_bio19", "mask_layer/current_habitat.shp")
df_lgm <- make_prediction_df("lgm_bio19", "mask_layer/lgm_habitat.shp")

df_pca_current <- cbind(
  df_current[, c("Longitude", "Latitude")],
  predict(gf, df_current[, vars])
)
df_pca_current$time <- "current"

df_pca_lgm <- cbind(
  df_lgm[, c("Longitude", "Latitude")],
  predict(gf, df_lgm[, vars])
)
df_pca_lgm$time <- "lgm"

df_pca <- rbind(df_pca_current, df_pca_lgm)

PC <- prcomp(df_pca[, vars])
pcx <- PC$x

pc1 <- pcx[, 1]
pc2 <- pcx[, 2]
pc3 <- pcx[, 3]

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

outlineABBC <- readOGR("mask_layer/current_background.shp", verbose = FALSE)
crs(outlineABBC) <- crs_ll
outlineABBC_lgm <- readOGR("mask_layer/lgm_background.shp", verbose = FALSE)
crs(outlineABBC_lgm) <- crs_ll

n_current <- nrow(df_pca_current)

cols <- rgb(scale01(pc1), scale01(pc2), scale01(pc3))
cols_current <- cols[seq_len(n_current)]
cols_lgm <- cols[(n_current + 1):length(cols)]

vec <- vars
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec
arrow1 <- PC$rotation[vind, 1]
arrow2 <- PC$rotation[vind, 2]
arrow_scale <- 2

pdf("gf_pca_combine.pdf", height = 4, width = 5)
par(mfrow = c(1, 1))
plot(pcx[, 1:2], pch = ".", cex = 4, col = cols, asp = 1)
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
jit <- 0.03
text(
  arrow1 / arrow_scale + jit * sign(arrow1),
  arrow2 / arrow_scale + jit * sign(arrow2),
  labels = vec,
  cex = 0.7
)
dev.off()

pdf("gf_current.pdf", height = 5, width = 7)
par(mfrow = c(1, 2))
plot(df_pca_current[, c("Longitude", "Latitude")], pch = ".", cex = 0.1, asp = 1, col = cols_current)
plot(outlineABBC, add = TRUE, lwd = 1)
plot(df_pca_lgm[, c("Longitude", "Latitude")], pch = ".", cex = 0.1, asp = 1, col = cols_lgm)
plot(outlineABBC_lgm, add = TRUE, lwd = 1)
dev.off()

cat("Wrote gf_pca_combine.pdf and gf_current.pdf using R=PC1, G=PC2, B=PC3\n")
