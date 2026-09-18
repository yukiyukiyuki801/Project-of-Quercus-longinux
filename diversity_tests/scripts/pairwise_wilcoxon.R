#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: pairwise_wilcoxon.R INPUT.tsv OUTPUT.tsv")
}

dat <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
groups <- c("kuoi", "lativiolaciifolia", "longinux")
metrics <- c("He", "Ho", "FIS", "pi")
pairs <- combn(groups, 2, simplify = FALSE)

rows <- list()
index <- 1
for (metric in metrics) {
  for (pair in pairs) {
    x <- dat[dat$variety == pair[1], metric]
    y <- dat[dat$variety == pair[2], metric]
    test <- wilcox.test(
      x,
      y,
      alternative = "two.sided",
      paired = FALSE,
      exact = TRUE
    )
    rows[[index]] <- data.frame(
      metric = metric,
      group_1 = pair[1],
      group_2 = pair[2],
      n_1 = length(x),
      n_2 = length(y),
      W = unname(test$statistic),
      p_raw = test$p.value
    )
    index <- index + 1
  }
}

result <- do.call(rbind, rows)
result$p_holm <- ave(
  result$p_raw,
  result$metric,
  FUN = function(p) p.adjust(p, method = "holm")
)
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.table(result, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
