## CesaR — full cohort demo on the bundled cesar_demo dataset
##
## Trains on the 10 normals and detects MET / ERBB2 CNVs across all 28
## samples of the spike-in dilution series, then prints the per-condition
## dose-response. Uses ONLY the in-package dataset — no external paths.
##
## Run from R:
##   source(system.file("scripts/run_demo_full.R", package = "CesaR"))
## or from the shell:
##   Rscript -e "source(system.file('scripts/run_demo_full.R', package='CesaR'))"

suppressPackageStartupMessages(library(CesaR))

cat("================ CesaR full-cohort demo ================\n")
cat(sprintf("CesaR version: %s\n\n", packageVersion("CesaR")))

## ---------------------------------------------------------------------------
## 1. Load the bundled dataset
## ---------------------------------------------------------------------------
data("cesar_demo", package = "CesaR")
meta <- cesar_demo$sample_meta
cat(sprintf("Loaded cesar_demo: %d samples x %d segments\n",
            nrow(cesar_demo$depth_matrix), ncol(cesar_demo$depth_matrix)))
cat("Samples per condition:\n")
print(table(meta$condition))

## ---------------------------------------------------------------------------
## 2. Train on the 10 normals (same-gene anchor exclusion ON)
## ---------------------------------------------------------------------------
normals <- meta$condition == "normal"
cat(sprintf("\n[train] %d normal samples -> cesar_train()\n", sum(normals)))
model <- cesar_train(
  depth_matrix      = cesar_demo$depth_matrix[normals, ],
  bed               = cesar_demo$bed,
  min_mean_depth    = 50,
  min_anchors       = 3,
  max_anchors       = 25,
  exclude_same_gene = TRUE,
  verbose           = FALSE
)
print(model)

## ---------------------------------------------------------------------------
## 3. Detect on every sample, pull MET + ERBB2 gene-level calls
## ---------------------------------------------------------------------------
cat("\n[detect] scoring all 28 samples ...\n")

## Pull (copy_ratio, confidence) for one gene from a cesar_result.
## Returns a named list so callers can use $copy / $conf instead of [1]/[2].
get_gene <- function(res, g) {
  hit <- res$genes[res$genes$gene == g, ]
  if (nrow(hit) == 0L) {
    list(copy = NA_real_, conf = NA_real_)
  } else {
    list(copy = hit$copy_ratio, conf = hit$confidence)
  }
}

n_samples <- nrow(cesar_demo$depth_matrix)
flat <- data.frame(
  sample     = meta$sample_id,
  condition  = meta$condition,
  MET_copy   = NA_real_,
  MET_conf   = NA_real_,
  ERBB2_copy = NA_real_,
  ERBB2_conf = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_len(n_samples)) {
  res <- cesar_detect(model, cesar_demo$depth_matrix[i, ],
                      sample_id = meta$sample_id[i])
  met   <- get_gene(res, "MET")
  erbb2 <- get_gene(res, "ERBB2")
  flat$MET_copy[i]   <- round(met$copy,   3)
  flat$MET_conf[i]   <- round(met$conf,   1)
  flat$ERBB2_copy[i] <- round(erbb2$copy, 3)
  flat$ERBB2_conf[i] <- round(erbb2$conf, 1)
}

## ---------------------------------------------------------------------------
## 4. Per-condition dose-response (median copy_ratio across replicates)
## ---------------------------------------------------------------------------
cond_order <- c("normal", "met2.1875ERBB2.625",
                "met2.375ERBB3.25", "met2.75ERBB4.75")

agg <- data.frame(
  condition  = cond_order,
  n          = NA_integer_,
  MET_copy   = NA_real_,
  MET_conf   = NA_real_,
  ERBB2_copy = NA_real_,
  ERBB2_conf = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_along(cond_order)) {
  d <- flat[flat$condition == cond_order[i], ]
  agg$n[i]          <- nrow(d)
  agg$MET_copy[i]   <- round(median(d$MET_copy),   3)
  agg$MET_conf[i]   <- round(median(d$MET_conf),   1)
  agg$ERBB2_copy[i] <- round(median(d$ERBB2_copy), 3)
  agg$ERBB2_conf[i] <- round(median(d$ERBB2_conf), 1)
}

cat("\n================ Per-condition dose-response ================\n")
cat("(median across replicates; truth multipliers in cesar_demo$condition_meta)\n\n")
print(agg, row.names = FALSE)

cat("\nDemo complete. The monotonic rise in MET_copy / ERBB2_copy with\n")
cat("spike-in level — and conf ~ 0 in normals — is the expected result.\n")
