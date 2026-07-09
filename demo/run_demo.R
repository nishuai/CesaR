## Cesar demo — runs end-to-end and writes everything needed for the .md
##
## Outputs:
##   demo/figures/*.png    — per-sample CNV plots
##   demo/results/*.txt    — captured stdout of each step
##   demo/results/per_sample.tsv  — flat results table
##
## Run from the project root with:
##   Rscript demo/run_demo.R

suppressPackageStartupMessages({ library(Cesar) })

DEMO_DIR    <- "demo"
FIG_DIR     <- file.path(DEMO_DIR, "figures")
RESULT_DIR  <- file.path(DEMO_DIR, "results")
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## 1. Locate inputs that ship with the installed package
## ---------------------------------------------------------------------------
bed <- system.file("extdata/panel_30k.bed", package = "Cesar")
pon <- system.file("extdata/pon",           package = "Cesar")
pos <- system.file("extdata/test_samples/positive_metERBB2.depth.gz",
                   package = "Cesar")
neg <- system.file("extdata/test_samples/negative_normal.depth.gz",
                   package = "Cesar")

cat("Cesar demo — input fixtures shipped with the package\n")
cat("---------------------------------------------------\n")
cat(sprintf("  panel BED : %s\n", bed))
cat(sprintf("  PoN dir   : %s  (%d files)\n", pon, length(list.files(pon))))
cat(sprintf("  positive  : %s\n", basename(pos)))
cat(sprintf("  negative  : %s\n", basename(neg)))

## ---------------------------------------------------------------------------
## 2. Train a model from the bundled panel of normals
## ---------------------------------------------------------------------------
cap <- function(expr, file) {
  con <- file(file, open = "w", encoding = "UTF-8")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) })
  invisible(force(expr))
}

cap({
  cat("\n==================== TRAIN ====================\n")
  t0 <- Sys.time()
  model <<- cesar_train(pon_dir = pon, bed_file = bed, verbose = TRUE)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("\nTraining time: %.2f s\n\n", dt))
  summary(model)            # prints inside, no extra wrap
  invisible(NULL)
}, file.path(RESULT_DIR, "01_train.txt"))

## ---------------------------------------------------------------------------
## 3. Detect CNVs in the positive and the held-out negative
## ---------------------------------------------------------------------------
detect_one <- function(sample_path, label) {
  cap({
    cat(sprintf("\n==================== DETECT — %s ====================\n",
                label))
    res <- cesar_detect(model, sample_path)
    summary(res)
    cat("\nTop 5 segments by -log10(p):\n")
    s <- res$segments
    print(head(s[order(-s$neglog10_p), ], 5L), row.names = FALSE)
    invisible(NULL)
  }, file.path(RESULT_DIR, sprintf("02_detect_%s.txt", label)))
  res <- cesar_detect(model, sample_path)
  res
}

res_pos <- detect_one(pos, "positive")
res_neg <- detect_one(neg, "negative")

## ---------------------------------------------------------------------------
## 4. Render per-sample CNV plot (PNG, embeddable in the .md)
## ---------------------------------------------------------------------------
save_plot <- function(res, name) {
  out <- file.path(FIG_DIR, paste0(name, ".png"))
  grDevices::png(out, width = 1400, height = 600, res = 130)
  graphics::par(mar = c(4, 4, 3, 1))
  plot(res)
  grDevices::dev.off()
  out
}
fig_pos <- save_plot(res_pos, "positive_metERBB2")
fig_neg <- save_plot(res_neg, "negative_normal")

## ---------------------------------------------------------------------------
## 5. Flat per-sample comparison table
## ---------------------------------------------------------------------------
gene_row <- function(res, sample_label) {
  pick <- function(g) {
    h <- res$genes[res$genes$gene == g, ]
    if (nrow(h) == 0L) c(NA_real_, NA_real_)
    else c(h$copy_ratio, h$confidence)
  }
  m <- pick("MET"); e <- pick("ERBB2"); g <- pick("EGFR")
  data.frame(
    sample     = sample_label,
    MET_copy   = round(m[1], 3),  MET_conf   = round(m[2], 2),
    ERBB2_copy = round(e[1], 3),  ERBB2_conf = round(e[2], 2),
    EGFR_copy  = round(g[1], 3),  EGFR_conf  = round(g[2], 2),
    stringsAsFactors = FALSE
  )
}
flat <- rbind(gene_row(res_pos, "positive (MET~2.75x, ERBB2~4.75x amp)"),
              gene_row(res_neg, "negative (held-out normal)"))
write.table(flat, file.path(RESULT_DIR, "per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDemo done. Files written to ", normalizePath(DEMO_DIR), "\n", sep = "")
print(flat, row.names = FALSE)
