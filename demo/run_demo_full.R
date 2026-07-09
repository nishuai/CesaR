## Cesar full-cohort demo — runs train + detect across every sample
## in the bundled cesar_demo dataset (10 normals + 18 spike-ins) and
## writes everything needed for the .md.
##
## Run from the project root with:
##   Rscript demo/run_demo_full.R
##
## Outputs:
##   demo/figures/dose_response.png   — per-condition copy_ratio boxplot
##   demo/figures/heatmap.png         — per-segment z-score heatmap
##   demo/figures/sample_*.png        — one CNV plot per representative sample
##   demo/results/full_per_sample.tsv — flat table, all 28 samples
##   demo/results/full_train.txt      — captured training stdout

suppressPackageStartupMessages({ library(Cesar) })

DEMO_DIR    <- "demo"
FIG_DIR     <- file.path(DEMO_DIR, "figures")
RESULT_DIR  <- file.path(DEMO_DIR, "results")
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

cap <- function(expr, file) {
  con <- file(file, open = "w", encoding = "UTF-8")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) })
  invisible(force(expr))
}

## ---------------------------------------------------------------------------
## 1. Load the bundled dataset
## ---------------------------------------------------------------------------
data("cesar_demo", package = "Cesar")
meta <- cesar_demo$sample_meta
cat(sprintf("Loaded cesar_demo: %d samples x %d segments  (built %s)\n",
            nrow(cesar_demo$depth_matrix),
            ncol(cesar_demo$depth_matrix),
            cesar_demo$built_on))
print(table(condition = meta$condition))

## ---------------------------------------------------------------------------
## 2. Hold out 3 normals as negative controls; train on the rest
## ---------------------------------------------------------------------------
set.seed(11)
normal_idx <- which(meta$condition == "normal")
holdout    <- sample(normal_idx, 3L)
pon_idx    <- setdiff(normal_idx, holdout)

cap({
  cat("\n==================== TRAIN ====================\n")
  cat(sprintf("PoN samples : %d   held-out negatives: %d\n",
              length(pon_idx), length(holdout)))
  t0 <- Sys.time()
  model <<- cesar_train(
    depth_matrix = cesar_demo$depth_matrix[pon_idx, ],
    bed          = cesar_demo$bed,
    verbose      = TRUE
  )
  cat(sprintf("\nTraining time: %.2f s\n\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  summary(model)
  invisible(NULL)
}, file.path(RESULT_DIR, "full_train.txt"))

## ---------------------------------------------------------------------------
## 3. Detect on EVERY sample, including the held-out normals
## ---------------------------------------------------------------------------
detect_idx <- c(holdout, which(meta$condition != "normal"))
detect_idx <- sort(detect_idx)

results <- lapply(detect_idx, function(i) {
  res <- cesar_detect(
    model,
    cesar_demo$depth_matrix[i, ],
    sample_id = sprintf("%s [%s]", meta$sample_id[i], meta$condition[i])
  )
  list(idx = i, res = res)
})

## Flatten to a per-sample summary table.
gene_pick <- function(res, g) {
  h <- res$genes[res$genes$gene == g, ]
  if (nrow(h) == 0L) return(c(NA_real_, NA_real_))
  c(h$copy_ratio, h$confidence)
}
flat <- do.call(rbind, lapply(results, function(r) {
  i <- r$idx
  m <- gene_pick(r$res, "MET")
  e <- gene_pick(r$res, "ERBB2")
  g <- gene_pick(r$res, "EGFR")
  cnd <- meta$condition[i]
  is_holdout <- (cnd == "normal" && i %in% holdout)
  data.frame(
    sample      = meta$sample_id[i],
    condition   = if (is_holdout) "normal_holdout" else cnd,
    MET_copy    = round(m[1], 3), MET_conf    = round(m[2], 2),
    ERBB2_copy  = round(e[1], 3), ERBB2_conf  = round(e[2], 2),
    EGFR_copy   = round(g[1], 3), EGFR_conf   = round(g[2], 2),
    stringsAsFactors = FALSE
  )
}))
flat <- flat[order(flat$condition, flat$sample), ]
write.table(flat, file.path(RESULT_DIR, "full_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Per-condition median + IQR for the markdown summary table.
agg_fun <- function(v) sprintf("%.2f [%.2f, %.2f]",
                               median(v, na.rm = TRUE),
                               quantile(v, 0.25, na.rm = TRUE),
                               quantile(v, 0.75, na.rm = TRUE))
agg <- do.call(rbind, by(flat, flat$condition, function(d)
  data.frame(
    condition  = d$condition[1],
    n          = nrow(d),
    MET_copy   = agg_fun(d$MET_copy),
    MET_conf   = agg_fun(d$MET_conf),
    ERBB2_copy = agg_fun(d$ERBB2_copy),
    ERBB2_conf = agg_fun(d$ERBB2_conf),
    stringsAsFactors = FALSE
  )))
agg <- agg[order(c(normal_holdout = 1, met2.1875ERBB2.625 = 2,
                   met2.375ERBB3.25 = 3, met2.75ERBB4.75 = 4)[agg$condition]), ]

cap({
  cat("\n==================== PER-SAMPLE TABLE ====================\n")
  print(flat, row.names = FALSE)
  cat("\n==================== PER-CONDITION SUMMARY ====================\n")
  print(agg, row.names = FALSE)
  invisible(NULL)
}, file.path(RESULT_DIR, "full_summary.txt"))

## ---------------------------------------------------------------------------
## 4. Plot 1: dose-response boxplot (copy_ratio vs condition, MET & ERBB2)
## ---------------------------------------------------------------------------
plot_dose_response <- function() {
  cond_order <- c("normal_holdout", "met2.1875ERBB2.625",
                  "met2.375ERBB3.25", "met2.75ERBB4.75")
  cond_labels <- c("normal\n(n=3)", "MET 2.19x\nERBB2 0.63x",
                   "MET 2.38x\nERBB2 3.25x", "MET 2.75x\nERBB2 4.75x")
  flat$condition <- factor(flat$condition, levels = cond_order)
  graphics::par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  cols <- c("grey80", "#fde0a6", "#f0a668", "#d35400")

  graphics::boxplot(MET_copy ~ condition, data = flat,
                    names = cond_labels, las = 1, col = cols,
                    main = "MET — copy_ratio per sample",
                    ylab = "copy_ratio  (mu_train / r_obs)",
                    xlab = "")
  graphics::abline(h = 1, col = "darkgrey", lty = 3)
  graphics::points(jitter(as.integer(flat$condition), 0.5),
                   flat$MET_copy, pch = 19, col = "black", cex = 0.7)

  graphics::boxplot(ERBB2_copy ~ condition, data = flat,
                    names = cond_labels, las = 1, col = cols,
                    main = "ERBB2 — copy_ratio per sample",
                    ylab = "copy_ratio  (mu_train / r_obs)",
                    xlab = "")
  graphics::abline(h = 1, col = "darkgrey", lty = 3)
  graphics::points(jitter(as.integer(flat$condition), 0.5),
                   flat$ERBB2_copy, pch = 19, col = "black", cex = 0.7)
}
out_dr <- file.path(FIG_DIR, "dose_response.png")
grDevices::png(out_dr, width = 1600, height = 700, res = 130)
plot_dose_response()
grDevices::dev.off()

## ---------------------------------------------------------------------------
## 5. Plot 2: per-segment z-score heatmap, samples sorted by condition
## ---------------------------------------------------------------------------
zmat <- do.call(rbind, lapply(results, function(r) r$res$segments$z_score))
rownames(zmat) <- vapply(results, function(r)
  sprintf("%s [%s]",
          meta$sample_id[r$idx],
          if (meta$condition[r$idx] == "normal" && r$idx %in% holdout)
            "normal_holdout" else meta$condition[r$idx]),
  character(1))

## Order samples by condition so the heatmap reads top-to-bottom.
cond_label <- ifelse(grepl("normal_holdout", rownames(zmat)), 1L,
              ifelse(grepl("met2.1875ERBB2.625", rownames(zmat)), 2L,
              ifelse(grepl("met2.375ERBB3.25",   rownames(zmat)), 3L,
              ifelse(grepl("met2.75ERBB4.75",    rownames(zmat)), 4L, 5L))))
zmat <- zmat[order(cond_label), , drop = FALSE]

## Squash z-scores to ±10 for legibility.
zclip <- pmax(pmin(zmat, 10), -10)

out_hm <- file.path(FIG_DIR, "heatmap.png")
grDevices::png(out_hm, width = 1800, height = 900, res = 130)
graphics::par(mar = c(5, 12, 4, 4))
breaks <- seq(-10, 10, length.out = 101)
pal <- grDevices::colorRampPalette(c("#053061", "#2166ac", "#92c5de",
                                     "white",
                                     "#f4a582", "#d6604d", "#67001f"))(100)
graphics::image(seq_len(ncol(zclip)), seq_len(nrow(zclip)),
                t(zclip), col = pal, breaks = breaks,
                xlab = "Segment index", ylab = "",
                axes = FALSE,
                main = "Per-segment z-score across the dilution series")
graphics::axis(1, at = seq(0, ncol(zclip), 25))
graphics::axis(2, at = seq_len(nrow(zclip)),
               labels = rownames(zclip), las = 2, cex.axis = 0.55)
## Mark the gene-of-interest segment indices.
seg_genes <- cesar_demo$bed$gene
mark <- function(g, col) {
  idx <- which(seg_genes == g)
  if (length(idx) > 0L) {
    graphics::abline(v = c(min(idx) - 0.5, max(idx) + 0.5),
                     col = col, lty = 2)
    graphics::mtext(g, side = 3, at = mean(idx), col = col, cex = 0.9)
  }
}
mark("EGFR",  "#1b7837")
mark("MET",   "#762a83")
mark("ERBB2", "#b2182b")
grDevices::dev.off()

## ---------------------------------------------------------------------------
## 6. One representative CNV plot per condition (Cesar's built-in plot)
## ---------------------------------------------------------------------------
rep_per_condition <- function() {
  for (cnd in c("normal_holdout", "met2.1875ERBB2.625",
                "met2.375ERBB3.25", "met2.75ERBB4.75")) {
    rows <- which(grepl(cnd, vapply(results,
                          function(r) sprintf("%s",
                            if (meta$condition[r$idx] == "normal" &&
                                r$idx %in% holdout) "normal_holdout"
                            else meta$condition[r$idx]),
                          character(1)), fixed = TRUE))
    if (length(rows) == 0L) next
    rep_res <- results[[rows[1]]]$res
    out <- file.path(FIG_DIR, sprintf("sample_%s.png", cnd))
    grDevices::png(out, width = 1400, height = 500, res = 130)
    graphics::par(mar = c(4, 4, 3, 1))
    plot(rep_res)
    grDevices::dev.off()
  }
}
rep_per_condition()

cat("\nFull demo done. Files written to ",
    normalizePath(DEMO_DIR), "\n", sep = "")
print(agg, row.names = FALSE)
