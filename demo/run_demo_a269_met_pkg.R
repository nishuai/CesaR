## Cesar A269 MET-amplification demo (package API)
##
## Same problem as run_demo_a269_met.R, but using the public Cesar
## package API (cesar_train + cesar_detect) instead of the hand-rolled
## anchor fitting.
##
## Truth: A269010002 / 003 / 004 carry MET amplifications;
## the remaining 6 (005 / 006 / 007 / 009 / 010 / 011) are MET-negative
## and serve as the in-cohort PoN.
##
## Note on copy_ratio: cesar_train() chooses anchors purely by
## correlation across the PoN — it does NOT exclude same-gene anchors.
## When the test sample has a whole-gene amplification, some of its own
## MET segments end up acting as anchors for each other; the anchor
## mean rises with the segment depth and r_obs / copy_ratio gets pulled
## toward 1. The qualitative call (conf, sign of copy_ratio - 1) stays
## correct; the magnitude is conservative compared to the hand-rolled
## demo (run_demo_a269_met.R), which excludes same-gene anchors.
##
## Run from the project root with:
##   Rscript demo/run_demo_a269_met_pkg.R

library(Cesar)

setwd('../')
getwd()
source("R/bed_depth_in_pileup.R")

A269_DIR   <- "A269"
BED_FILE   <- "segmented_bed_A269_hg38_annotated.bed"
DEMO_DIR   <- "demo"
RESULT_DIR <- file.path(DEMO_DIR, "results_a269_met_pkg")
FIG_DIR    <- file.path(DEMO_DIR, "figures_a269_met_pkg")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)

TARGET_GENE  <- "MET"
MET_POSITIVE <- c("A269010002", "A269010003", "A269010004")
MET_NEGATIVE <- c("A269010005", "A269010006", "A269010007",
                  "A269010009", "A269010010", "A269010011")

cap <- function(expr, file) {
  con <- file(file, open = "w", encoding = "UTF-8")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) })
  invisible(force(expr))
}

cap({
  cat("================ A269 MET demo (Cesar package) ================\n")
  cat("Cesar version: ", as.character(packageVersion("Cesar")), "\n", sep = "")
  cat("MET-positive (3): ", paste(MET_POSITIVE, collapse = ", "), "\n", sep = "")
  cat("MET-negative (6): ", paste(MET_NEGATIVE, collapse = ", "), "\n", sep = "")

  bed <- read.table(BED_FILE, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE,
                    col.names = c("chr", "start", "end", "gene"))
  freq_files   <- list.files(A269_DIR, pattern = "\\.freq$", full.names = TRUE)
  sample_names <- sub("F01\\.sort\\.bam\\.mpileup\\.freq$", "",
                      basename(freq_files))

  cat(sprintf("\nReading %d samples x %d segments from %s/\n",
              length(freq_files), nrow(bed), A269_DIR))
  mat <- matrix(0, nrow = length(freq_files), ncol = nrow(bed))
  for (i in seq_along(freq_files)) {
    cat("  ", sample_names[i], "\n")
    tab <- read.table(freq_files[i], sep = "\t", header = TRUE)
    mat[i, ] <- bed_depth_in_pileup(tab, bed)
  }
  rownames(mat) <- sample_names

  pon_mat <- mat[MET_NEGATIVE, , drop = FALSE]
  cat(sprintf("\nTraining cesar_train() on %d MET-negative samples ...\n",
              nrow(pon_mat)))
  ## Cesar 1.1.0+: same-gene anchor exclusion ON (MET = 25/1645 segs, amplified
  ## in 3/9 samples — without exclusion the anchor mean tracks the target).
  model <- cesar_train(
    depth_matrix      = pon_mat,
    bed               = bed,
    min_mean_depth    = 50,
    min_anchors       = 3,
    max_anchors       = 25,
    verbose           = TRUE
  )
  saveRDS(model, file.path(RESULT_DIR, "model.rds"))
  print(model)

  target_idx <- which(bed$gene == TARGET_GENE)
  cat(sprintf("\n%s segments in panel: %d\n", TARGET_GENE, length(target_idx)))

  cat("\nRunning cesar_detect() on all 9 samples ...\n")
  results <- lapply(sample_names, function(s)
    cesar_detect(model, mat[s, ], sample_id = s))
  names(results) <- sample_names

  ## Long per-gene table
  gene_long <- do.call(rbind, lapply(results, function(r) {
    g <- r$genes; g$sample <- r$sample; g
  }))
  gene_long$is_target_pos <- gene_long$sample %in% MET_POSITIVE
  write.table(
    gene_long[, c("sample", "gene", "is_target_pos",
                  "n_segments", "copy_ratio", "confidence")],
    file.path(RESULT_DIR, "all_genes_per_sample.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  ## Segment-level rows for MET
  target_seg <- do.call(rbind, lapply(results, function(r) {
    s <- r$segments[r$segments$gene == TARGET_GENE, ]
    s$sample        <- r$sample
    s$is_target_pos <- r$sample %in% MET_POSITIVE
    s$cn_abs        <- 2 * s$copy_ratio
    s$seg_idx       <- seq_len(nrow(s))
    s
  }))
  write.table(
    target_seg[, c("sample", "is_target_pos", "seg_idx",
                   "chr", "start", "end", "gene",
                   "depth", "anchor_ratio", "copy_ratio",
                   "z_score", "p_value", "neglog10_p", "cn_abs")],
    file.path(RESULT_DIR, "MET_segment_level.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  ## Per-sample summary
  summary_tbl <- do.call(rbind, lapply(results, function(r) {
    g <- r$genes[r$genes$gene == TARGET_GENE, ]
    s <- r$segments[r$segments$gene == TARGET_GENE, ]
    cn_abs <- 2 * s$copy_ratio
    data.frame(
      sample        = r$sample,
      is_target_pos = r$sample %in% MET_POSITIVE,
      n_seg         = nrow(s),
      copy_ratio    = round(g$copy_ratio, 3),
      cn_abs_median = round(median(cn_abs, na.rm = TRUE), 3),
      cn_abs_max    = round(max(cn_abs,    na.rm = TRUE), 3),
      cn_abs_min    = round(min(cn_abs,    na.rm = TRUE), 3),
      conf          = round(g$confidence, 2),
      max_neglog10p = round(max(s$neglog10_p, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    )
  }))
  summary_tbl <- summary_tbl[order(!summary_tbl$is_target_pos,
                                    summary_tbl$sample), ]
  write.table(summary_tbl, file.path(RESULT_DIR, "MET_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n================ MET summary ================\n")
  print(summary_tbl, row.names = FALSE)

  ## Plots
  out_bar <- file.path(FIG_DIR, "met_per_sample.png")
  grDevices::png(out_bar, width = 1500, height = 700, res = 130)
  graphics::par(mar = c(7, 5, 4, 1))
  d <- summary_tbl
  cols <- ifelse(d$is_target_pos, "#d6604d", "#92c5de")
  ymax <- max(2 * d$copy_ratio, na.rm = TRUE) * 1.15
  bp <- graphics::barplot(2 * d$copy_ratio, names.arg = d$sample,
                          col = cols, las = 2, ylim = c(0, ymax),
                          ylab = "MET gene-level CN  (= 2 * copy_ratio from cesar_detect)",
                          main = "A269 MET demo (pkg API) - per-sample copy number")
  graphics::abline(h = 2, col = "grey50", lty = 3)
  graphics::text(bp, 2 * d$copy_ratio + ymax * 0.025,
                 labels = sprintf("conf=%.1f", d$conf), cex = 0.78)
  graphics::legend("topright", fill = c("#d6604d", "#92c5de"),
                   legend = c("MET amplified (truth)", "MET negative (PoN)"),
                   bty = "n")
  grDevices::dev.off()

  cn_mat <- do.call(rbind, lapply(sample_names, function(s) {
    rows <- target_seg[target_seg$sample == s, ]
    rows <- rows[order(rows$seg_idx), ]
    rows$cn_abs
  }))
  rownames(cn_mat) <- sample_names
  ord_h <- order(!sample_names %in% MET_POSITIVE, sample_names)
  cn_mat <- cn_mat[ord_h, , drop = FALSE]
  cn_clip <- pmax(pmin(cn_mat, 5), 0)

  out_hm <- file.path(FIG_DIR, "met_segment_heatmap.png")
  grDevices::png(out_hm, width = 1500, height = 700, res = 130)
  graphics::par(mar = c(5, 8, 4, 6))
  pal <- grDevices::colorRampPalette(
    c("#053061", "#2166ac", "white", "#d6604d", "#67001f"))(100)
  breaks <- seq(0, 5, length.out = 101)
  graphics::image(seq_len(ncol(cn_clip)), seq_len(nrow(cn_clip)),
                  t(cn_clip), col = pal, breaks = breaks,
                  axes = FALSE, xlab = "MET segment index", ylab = "",
                  main = "A269 MET demo (pkg API) - per-segment CN (clipped at 5)")
  graphics::axis(1, at = seq(1, ncol(cn_clip), by = 2))
  graphics::axis(2, at = seq_len(nrow(cn_clip)),
                 labels = paste0(rownames(cn_clip),
                                 ifelse(rownames(cn_clip) %in% MET_POSITIVE,
                                        " *", "")),
                 las = 2, cex.axis = 0.85)
  graphics::box()
  pu <- graphics::par("usr")
  bar_x <- pu[2] + (pu[2] - pu[1]) * 0.02
  bar_w <- (pu[2] - pu[1]) * 0.025
  bar_y <- seq(pu[3], pu[4], length.out = 101)
  for (k in seq_len(100))
    graphics::rect(bar_x, bar_y[k], bar_x + bar_w, bar_y[k + 1],
                   col = pal[k], border = NA, xpd = TRUE)
  graphics::text(bar_x + bar_w * 1.5,
                 seq(pu[3], pu[4], length.out = 6),
                 labels = sprintf("%.0f", seq(0, 5, length.out = 6)),
                 xpd = TRUE, cex = 0.8, adj = 0)
  grDevices::dev.off()

  out_sc <- file.path(FIG_DIR, "met_vs_conf.png")
  grDevices::png(out_sc, width = 1100, height = 800, res = 130)
  graphics::par(mar = c(5, 5, 4, 1))
  conf_clip <- pmin(d$conf, 60)
  graphics::plot(conf_clip, d$copy_ratio, pch = 19, cex = 1.6,
                 xlim = c(-2, 68),
                 col = ifelse(d$is_target_pos, "#d6604d", "#92c5de"),
                 xlab = "MET confidence  (mean -log10(p) over MET segs, clipped at 60)",
                 ylab = "MET copy_ratio  (mu_train / r_obs; >1 amp, <1 loss)",
                 main = "A269 MET demo (pkg API) - separation of positives from PoN")
  graphics::abline(h = 1, col = "grey50", lty = 3)
  graphics::abline(v = 2, col = "grey50", lty = 3)
  lab_pos <- ifelse(conf_clip > 50, 2, 4)
  graphics::text(conf_clip, d$copy_ratio, labels = d$sample,
                 pos = lab_pos, offset = 0.5, cex = 0.75)
  graphics::legend("topleft", pch = 19, col = c("#d6604d", "#92c5de"),
                   legend = c("MET amplified", "MET negative"), bty = "n")
  grDevices::dev.off()

  cat(sprintf("\nDemo done. Tables in %s/\n            Figures in %s/\n",
              normalizePath(RESULT_DIR, winslash = "/"),
              normalizePath(FIG_DIR,    winslash = "/")))
  invisible(NULL)
}, file.path(RESULT_DIR, "log.txt"))
