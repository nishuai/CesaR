## Cesar A269 ERBB2-amplification demo (package API)
##
## Same problem as run_demo_a269_erbb2.R, but using the public Cesar
## package API (cesar_train + cesar_detect) instead of the hand-rolled
## anchor fitting.
##
## Truth: A269010005 / 006 / 007 / 009 carry ERBB2 amplifications;
## the remaining 5 (002 / 003 / 004 / 010 / 011) are ERBB2-negative and
## serve as the in-cohort panel-of-normals for the ERBB2 locus.
##
## Note on copy_ratio: cesar_train() chooses anchors purely by
## correlation across the PoN — it does NOT exclude same-gene anchors.
## When the test sample has a whole-gene amplification, some of its own
## ERBB2 segments end up acting as anchors for each other; the anchor
## mean rises with the segment depth and r_obs / copy_ratio gets pulled
## toward 1. The qualitative call (conf, sign of copy_ratio - 1) stays
## correct; the magnitude is conservative.
##
## Run from the project root with:
##   Rscript demo/run_demo_a269_erbb2_pkg.R
##
## Outputs:
##   demo/results_a269_erbb2_pkg/ERBB2_summary.tsv
##   demo/results_a269_erbb2_pkg/ERBB2_segment_level.tsv
##   demo/results_a269_erbb2_pkg/all_genes_per_sample.tsv
##   demo/results_a269_erbb2_pkg/log.txt
##   demo/results_a269_erbb2_pkg/model.rds
##   demo/figures_a269_erbb2_pkg/erbb2_per_sample.png
##   demo/figures_a269_erbb2_pkg/erbb2_segment_heatmap.png
##   demo/figures_a269_erbb2_pkg/erbb2_vs_conf.png

suppressPackageStartupMessages({
  library(Cesar)
})
source("R/bed_depth_in_pileup.R")

## ---------------------------------------------------------------------------
## 0. Paths and ground truth
## ---------------------------------------------------------------------------
A269_DIR   <- "A269"
BED_FILE   <- "segmented_bed_A269_hg38_annotated.bed"
DEMO_DIR   <- "demo"
RESULT_DIR <- file.path(DEMO_DIR, "results_a269_erbb2_pkg")
FIG_DIR    <- file.path(DEMO_DIR, "figures_a269_erbb2_pkg")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)

TARGET_GENE    <- "ERBB2"
ERBB2_POSITIVE <- c("A269010005", "A269010006", "A269010007", "A269010009")
ERBB2_NEGATIVE <- c("A269010002", "A269010003", "A269010004",
                    "A269010010", "A269010011")

cap <- function(expr, file) {
  con <- file(file, open = "w", encoding = "UTF-8")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) })
  invisible(force(expr))
}

cap({
  cat("================ A269 ERBB2 demo (Cesar package) ================\n")
  cat("Cesar version: ", as.character(packageVersion("Cesar")), "\n", sep = "")
  cat("ERBB2-positive (4): ", paste(ERBB2_POSITIVE, collapse = ", "), "\n", sep = "")
  cat("ERBB2-negative (5): ", paste(ERBB2_NEGATIVE, collapse = ", "), "\n", sep = "")

  ## -------------------------------------------------------------------------
  ## 1. Build the depth matrix once, reuse for train + detect
  ## -------------------------------------------------------------------------
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

  ## -------------------------------------------------------------------------
  ## 2. Train via the package API on the 5 ERBB2-negative samples
  ## -------------------------------------------------------------------------
  pon_mat <- mat[ERBB2_NEGATIVE, , drop = FALSE]
  cat(sprintf("\nTraining cesar_train() on %d ERBB2-negative samples ...\n",
              nrow(pon_mat)))

  ## Cesar 1.1.0+: turn on same-gene anchor exclusion. With ERBB2 contributing
  ## 27/1645 segments and being amplified in 4/9 samples, we need to keep the
  ## anchor pool ERBB2-free; otherwise a positive sample's own ERBB2 mates
  ## become its highest-correlation anchors and the ratio is attenuated.
  model <- cesar_train(
    depth_matrix      = pon_mat,
    bed               = bed,
    min_mean_depth    = 50,
    min_anchors       = 3,
    max_anchors       = 25,
    exclude_same_gene = TRUE,
    verbose           = TRUE
  )
  saveRDS(model, file.path(RESULT_DIR, "model.rds"))
  print(model)

  target_idx <- which(bed$gene == TARGET_GENE)
  cat(sprintf("\n%s segments in panel: %d\n", TARGET_GENE, length(target_idx)))

  ## -------------------------------------------------------------------------
  ## 3. Detect via the package API on every sample
  ## -------------------------------------------------------------------------
  cat("\nRunning cesar_detect() on all 9 samples ...\n")
  results <- lapply(sample_names, function(s)
    cesar_detect(model, mat[s, ], sample_id = s))
  names(results) <- sample_names

  ## Per-gene aggregate across all samples (long table).
  gene_long <- do.call(rbind, lapply(results, function(r) {
    g <- r$genes
    g$sample <- r$sample
    g
  }))
  gene_long$is_target_pos <- gene_long$sample %in% ERBB2_POSITIVE
  write.table(
    gene_long[, c("sample", "gene", "is_target_pos",
                  "n_segments", "copy_ratio", "confidence")],
    file.path(RESULT_DIR, "all_genes_per_sample.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  ## Segment-level rows for the target gene.
  target_seg <- do.call(rbind, lapply(results, function(r) {
    s <- r$segments[r$segments$gene == TARGET_GENE, ]
    s$sample        <- r$sample
    s$is_target_pos <- r$sample %in% ERBB2_POSITIVE
    s$cn_abs        <- 2 * s$copy_ratio
    s$seg_idx       <- seq_len(nrow(s))
    s
  }))
  write.table(
    target_seg[, c("sample", "is_target_pos", "seg_idx",
                   "chr", "start", "end", "gene",
                   "depth", "anchor_ratio", "copy_ratio",
                   "z_score", "p_value", "neglog10_p", "cn_abs")],
    file.path(RESULT_DIR, "ERBB2_segment_level.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

  ## Per-sample target-gene summary.
  summary_tbl <- do.call(rbind, lapply(results, function(r) {
    g <- r$genes[r$genes$gene == TARGET_GENE, ]
    s <- r$segments[r$segments$gene == TARGET_GENE, ]
    cn_abs <- 2 * s$copy_ratio
    data.frame(
      sample        = r$sample,
      is_target_pos = r$sample %in% ERBB2_POSITIVE,
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
  write.table(summary_tbl, file.path(RESULT_DIR, "ERBB2_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n================ ERBB2 summary ================\n")
  print(summary_tbl, row.names = FALSE)

  ## -------------------------------------------------------------------------
  ## 4. Plots — same look as the hand-rolled demo so they're swappable
  ## -------------------------------------------------------------------------
  ## (a) per-sample bar — copy_ratio (gene-level, what the package returns).
  ##     Also annotate gene-level CN (= 2 * copy_ratio) for clinicians.
  out_bar <- file.path(FIG_DIR, "erbb2_per_sample.png")
  grDevices::png(out_bar, width = 1500, height = 700, res = 130)
  graphics::par(mar = c(7, 5, 4, 1))
  d <- summary_tbl
  cols <- ifelse(d$is_target_pos, "#d6604d", "#92c5de")
  ymax <- max(2 * d$copy_ratio, na.rm = TRUE) * 1.15
  bp <- graphics::barplot(2 * d$copy_ratio, names.arg = d$sample,
                          col = cols, las = 2, ylim = c(0, ymax),
                          ylab = sprintf("%s gene-level CN  (= 2 * copy_ratio from cesar_detect)",
                                          TARGET_GENE),
                          main = sprintf("A269 %s demo (pkg API) - per-sample copy number",
                                          TARGET_GENE))
  graphics::abline(h = 2, col = "grey50", lty = 3)
  graphics::text(bp, 2 * d$copy_ratio + ymax * 0.025,
                 labels = sprintf("conf=%.1f", d$conf), cex = 0.78)
  graphics::legend("topright", fill = c("#d6604d", "#92c5de"),
                   legend = c(sprintf("%s amplified (truth)", TARGET_GENE),
                              sprintf("%s negative (PoN)", TARGET_GENE)),
                   bty = "n")
  grDevices::dev.off()

  ## (b) per-segment heatmap, ERBB2 segments x 9 samples
  cn_mat <- do.call(rbind, lapply(sample_names, function(s) {
    rows <- target_seg[target_seg$sample == s, ]
    rows <- rows[order(rows$seg_idx), ]
    rows$cn_abs
  }))
  rownames(cn_mat) <- sample_names
  ord_h <- order(!sample_names %in% ERBB2_POSITIVE, sample_names)
  cn_mat <- cn_mat[ord_h, , drop = FALSE]
  cn_clip <- pmax(pmin(cn_mat, 5), 0)

  out_hm <- file.path(FIG_DIR, "erbb2_segment_heatmap.png")
  grDevices::png(out_hm, width = 1500, height = 700, res = 130)
  graphics::par(mar = c(5, 8, 4, 6))
  pal <- grDevices::colorRampPalette(
    c("#053061", "#2166ac", "white", "#d6604d", "#67001f"))(100)
  breaks <- seq(0, 5, length.out = 101)
  graphics::image(seq_len(ncol(cn_clip)), seq_len(nrow(cn_clip)),
                  t(cn_clip), col = pal, breaks = breaks,
                  axes = FALSE,
                  xlab = sprintf("%s segment index", TARGET_GENE), ylab = "",
                  main = sprintf("A269 %s demo (pkg API) - per-segment CN (clipped at 5)",
                                  TARGET_GENE))
  graphics::axis(1, at = seq(1, ncol(cn_clip), by = 2))
  graphics::axis(2, at = seq_len(nrow(cn_clip)),
                 labels = paste0(rownames(cn_clip),
                                 ifelse(rownames(cn_clip) %in% ERBB2_POSITIVE,
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

  ## (c) copy_ratio vs confidence scatter
  out_sc <- file.path(FIG_DIR, "erbb2_vs_conf.png")
  grDevices::png(out_sc, width = 1100, height = 800, res = 130)
  graphics::par(mar = c(5, 5, 4, 1))
  conf_clip <- pmin(d$conf, 60)
  graphics::plot(conf_clip, d$copy_ratio,
                 pch = 19, cex = 1.6,
                 xlim = c(-2, 68),
                 col = ifelse(d$is_target_pos, "#d6604d", "#92c5de"),
                 xlab = "ERBB2 confidence  (mean -log10(p) over ERBB2 segs, clipped at 60)",
                 ylab = "ERBB2 copy_ratio  (mu_train / r_obs; >1 amp, <1 loss)",
                 main = sprintf("A269 %s demo (pkg API) - separation of positives from PoN",
                                 TARGET_GENE))
  graphics::abline(h = 1,   col = "grey50", lty = 3)
  graphics::abline(v = 2,   col = "grey50", lty = 3)
  lab_pos <- ifelse(conf_clip > 50, 2, 4)
  graphics::text(conf_clip, d$copy_ratio,
                 labels = d$sample, pos = lab_pos,
                 offset = 0.5, cex = 0.75)
  graphics::legend("topleft", pch = 19,
                   col = c("#d6604d", "#92c5de"),
                   legend = c(sprintf("%s amplified", TARGET_GENE),
                              sprintf("%s negative", TARGET_GENE)),
                   bty = "n")
  grDevices::dev.off()

  cat(sprintf("\nDemo done. Tables in %s/\n            Figures in %s/\n",
              normalizePath(RESULT_DIR, winslash = "/"),
              normalizePath(FIG_DIR,    winslash = "/")))
  invisible(NULL)
}, file.path(RESULT_DIR, "log.txt"))
