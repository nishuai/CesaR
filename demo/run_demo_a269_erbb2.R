## Cesar A269 ERBB2-amplification demo
##
## Real-cohort variant of the MET demo, retargeted to ERBB2.
## Truth: A269010005 / 006 / 007 / 009 carry ERBB2 amplifications;
## the remaining 5 (002 / 003 / 004 / 010 / 011) are ERBB2-negative
## and serve as the in-cohort panel-of-normals for the ERBB2 locus.
##
## Run from the project root with:
##   Rscript demo/run_demo_a269_erbb2.R
##
## Inputs:
##   A269/A269010*.sort.bam.mpileup.freq            -- 9 mpileup freq files
##   segmented_bed_A269_hg38_annotated.bed          -- 1645 segments, 16 genes
##                                                    (ERBB2 = 27 segments)
##
## Outputs:
##   demo/results_a269_erbb2/ERBB2_summary.tsv
##   demo/results_a269_erbb2/ERBB2_segment_level.tsv
##   demo/results_a269_erbb2/all_genes_per_sample.tsv
##   demo/results_a269_erbb2/log.txt
##   demo/results_a269_erbb2/model_anchors.rda
##   demo/results_a269_erbb2/model_parameters.rda
##   demo/figures_a269_erbb2/erbb2_per_sample.png
##   demo/figures_a269_erbb2/erbb2_segment_heatmap.png
##   demo/figures_a269_erbb2/erbb2_vs_conf.png

suppressPackageStartupMessages({
  library(MASS)
})
source("R/bed_depth_in_pileup.R")

## ---------------------------------------------------------------------------
## 0. Paths and ground truth
## ---------------------------------------------------------------------------
A269_DIR   <- "A269"
BED_FILE   <- "segmented_bed_A269_hg38_annotated.bed"
DEMO_DIR   <- "demo"
RESULT_DIR <- file.path(DEMO_DIR, "results_a269_erbb2")
FIG_DIR    <- file.path(DEMO_DIR, "figures_a269_erbb2")
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
  cat("================ A269 ERBB2 demo ================\n")
  cat("ERBB2-positive (4, held out for testing):\n  ",
      paste(ERBB2_POSITIVE, collapse = ", "), "\n", sep = "")
  cat("ERBB2-negative (5, used as PoN):\n  ",
      paste(ERBB2_NEGATIVE, collapse = ", "), "\n", sep = "")

  ## -------------------------------------------------------------------------
  ## 1. Read coverage matrix from per-sample .freq files
  ## -------------------------------------------------------------------------
  segment_bed <- read.table(BED_FILE, sep = "\t", header = FALSE,
                            stringsAsFactors = FALSE)
  freq_files  <- list.files(A269_DIR, pattern = "\\.freq$")
  sample_names <- sub("F01.sort.bam.mpileup.freq", "", freq_files)

  cat(sprintf("\nReading %d samples x %d segments from %s/\n",
              length(freq_files), nrow(segment_bed), A269_DIR))
  mat <- matrix(0, nrow = length(freq_files), ncol = nrow(segment_bed))
  for (i in seq_along(freq_files)) {
    cat("  ", sample_names[i], "\n")
    tab <- read.table(file.path(A269_DIR, freq_files[i]),
                      sep = "\t", header = TRUE)
    mat[i, ] <- bed_depth_in_pileup(tab, segment_bed)
  }
  rownames(mat) <- sample_names
  gene_names <- gsub("\\|.*", "", segment_bed[, 4])

  ## -------------------------------------------------------------------------
  ## 2. Train CESAR anchors on the 5 ERBB2-negative samples
  ##    Same-gene anchors are excluded so that a positive ERBB2 sample cannot
  ##    pull its own anchors with it during detection.
  ## -------------------------------------------------------------------------
  pon_idx <- which(sample_names %in% ERBB2_NEGATIVE)
  cat(sprintf("\nTraining on %d ERBB2-negative samples ...\n", length(pon_idx)))

  train <- mat[pon_idx, , drop = FALSE]
  min_ac <- 3; max_ac <- 25; min_cov <- 50
  low_cov <- colMeans(train) < min_cov
  train[, low_cov] <- 0

  cov_matrix <- suppressWarnings(cor(train))
  cov_matrix[is.na(cov_matrix)] <- 0
  top_k <- min(max_ac + 1, ncol(cov_matrix))
  max_cov_idx <- apply(cov_matrix, 2,
                       function(x) head(order(x, decreasing = TRUE), top_k))

  filter_anchors_by_gene <- function(idx_vec, target_gene) {
    idx_vec[gene_names[idx_vec] != target_gene]
  }

  gene_anchor_sd <- list()
  for (n_anchor in (min_ac + 1):max_ac) {
    sds <- sapply(seq_len(ncol(max_cov_idx)), function(i) {
      cand <- filter_anchors_by_gene(max_cov_idx[, i], gene_names[i])
      if (length(cand) < n_anchor) return(NA_real_)
      chosen <- cand[1:n_anchor]
      ratio <- rowMeans(train[, chosen[-1], drop = FALSE]) / train[, chosen[1]]
      ratio <- ratio[is.finite(ratio)]
      if (length(ratio) < 2 || mean(ratio) == 0) return(NA_real_)
      sd(ratio) / mean(ratio)
    })
    gene_anchor_sd[[n_anchor]] <- sds
  }
  gene_anchor_sd <- do.call(rbind, gene_anchor_sd)
  anchors <- apply(gene_anchor_sd, 2,
                   function(x) if (all(is.na(x))) NA_integer_ else which.min(x))

  model_anchors <- list(); model_parameters <- list()
  for (i in seq_along(anchors)) {
    if (is.na(anchors[i])) {
      model_anchors[[i]] <- NA
      model_parameters[[i]] <- c(NA_real_, NA_real_)
      next
    }
    cand <- filter_anchors_by_gene(max_cov_idx[, i], gene_names[i])
    idx  <- cand[2:(anchors[i] + min_ac)]
    model_anchors[[i]] <- idx
    a_mean <- rowMeans(train[, idx, drop = FALSE])
    ratio <- a_mean / train[, i]
    ratio <- ratio[is.finite(ratio) & ratio > 0]
    model_parameters[[i]] <- if (length(ratio) < 2) c(NA_real_, NA_real_)
                              else c(mean(ratio), sd(ratio))
  }
  saveRDS(model_anchors,    file.path(RESULT_DIR, "model_anchors.rda"))
  saveRDS(model_parameters, file.path(RESULT_DIR, "model_parameters.rda"))

  target_idx <- which(gene_names == TARGET_GENE)
  cat(sprintf("%s segments fitted: %d\n", TARGET_GENE, length(target_idx)))

  ## -------------------------------------------------------------------------
  ## 3. Detect on all 9 samples
  ## -------------------------------------------------------------------------
  cat("\nRunning detection on all 9 samples ...\n")
  results_long  <- list()
  target_seg_rows <- list()
  for (s in seq_along(sample_names)) {
    test_depth <- mat[s, ]
    cnv_seg  <- numeric(length(test_depth))
    pval_seg <- numeric(length(test_depth))
    for (i in seq_along(test_depth)) {
      if (is.na(model_anchors[[i]][1]) || test_depth[i] == 0 ||
          is.na(model_parameters[[i]][1]) || is.na(model_parameters[[i]][2]) ||
          model_parameters[[i]][2] == 0) {
        cnv_seg[i] <- NA_real_; pval_seg[i] <- NA_real_; next
      }
      cur_ratio <- mean(mat[s, model_anchors[[i]]]) / test_depth[i]
      z <- (cur_ratio - model_parameters[[i]][1]) / model_parameters[[i]][2]
      cnv_seg[i]  <- model_parameters[[i]][1] / cur_ratio
      pval_seg[i] <- max(2 * pnorm(-abs(z)), 1e-100)
    }
    df <- data.frame(gene = gene_names, cn_rel = cnv_seg, pval = pval_seg,
                     seg_idx = seq_along(cnv_seg), stringsAsFactors = FALSE)
    agg <- do.call(rbind, lapply(split(df, df$gene), function(d) {
      cn <- d$cn_rel[is.finite(d$cn_rel)]
      pv <- d$pval[is.finite(d$pval)]
      data.frame(
        gene          = d$gene[1],
        n_seg         = length(cn),
        cn_rel_median = if (length(cn)) median(cn) else NA_real_,
        cn_rel_mean   = if (length(cn)) mean(cn)   else NA_real_,
        cn_rel_max    = if (length(cn)) max(cn)    else NA_real_,
        cn_rel_min    = if (length(cn)) min(cn)    else NA_real_,
        pval_min      = if (length(pv)) min(pv)    else NA_real_,
        stringsAsFactors = FALSE)
    }))
    agg$cn_abs_median <- 2 * agg$cn_rel_median
    agg$cn_abs_mean   <- 2 * agg$cn_rel_mean
    agg$cn_abs_max    <- 2 * agg$cn_rel_max
    agg$cn_abs_min    <- 2 * agg$cn_rel_min
    agg$conf         <- -log10(agg$pval_min)
    agg$sample       <- sample_names[s]
    agg$is_target_pos <- sample_names[s] %in% ERBB2_POSITIVE
    results_long[[s]] <- agg

    tg_df <- df[df$gene == TARGET_GENE, ]
    tg_df$sample        <- sample_names[s]
    tg_df$is_target_pos <- sample_names[s] %in% ERBB2_POSITIVE
    tg_df$cn_abs        <- 2 * tg_df$cn_rel
    target_seg_rows[[s]] <- tg_df[, c("sample", "is_target_pos", "seg_idx",
                                      "gene", "cn_rel", "cn_abs", "pval")]
  }

  ## -------------------------------------------------------------------------
  ## 4. Write tables
  ## -------------------------------------------------------------------------
  long <- do.call(rbind, results_long)
  write.table(long[, c("sample", "gene", "is_target_pos", "n_seg",
                       "cn_rel_median", "cn_rel_mean", "cn_rel_max", "cn_rel_min",
                       "cn_abs_median", "cn_abs_mean", "cn_abs_max", "cn_abs_min",
                       "pval_min", "conf")],
              file.path(RESULT_DIR, "all_genes_per_sample.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  target_seg <- do.call(rbind, target_seg_rows)
  write.table(target_seg, file.path(RESULT_DIR, "ERBB2_segment_level.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  target_only <- long[long$gene == TARGET_GENE,
                      c("sample", "is_target_pos", "n_seg",
                        "cn_abs_median", "cn_abs_mean", "cn_abs_max", "cn_abs_min",
                        "pval_min", "conf")]
  target_only <- target_only[order(!target_only$is_target_pos,
                                    target_only$sample), ]
  write.table(target_only, file.path(RESULT_DIR, "ERBB2_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n================ ERBB2 summary ================\n")
  print(target_only, row.names = FALSE)

  ## -------------------------------------------------------------------------
  ## 5. Plots
  ## -------------------------------------------------------------------------
  ## (a) per-sample ERBB2 CN bars, sorted by truth then sample
  out_bar <- file.path(FIG_DIR, "erbb2_per_sample.png")
  grDevices::png(out_bar, width = 1500, height = 700, res = 130)
  graphics::par(mar = c(7, 5, 4, 1))
  ord <- order(!target_only$is_target_pos, target_only$sample)
  d   <- target_only[ord, ]
  cols <- ifelse(d$is_target_pos, "#d6604d", "#92c5de")
  ymax <- max(d$cn_abs_max, na.rm = TRUE) * 1.1
  bp <- graphics::barplot(d$cn_abs_median, names.arg = d$sample,
                          col = cols, las = 2, ylim = c(0, ymax),
                          ylab = sprintf("%s absolute copy number (median across %d segs)",
                                          TARGET_GENE, length(target_idx)),
                          main = sprintf("A269 %s demo - per-sample copy number",
                                          TARGET_GENE))
  graphics::arrows(bp, d$cn_abs_min, bp, d$cn_abs_max,
                   angle = 90, code = 3, length = 0.05, col = "grey30")
  graphics::abline(h = 2, col = "grey50", lty = 3)
  graphics::text(bp, d$cn_abs_max + 0.15,
                 labels = sprintf("conf=%.1f", d$conf), cex = 0.8)
  graphics::legend("topright", fill = c("#d6604d", "#92c5de"),
                   legend = c(sprintf("%s amplified (truth)", TARGET_GENE),
                              sprintf("%s negative (PoN)", TARGET_GENE)),
                   bty = "n")
  grDevices::dev.off()

  ## (b) per-segment heatmap, ERBB2 segments x 9 samples
  cn_mat <- do.call(rbind, lapply(seq_along(sample_names), function(s) {
    rows <- target_seg[target_seg$sample == sample_names[s], ]
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
                  main = sprintf("A269 %s demo - per-segment absolute CN (clipped at 5)",
                                  TARGET_GENE))
  graphics::axis(1, at = seq(1, ncol(cn_clip), by = 2))
  graphics::axis(2, at = seq_len(nrow(cn_clip)),
                 labels = paste0(rownames(cn_clip),
                                 ifelse(rownames(cn_clip) %in% ERBB2_POSITIVE,
                                        " *", "")),
                 las = 2, cex.axis = 0.85)
  graphics::box()
  ## Color bar
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

  ## (c) ERBB2 CN_median vs confidence scatter
  out_sc <- file.path(FIG_DIR, "erbb2_vs_conf.png")
  grDevices::png(out_sc, width = 1100, height = 800, res = 130)
  graphics::par(mar = c(5, 5, 4, 1))
  conf_clip <- pmin(target_only$conf, 60)
  ## Give the x-axis breathing room so labels don't clip on either edge.
  graphics::plot(conf_clip, target_only$cn_abs_median,
                 pch = 19, cex = 1.6,
                 xlim = c(-2, 68),
                 col = ifelse(target_only$is_target_pos, "#d6604d", "#92c5de"),
                 xlab = "ERBB2 confidence  (-log10(min p), clipped at 60)",
                 ylab = sprintf("%s absolute CN (median of %d segs)",
                                 TARGET_GENE, length(target_idx)),
                 main = sprintf("A269 %s demo - separation of positives from PoN",
                                 TARGET_GENE))
  graphics::abline(h = 2, col = "grey50", lty = 3)
  graphics::abline(v = 2, col = "grey50", lty = 3)
  ## Flip labels for samples sitting at the right edge so they don't clip.
  lab_pos <- ifelse(conf_clip > 50, 2, 4)
  graphics::text(conf_clip, target_only$cn_abs_median,
                 labels = target_only$sample, pos = lab_pos,
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
