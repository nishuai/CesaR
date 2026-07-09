## Cesar A269 EGFR exon-19 deletion (19del) demo
##
## Detect EGFR 19del from per-base mpileup `.freq` files by exploiting an
## indirect indel footprint: samtools mpileup counts an indel-spanning read
## in DEPTH but not in REF / ALT base columns, so
##    deficit = DEPTH - (REF+ + REF- + A+ + A- + C+ + C- + T+ + T- + G+ + G-)
## equals the number of reads supporting *any* indel at that base.
##
## A 9-24 bp deletion shows up as a contiguous run of bases where every
## position has the same non-zero deficit (the indel-supporting reads
## carry the same sequence drop-out across the deleted span).
##
## Cohort (9 ctDNA samples on the QuarStar 94-gene panel):
##   - EGFR 19del positive (truth):  A269010010, A269010011
##   - EGFR 19del negative (PoN):    A269010002..007, A269010009  (n=7; 008 absent)
##
## Run from the project root with:
##   Rscript demo/run_demo_a269_egfr19del.R
##
## Outputs:
##   demo/results_a269_egfr19del/exon19_summary.tsv
##   demo/results_a269_egfr19del/exon19_per_base.tsv
##   demo/results_a269_egfr19del/log.txt
##   demo/figures_a269_egfr19del/exon19_deficit_traces.png
##   demo/figures_a269_egfr19del/exon19_run_length_bar.png

DEMO_DIR   <- "demo"
RESULT_DIR <- file.path(DEMO_DIR, "results_a269_egfr19del")
FIG_DIR    <- file.path(DEMO_DIR, "figures_a269_egfr19del")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,    showWarnings = FALSE, recursive = TRUE)

A269_DIR  <- "A269"
## EGFR exon 19, hg38: chr7:55174722-55174820 (the 19del hotspot proper),
## flanked by ~20 bp on each side from segment chr7:55174704-55174823.
EGFR_EX19_CHR   <- "chr7"
EGFR_EX19_START <- 55174704L
EGFR_EX19_END   <- 55174823L
HOTSPOT_START   <- 55174722L
HOTSPOT_END     <- 55174820L

POSITIVES <- c("A269010010", "A269010011")
NEGATIVES <- c("A269010002", "A269010003", "A269010004",
               "A269010005", "A269010006", "A269010007", "A269010009")

cap <- function(expr, file) {
  con <- file(file, open = "w", encoding = "UTF-8")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(); close(con) })
  invisible(force(expr))
}

cap({
  cat("================ A269 EGFR 19del demo ================\n")
  cat("Hotspot     : ", EGFR_EX19_CHR, ":", HOTSPOT_START, "-", HOTSPOT_END,
      " (hg38, EGFR exon 19)\n", sep = "")
  cat("Positives   : ", paste(POSITIVES, collapse = ", "), "\n", sep = "")
  cat("Negatives   : ", paste(NEGATIVES, collapse = ", "), "\n", sep = "")

  ## -------------------------------------------------------------------------
  ## 1. Read each sample's exon-19 base-level slice
  ## -------------------------------------------------------------------------
  freq_files   <- list.files(A269_DIR, pattern = "\\.freq$", full.names = TRUE)
  sample_names <- sub("F01\\.sort\\.bam\\.mpileup\\.freq$", "",
                      basename(freq_files))

  read_exon19 <- function(path) {
    df <- read.table(path, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
    keep <- df[[1]] == EGFR_EX19_CHR &
            df[[2]] >= EGFR_EX19_START & df[[2]] <= EGFR_EX19_END
    df <- df[keep, , drop = FALSE]
    pos       <- df[[2]]
    depth     <- df[[3]]
    base_sum  <- rowSums(df[, 5:14, drop = FALSE])
    deficit   <- depth - base_sum
    deficit[deficit < 0] <- 0L
    data.frame(pos = pos, depth = depth, base_sum = base_sum,
               deficit = deficit, vaf = ifelse(depth > 0, deficit/depth, 0),
               stringsAsFactors = FALSE)
  }

  per_sample <- setNames(lapply(freq_files, read_exon19), sample_names)

  ## Per-base table for the markdown / downstream plots
  per_base <- do.call(rbind, lapply(sample_names, function(s) {
    d <- per_sample[[s]]
    data.frame(sample = s, status = ifelse(s %in% POSITIVES, "positive",
                                    ifelse(s %in% NEGATIVES, "negative",
                                           "unlabelled")),
               pos = d$pos, depth = d$depth, deficit = d$deficit,
               vaf  = round(d$vaf, 4),
               stringsAsFactors = FALSE)
  }))
  write.table(per_base, file.path(RESULT_DIR, "exon19_per_base.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  ## -------------------------------------------------------------------------
  ## 2. Per-sample summary: longest contiguous run of (deficit >= cutoff)
  ##    inside the hotspot, and the median VAF over that run.
  ##
  ## Cutoff strategy: build a per-base background from the 7 negatives,
  ## then call any position where deficit > max(2, 99th percentile of
  ## negatives at that position) a "hit". This adapts to local coverage
  ## differences without hard-coding a VAF threshold.
  ## -------------------------------------------------------------------------
  pos_axis <- per_sample[[sample_names[1]]]$pos
  neg_def_mat <- do.call(rbind, lapply(NEGATIVES,
                                       function(s) per_sample[[s]]$deficit))
  per_pos_p99 <- apply(neg_def_mat, 2,
                       function(x) quantile(x, 0.99, names = FALSE, na.rm = TRUE))
  per_pos_cutoff <- pmax(2, per_pos_p99 + 1)  # require at least 2 reads + headroom

  longest_run <- function(hit_vec) {
    if (!any(hit_vec)) return(c(start = NA_integer_, end = NA_integer_,
                                 length = 0L))
    rl <- rle(hit_vec)
    runs_true <- which(rl$values)
    if (!length(runs_true)) return(c(start = NA_integer_, end = NA_integer_,
                                      length = 0L))
    cum <- cumsum(rl$lengths)
    start <- cum - rl$lengths + 1
    best  <- runs_true[which.max(rl$lengths[runs_true])]
    c(start = start[best], end = cum[best], length = rl$lengths[best])
  }

  hotspot_mask <- pos_axis >= HOTSPOT_START & pos_axis <= HOTSPOT_END

  summary_rows <- lapply(sample_names, function(s) {
    d <- per_sample[[s]]
    hit <- d$deficit > per_pos_cutoff
    hit_in_hotspot <- hit & hotspot_mask
    run <- longest_run(hit_in_hotspot)
    if (run["length"] > 0L) {
      run_pos_start <- d$pos[run["start"]]
      run_pos_end   <- d$pos[run["end"]]
      run_def       <- d$deficit[run["start"]:run["end"]]
      run_dep       <- d$depth  [run["start"]:run["end"]]
      run_vaf       <- median(run_def / run_dep)
      run_def_med   <- median(run_def)
      run_depth_med <- median(run_dep)
    } else {
      run_pos_start <- NA_integer_; run_pos_end <- NA_integer_
      run_vaf <- 0; run_def_med <- 0; run_depth_med <- median(d$depth[hotspot_mask])
    }
    data.frame(
      sample          = s,
      status          = ifelse(s %in% POSITIVES, "positive",
                        ifelse(s %in% NEGATIVES, "negative", "unlabelled")),
      hotspot_depth_median = round(median(d$depth[hotspot_mask]), 0),
      max_deficit_in_hotspot = max(d$deficit[hotspot_mask]),
      run_length_bp   = unname(run["length"]),
      run_start       = run_pos_start,
      run_end         = run_pos_end,
      run_deficit_med = run_def_med,
      run_vaf         = round(run_vaf, 4),
      call_19del      = unname(run["length"]) >= 9L,  # 19del: 9-24 bp typical
      stringsAsFactors = FALSE
    )
  })
  summary_tbl <- do.call(rbind, summary_rows)
  summary_tbl <- summary_tbl[order(summary_tbl$status != "positive",
                                    summary_tbl$sample), ]
  write.table(summary_tbl, file.path(RESULT_DIR, "exon19_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("\n================ EGFR 19del summary ================\n")
  print(summary_tbl, row.names = FALSE)

  ## -------------------------------------------------------------------------
  ## 3. Plots
  ## -------------------------------------------------------------------------
  ## (a) per-sample deficit trace across exon 19 (overlay)
  out_tr <- file.path(FIG_DIR, "exon19_deficit_traces.png")
  grDevices::png(out_tr, width = 1600, height = 800, res = 130)
  graphics::par(mar = c(5, 5, 4, 1))
  ymax <- max(sapply(per_sample, function(d) max(d$vaf))) * 1.1
  graphics::plot(NA, xlim = c(EGFR_EX19_START, EGFR_EX19_END),
                 ylim = c(0, ymax),
                 xlab = "chr7 position (hg38)",
                 ylab = "indel-supporting fraction  (deficit / DEPTH)",
                 main = "EGFR exon 19 - per-base indel signal across A269 cohort")
  graphics::abline(v = c(HOTSPOT_START, HOTSPOT_END),
                   col = "grey70", lty = 3)
  graphics::mtext(sprintf("hotspot %d-%d", HOTSPOT_START, HOTSPOT_END),
                  side = 3, at = (HOTSPOT_START + HOTSPOT_END) / 2,
                  cex = 0.85, col = "grey40")
  for (s in sample_names) {
    col <- if (s %in% POSITIVES) "#d6604d"
           else                  "#cccccc"
    lwd <- if (s %in% POSITIVES) 2.4 else 1.2
    graphics::lines(per_sample[[s]]$pos, per_sample[[s]]$vaf,
                    col = col, lwd = lwd)
  }
  graphics::legend("topleft", lwd = c(2.4, 1.2),
                   col = c("#d6604d", "#cccccc"),
                   legend = c("EGFR 19del positive (truth)",
                              "EGFR 19del negative (PoN)"),
                   bty = "n")
  ## Tag the called runs with sample names
  for (i in seq_len(nrow(summary_tbl))) {
    r <- summary_tbl[i, ]
    if (isTRUE(r$call_19del)) {
      graphics::text(r$run_start, r$run_vaf + ymax * 0.04,
                     labels = sprintf("%s  %d bp  VAF %.1f%%",
                                       r$sample, r$run_length_bp,
                                       100 * r$run_vaf),
                     col = "#67001f", pos = 4, cex = 0.85)
      graphics::segments(r$run_start, r$run_vaf, r$run_end, r$run_vaf,
                         col = "#67001f", lwd = 3)
    }
  }
  grDevices::dev.off()

  ## (b) bar of longest contiguous run length, colored by status
  out_bar <- file.path(FIG_DIR, "exon19_run_length_bar.png")
  grDevices::png(out_bar, width = 1500, height = 600, res = 130)
  graphics::par(mar = c(7, 5, 4, 1))
  ord <- order(summary_tbl$status != "positive", summary_tbl$sample)
  d   <- summary_tbl[ord, ]
  cols <- ifelse(d$status == "positive", "#d6604d", "#92c5de")
  bp <- graphics::barplot(d$run_length_bp, names.arg = d$sample,
                          col = cols, las = 2,
                          ylim = c(0, max(d$run_length_bp, 24) * 1.25),
                          ylab = "Longest contiguous indel-supporting run (bp)",
                          main = "A269 EGFR 19del demo - per-sample call length")
  graphics::abline(h = 9, col = "grey50", lty = 3)
  graphics::mtext("19del minimum length (9 bp)", side = 4, at = 9,
                  las = 1, cex = 0.75, col = "grey40", line = -2)
  graphics::text(bp, d$run_length_bp + 1.5,
                 labels = ifelse(d$run_length_bp > 0L,
                                 sprintf("%d bp\nVAF %.1f%%",
                                          d$run_length_bp, 100 * d$run_vaf),
                                 "ND"),
                 cex = 0.78)
  graphics::legend("topright", fill = c("#d6604d", "#92c5de"),
                   legend = c("EGFR 19del positive (truth)",
                              "EGFR 19del negative (PoN)"),
                   bty = "n")
  grDevices::dev.off()

  cat(sprintf("\nDemo done. Tables in %s/\n            Figures in %s/\n",
              normalizePath(RESULT_DIR, winslash = "/"),
              normalizePath(FIG_DIR,    winslash = "/")))
  invisible(NULL)
}, file.path(RESULT_DIR, "log.txt"))
