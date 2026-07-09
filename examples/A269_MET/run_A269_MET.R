###############################################################
## run_A269_MET.R
##
## CesaR worked example — MET amplification in the A269 cohort.
## Demonstrates the full 4-step pipeline on real, full-size pileup data:
##
##   1. build_master_depth()  on 6 MET-negative samples
##   2. cesar_segment()       (CBS) on the collapsed hg38 panel
##   3. cesar_train()         on the MET-negative coverage matrix
##   4. cesar_detect()        on the 3 MET-positive samples
##
## Sample roles (per investigator annotation):
##   002 / 003 / 004                       : MET-amplified, test samples
##   005 / 006 / 007 / 009 / 010 / 011     : MET-negative, PoN
##
## Run from the repository root:
##   Rscript examples/A269_MET/run_A269_MET.R
###############################################################
library(CesaR)

## Resolve the example's own folder so this script works from any cwd.
EX_DIR    <- "examples/A269_MET"
DATA_DIR  <- file.path(EX_DIR, "data")
PILEUPS   <- file.path(DATA_DIR, "pileups")
## Pre-segmented hg38 panel from the legacy pipeline. We use it as the
## starting point: collapse adjacent same-gene rows back to a "raw" panel
## (each gene = one or a few long regions), then feed that through
## cesar_segment() to re-segment via CBS on the MET-negative master depth.
## The vendor's hg19 raw panel BED is NOT used -- its coordinates
## don't match the hg38 mpileup files.
SEG_BED   <- file.path(DATA_DIR, "segmented_bed_A269_hg38_annotated.bed")

MET_POSITIVE <- c("A269010002", "A269010003", "A269010004")
MET_NEGATIVE <- c("A269010005", "A269010006", "A269010007",
                  "A269010009", "A269010010", "A269010011")

## ---------------------------------------------------------------------------
## 1. Sample inventory
## ---------------------------------------------------------------------------
freq_files   <- list.files(PILEUPS, pattern = "\\.freq$", full.names = TRUE)
sample_names <- sub("F01\\.sort\\.bam\\.mpileup\\.freq$", "",
                    basename(freq_files))
n_samples    <- length(freq_files)

cat(sprintf("Total samples    : %d\n", n_samples))
cat(sprintf("MET-positive     : %s\n", paste(MET_POSITIVE, collapse = ", ")))
cat(sprintf("MET-negative PoN : %s\n", paste(MET_NEGATIVE, collapse = ", ")))

## ---------------------------------------------------------------------------
## 2. Read every pileup once (needed for master depth + coverage matrix)
## ---------------------------------------------------------------------------
cat("\n[reading pileups]\n")
pileups <- vector("list", n_samples)
names(pileups) <- sample_names
for (i in seq_len(n_samples)) {
  cat("  ", sample_names[i], "\n", sep = "")
  pileups[[i]] <- read_pileup(freq_files[i])
}

## ---------------------------------------------------------------------------
## 3. Step 1 — master depth from the 6 MET-negative samples
## ---------------------------------------------------------------------------
cat("\n[step 1] build_master_depth() on MET-negative PoN\n")
master <- build_master_depth(
  pon_pileups = pileups[MET_NEGATIVE],
  verbose     = FALSE
)
cat(sprintf("  -> %d positions, mean depth = %.0f\n",
            nrow(master), mean(master$DEPTH)))

## ---------------------------------------------------------------------------
## 4. Step 2 — re-segment the panel using CBS on the master depth
## ---------------------------------------------------------------------------
## Collapse adjacent same-gene rows of the hg38 segmented BED back to a
## "raw" panel — one or a few long regions per gene — so cesar_segment()
## has something to split. This mirrors the bonus walkthrough at the
## bottom of run_demo.R.
cat("\n[step 2] cesar_segment() on collapsed raw panel\n")
hg38_seg <- read_panel_bed(SEG_BED)

raw_rows <- list()
cur <- hg38_seg[1, ]
for (i in 2:nrow(hg38_seg)) {
  same_chr  <- hg38_seg$chr [i] == cur$chr
  same_gene <- hg38_seg$gene[i] == cur$gene
  adjacent  <- hg38_seg$start[i] - cur$end <= 1
  if (same_chr && same_gene && adjacent) {
    cur$end <- hg38_seg$end[i]
  } else {
    raw_rows[[length(raw_rows) + 1L]] <- cur
    cur <- hg38_seg[i, ]
  }
}
raw_rows[[length(raw_rows) + 1L]] <- cur
raw_bed <- do.call(rbind, raw_rows)
rownames(raw_bed) <- NULL

cat(sprintf("  collapsed:        %d hg38 segments -> %d raw regions\n",
            nrow(hg38_seg), nrow(raw_bed)))
cat(sprintf("  longest region:   %d bp\n",
            max(raw_bed$end - raw_bed$start)))

seg_bed <- cesar_segment(master, raw_bed, verbose = FALSE)
cat(sprintf("  after segment():  %d segments (+%d new breaks vs raw)\n",
            nrow(seg_bed), nrow(seg_bed) - nrow(raw_bed)))

## ---------------------------------------------------------------------------
## 5. Build the samples x segments depth matrix on the re-segmented BED
## ---------------------------------------------------------------------------
cat("\n[coverage] building matrix on segmented BED\n")
mat <- matrix(0, nrow = n_samples, ncol = nrow(seg_bed),
              dimnames = list(sample_names, NULL))
for (i in seq_len(n_samples)) {
  mat[i, ] <- bed_depth_in_pileup(pileups[[i]], seg_bed)
}

## ---------------------------------------------------------------------------
## 6. Step 3 — train on the 6 MET-negative samples
## ---------------------------------------------------------------------------
pon_rows <- sample_names %in% MET_NEGATIVE
cat(sprintf("\n[step 3] cesar_train() on %d MET-negative samples\n",
            sum(pon_rows)))
model <- cesar_train(
  depth_matrix      = mat[pon_rows, , drop = FALSE],
  bed               = seg_bed,
  min_mean_depth    = 50,
  min_anchors       = 3,
  max_anchors       = 25,
  exclude_same_gene = TRUE,
  verbose           = FALSE
)
print(model)

## ---------------------------------------------------------------------------
## 7. Step 4 — detect on the 3 MET-positive samples, one at a time
## ---------------------------------------------------------------------------
## Per-gene aggregator. Reports median / mean / max / min of copy_ratio
## and the smallest p-value across each gene's segments. Straight
## gene-by-gene loop (no tapply / no apply).
aggregate_genes <- function(res, sample_id, is_met_pos) {
  segs  <- res$segments
  genes <- sort(unique(segs$gene))
  n_g   <- length(genes)

  out <- data.frame(
    sample        = rep(sample_id, n_g),
    gene          = genes,
    is_MET_pos    = is_met_pos,
    n_seg         = integer(n_g),
    cn_rel_median = numeric(n_g),
    cn_rel_mean   = numeric(n_g),
    cn_rel_max    = numeric(n_g),
    cn_rel_min    = numeric(n_g),
    pval_min      = numeric(n_g),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes)) {
    rows <- segs[segs$gene == genes[i], ]
    cn   <- rows$copy_ratio[is.finite(rows$copy_ratio)]
    pv   <- rows$p_value   [is.finite(rows$p_value)]
    if (length(cn) == 0L) {
      out$n_seg[i]         <- 0L
      out$cn_rel_median[i] <- NA_real_
      out$cn_rel_mean[i]   <- NA_real_
      out$cn_rel_max[i]    <- NA_real_
      out$cn_rel_min[i]    <- NA_real_
    } else {
      out$n_seg[i]         <- length(cn)
      out$cn_rel_median[i] <- median(cn)
      out$cn_rel_mean[i]   <- mean(cn)
      out$cn_rel_max[i]    <- max(cn)
      out$cn_rel_min[i]    <- min(cn)
    }
    out$pval_min[i] <- if (length(pv) == 0L) NA_real_ else min(pv)
  }

  out$cn_abs_median <- 2 * out$cn_rel_median
  out$cn_abs_mean   <- 2 * out$cn_rel_mean
  out$cn_abs_max    <- 2 * out$cn_rel_max
  out$cn_abs_min    <- 2 * out$cn_rel_min
  out$conf          <- -log10(out$pval_min)
  out
}

cat("\n[step 4] cesar_detect() on MET-positive samples\n")

cat("  A269010002\n")
res_002 <- cesar_detect(model, mat["A269010002", ], sample_id = "A269010002")
agg_002 <- aggregate_genes(res_002, "A269010002", TRUE)

cat("  A269010003\n")
res_003 <- cesar_detect(model, mat["A269010003", ], sample_id = "A269010003")
agg_003 <- aggregate_genes(res_003, "A269010003", TRUE)

cat("  A269010004\n")
res_004 <- cesar_detect(model, mat["A269010004", ], sample_id = "A269010004")
agg_004 <- aggregate_genes(res_004, "A269010004", TRUE)

## ---------------------------------------------------------------------------
## 8. MET-only summary table
## ---------------------------------------------------------------------------
met_only <- rbind(agg_002[agg_002$gene == "MET", ],
                  agg_003[agg_003$gene == "MET", ],
                  agg_004[agg_004$gene == "MET", ])
met_only <- met_only[, c("sample", "is_MET_pos", "n_seg",
                         "cn_abs_median", "cn_abs_mean",
                         "cn_abs_max",    "cn_abs_min",
                         "pval_min",      "conf")]
rownames(met_only) <- NULL

cat("\n=== MET CesaR summary (positive samples) ===\n")
print(met_only, row.names = FALSE)
