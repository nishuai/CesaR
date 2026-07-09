#' Detect CNVs in a test sample using a trained CesaR model
#'
#' Computes anchor-recalibrated depth ratios for every BED segment in
#' \code{depth}, converts them to per-segment Z-scores against the
#' training-time mean/SD, and aggregates by gene.
#'
#' Unlike the older Cesar package, CesaR's detect step accepts only a
#' precomputed depth vector. Build it with [bed_depth_in_pileup()] applied
#' to the test sample's pileup.
#'
#' @param model A \code{"cesar_model"} object from [cesar_train()].
#' @param depth Numeric vector of length \code{nrow(model$bed)} giving the
#'   test sample's mean depth in each BED segment.
#' @param sample_id Optional name for the sample. Defaults to
#'   \code{"(unnamed)"}.
#' @param min_p Lower clamp on per-segment p-values before
#'   \eqn{-\log_{10}} transformation, to avoid \code{Inf}.
#'   Default \code{1e-100}.
#'
#' @return An object of class \code{"cesar_result"} -- a list with components
#'   \code{segments}, \code{genes}, \code{sample}.
#'
#' @seealso [cesar_train()]
#' @export
cesar_detect <- function(model, depth, sample_id = NULL, min_p = 1e-100) {
  if (!inherits(model, "cesar_model")) {
    stop("`model` must be a cesar_model from cesar_train().", call. = FALSE)
  }
  if (!is.numeric(depth) || is.matrix(depth)) {
    stop("`depth` must be a numeric vector. Use bed_depth_in_pileup() to ",
         "convert a pileup to a depth vector.", call. = FALSE)
  }
  bed <- model$bed
  n_seg <- nrow(bed)
  if (length(depth) != n_seg) {
    stop("length(depth) (", length(depth), ") must equal nrow(model$bed) (",
         n_seg, ").", call. = FALSE)
  }

  sample_name <- if (!is.null(sample_id)) sample_id else "(unnamed)"

  copy_ratio   <- numeric(n_seg)
  anchor_ratio <- numeric(n_seg)
  z_score      <- numeric(n_seg)
  p_values     <- numeric(n_seg)

  ## Per-segment ratio formula (matches Ni et al. and Cesar 1.x):
  ##   r_obs = mean(anchor_depth) / segment_depth
  ##   z     = (r_obs - mu_train) / sd_train
  ##   p     = 2 * pnorm(-|z|)         (two-sided)
  ##   copy  = mu_train / r_obs        (so copy=1 ↔ no change)
  for (i in seq_len(n_seg)) {
    anc <- model$anchors[[i]]
    par <- model$params[[i]]
    if (depth[i] == 0 || any(is.na(anc)) || any(is.na(par))) {
      anchor_ratio[i] <- NA_real_
      copy_ratio[i]   <- NA_real_
      z_score[i]      <- NA_real_
      p_values[i]     <- NA_real_
      next
    }
    r_obs <- mean(depth[anc]) / depth[i]
    anchor_ratio[i] <- r_obs
    if (par[2] == 0 || !is.finite(par[2])) {
      z_score[i]    <- NA_real_
      p_values[i]   <- NA_real_
      copy_ratio[i] <- par[1] / r_obs
      next
    }
    z_score[i]    <- (r_obs - par[1]) / par[2]
    p_values[i]   <- max(2 * stats::pnorm(-abs(z_score[i])), min_p)
    copy_ratio[i] <- par[1] / r_obs
  }

  segs <- data.frame(
    chr          = bed$chr,
    start        = bed$start,
    end          = bed$end,
    gene         = strip_gene_suffix(bed$gene),
    depth        = depth,
    anchor_ratio = anchor_ratio,
    copy_ratio   = copy_ratio,
    z_score      = z_score,
    p_value      = p_values,
    neglog10_p   = -log10(p_values),
    stringsAsFactors = FALSE
  )

  by_gene <- stats::aggregate(
    segs[, c("copy_ratio", "neglog10_p")],
    by  = list(gene = segs$gene),
    FUN = function(v) mean(v, na.rm = TRUE)
  )
  n_segs <- as.numeric(table(segs$gene)[by_gene$gene])
  genes <- data.frame(
    gene        = by_gene$gene,
    copy_ratio  = by_gene$copy_ratio,
    confidence  = by_gene$neglog10_p,
    n_segments  = n_segs,
    stringsAsFactors = FALSE
  )
  genes <- genes[order(genes$confidence, decreasing = TRUE,
                       na.last = TRUE), ]
  rownames(genes) <- NULL

  structure(
    list(segments = segs, genes = genes, sample = sample_name),
    class = "cesar_result"
  )
}

#' @export
print.cesar_result <- function(x, n = 5L, ...) {
  cat("<cesar_result>\n")
  cat("  Sample : ", x$sample, "\n", sep = "")
  cat("  Segments scored : ",
      sum(!is.na(x$segments$p_value)), "/", nrow(x$segments), "\n", sep = "")
  cat("\nTop ", min(n, nrow(x$genes)), " genes by confidence:\n", sep = "")
  print(utils::head(x$genes, n), row.names = FALSE)
  invisible(x)
}

#' @export
summary.cesar_result <- function(object, ...) {
  print(object, n = 10L)
  invisible(object)
}
