#' Train a CesaR CNV model from a precomputed PoN depth matrix
#'
#' Fits, for every BED segment, a personalized set of "anchor" segments
#' whose depth co-varies with that segment across the panel of normals.
#' At detection time, ratios computed against these anchors recover the
#' true copy-number signal even when raw coverage drifts non-linearly
#' across the panel due to GC, capture, or library-prep effects.
#'
#' Unlike the older Cesar package, CesaR's training step accepts only a
#' precomputed depth matrix; build it yourself with [bed_depth_in_pileup()]
#' for every PoN sample, or use a higher-level driver. This keeps the core
#' API free of file-IO and lets users plug CesaR into any upstream pipeline.
#'
#' @section Algorithm:
#' Compute the segment-to-segment Pearson correlation across the PoN; for
#' each target segment, take the most positively correlated other segments
#' as anchor candidates (with same-gene candidates excluded by default,
#' via \code{exclude_same_gene}). The anchor count k is chosen between
#' \code{min_anchors+1} and \code{max_anchors} to minimize the coefficient
#' of variation of the anchor-recalibrated ratio
#' \eqn{r_j = \mathrm{mean}(X[, anchors]) / X[, j]} across PoN samples.
#' The fitted normal mean and SD of \eqn{r_j} are stored as the per-segment
#' detection parameters.
#'
#' @param depth_matrix Numeric matrix of size
#'   \eqn{n_{samples} \times n_{segments}} containing segment-mean depths
#'   from the panel of normals. Build with [bed_depth_in_pileup()] applied
#'   to each PoN pileup, or supply your own.
#' @param bed A 4-column BED data frame (\code{chr, start, end, gene}) with
#'   one row per column of \code{depth_matrix}. Typically the output of
#'   [cesar_segment()].
#' @param min_mean_depth Numeric. Segments whose mean PoN depth is below
#'   this threshold are excluded from anchor candidacy. Default 200.
#' @param min_anchors,max_anchors Integer bounds on the number of anchor
#'   segments per target. Defaults 80 and 100.
#' @param exclude_same_gene Logical. Default \code{TRUE}. Restricts anchor
#'   candidates for each target segment to segments whose gene label
#'   differs from the target's. Without this, a whole-gene amplification
#'   that covers every segment of a gene makes the gene-mates each other's
#'   most-correlated PoN partners, the anchor mean tracks the target's
#'   amplified depth, and the recalibrated ratio collapses toward 1.
#' @param verbose Logical. Print progress messages.
#'
#' @return An object of class \code{"cesar_model"} -- a list with components
#'   \code{bed}, \code{anchors}, \code{params}, \code{depth_matrix} and
#'   \code{settings}. Use [cesar_detect()] to apply it to a test sample.
#'
#' @seealso [cesar_detect()], [cesar_segment()], [bed_depth_in_pileup()]
#' @export
cesar_train <- function(depth_matrix,
                        bed,
                        min_mean_depth    = 200,
                        min_anchors       = 80,
                        max_anchors       = 100,
                        exclude_same_gene = TRUE,
                        verbose           = TRUE) {
  if (!is.matrix(depth_matrix) && !is.data.frame(depth_matrix)) {
    stop("`depth_matrix` must be a numeric matrix (samples x segments).",
         call. = FALSE)
  }
  if (max_anchors <= min_anchors + 1L) {
    stop("`max_anchors` must exceed `min_anchors + 1`.", call. = FALSE)
  }
  if (!is.data.frame(bed) ||
      !all(c("chr", "start", "end", "gene") %in% colnames(bed))) {
    stop("`bed` must be a 4-column data frame with columns ",
         "chr/start/end/gene (e.g. from cesar_segment()).", call. = FALSE)
  }
  if (ncol(depth_matrix) != nrow(bed)) {
    stop("ncol(depth_matrix) (", ncol(depth_matrix),
         ") must equal nrow(bed) (", nrow(bed), ").", call. = FALSE)
  }
  if (nrow(depth_matrix) < 5L) {
    stop("Need at least 5 PoN samples; got ", nrow(depth_matrix),
         ".", call. = FALSE)
  }

  X <- unname(as.matrix(depth_matrix))
  pon_files <- rownames(depth_matrix)

  if (verbose) {
    message("CesaR: training from depth matrix (", nrow(X),
            " samples x ", ncol(X), " segments) ...")
  }

  fit_anchors_from_matrix(
    X, bed,
    pon_files         = pon_files,
    min_mean_depth    = min_mean_depth,
    min_anchors       = min_anchors,
    max_anchors       = max_anchors,
    exclude_same_gene = exclude_same_gene,
    verbose           = verbose
  )
}

#' @keywords internal
#' @noRd
fit_anchors_from_matrix <- function(X, bed, pon_files,
                                    min_mean_depth, min_anchors,
                                    max_anchors, exclude_same_gene = TRUE,
                                    verbose) {
  ## Mask chronically low-coverage segments so they can't be picked as
  ## anchors and don't get evaluated as targets.
  low_cov <- colMeans(X) < min_mean_depth
  X[, low_cov] <- 0
  if (verbose && any(low_cov)) {
    message("CesaR: masked ", sum(low_cov), "/", ncol(X),
            " segments below mean depth ", min_mean_depth)
  }

  if (verbose) message("CesaR: computing segment correlation matrix ...")
  cov_matrix <- suppressWarnings(stats::cor(X))
  cov_matrix[is.na(cov_matrix)] <- 0

  target_genes <- strip_gene_suffix(bed$gene)
  candidates <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    ord <- order(cov_matrix[, j], decreasing = TRUE)
    if (exclude_same_gene) {
      keep <- target_genes[ord] != target_genes[j] | ord == j
      ord <- ord[keep]
    }
    candidates[[j]] <- utils::head(ord, max_anchors)
  }
  if (exclude_same_gene && verbose) {
    pool_sizes <- vapply(candidates, length, integer(1))
    message("CesaR: same-gene anchor exclusion ON (",
            length(unique(target_genes)),
            " gene labels; pool size per target: median=",
            stats::median(pool_sizes - 1L),
            ", min=", min(pool_sizes - 1L), ")")
  }

  if (verbose) message("CesaR: selecting optimal anchor count per segment ...")
  cv_grid <- matrix(NA_real_, nrow = max_anchors, ncol = ncol(X))
  for (j in seq_len(ncol(X))) {
    cand <- candidates[[j]]
    if (length(cand) < min_anchors + 2L) next
    upper <- min(max_anchors, length(cand))
    for (k in (min_anchors + 1L):upper) {
      anchor_mean <- rowMeans(X[, cand[2:k], drop = FALSE])
      r <- anchor_mean / X[, cand[1]]
      m <- mean(r)
      cv_grid[k, j] <- if (!is.finite(m) || m == 0) NA_real_ else stats::sd(r) / m
    }
  }
  best_k <- apply(cv_grid, 2, function(col) {
    if (all(is.na(col))) NA_integer_ else which.min(col)
  })

  anchors_out <- vector("list", ncol(X))
  params_out  <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    if (is.na(best_k[j]) || X[1, j] == 0) {
      anchors_out[[j]] <- NA_integer_
      params_out[[j]]  <- c(mean = NA_real_, sd = NA_real_)
      next
    }
    anchor_idx <- candidates[[j]][2:best_k[j]]
    anchors_out[[j]] <- as.integer(anchor_idx)
    anchor_mean <- rowMeans(X[, anchor_idx, drop = FALSE])
    ratio <- anchor_mean / X[, j]
    fit <- tryCatch(MASS::fitdistr(ratio, "normal")$estimate,
                    error = function(e) c(mean = mean(ratio),
                                          sd   = stats::sd(ratio)))
    params_out[[j]] <- as.numeric(fit)
    names(params_out[[j]]) <- c("mean", "sd")
  }

  if (verbose) {
    message("CesaR: training complete. ",
            sum(!is.na(best_k)), "/", ncol(X), " segments fitted.")
  }

  structure(
    list(
      bed          = bed,
      anchors      = anchors_out,
      params       = params_out,
      depth_matrix = X,
      pon_files    = pon_files,
      settings     = list(min_mean_depth    = min_mean_depth,
                          min_anchors       = as.integer(min_anchors),
                          max_anchors       = as.integer(max_anchors),
                          exclude_same_gene = isTRUE(exclude_same_gene))
    ),
    class = "cesar_model"
  )
}

#' @export
print.cesar_model <- function(x, ...) {
  n_seg <- nrow(x$bed)
  n_fit <- sum(!vapply(x$params, function(p) all(is.na(p)), logical(1)))
  genes <- unique(strip_gene_suffix(x$bed$gene))
  cat("<cesar_model>\n")
  cat("  Panel BED          : ", n_seg, " segments across ",
      length(genes), " unique gene labels\n", sep = "")
  cat("  PoN samples used   : ", length(x$pon_files), "\n", sep = "")
  cat("  Segments fitted    : ", n_fit, "/", n_seg, "\n", sep = "")
  cat("  Anchor count range : (", x$settings$min_anchors, "+1, ",
      x$settings$max_anchors, "]\n", sep = "")
  cat("  Min mean depth     : ", x$settings$min_mean_depth, "\n", sep = "")
  cat("  Same-gene excluded : ",
      isTRUE(x$settings$exclude_same_gene), "\n", sep = "")
  invisible(x)
}

#' Strip RefSeq/Ensembl transcript suffix from a gene-label vector.
#' @keywords internal
#' @noRd
strip_gene_suffix <- function(g) {
  g <- as.character(g)
  vapply(strsplit(g, "\\|NM_|\\|EN", perl = FALSE),
         function(parts) parts[[1]], character(1))
}
