#' Build a per-position master depth profile from a panel of normals
#'
#' Averages read depth at each genomic position across a set of normal
#' samples, producing the \emph{master depth} that drives the panel
#' re-segmentation step (\code{\link{cesar_segment}}). The result is the
#' coverage profile a "typical" sample would have under no copy-number
#' alteration, with random per-base noise smoothed out by averaging
#' across PoN samples.
#'
#' @section Input modes:
#' Pick exactly one of the following:
#' \itemize{
#'   \item \code{pon_dir} -- a directory; CesaR lists files matching
#'     \code{pileup_pattern}, reads each via \code{\link{read_pileup}}.
#'   \item \code{pon_files} -- a character vector of explicit pileup paths.
#'   \item \code{pon_pileups} -- a list of already-loaded data frames; each
#'     must have \code{CHR} / \code{POSITION} / \code{DEPTH} as the first
#'     three columns.
#' }
#'
#' @section Position alignment:
#' If every input has identical \code{(CHR, POSITION)} rows in the same
#' order (the common case for matched PoN files from one panel run), CesaR
#' takes a fast \code{rowMeans()} path. Otherwise it falls back to an inner
#' merge so that misaligned positions do not silently corrupt the average.
#'
#' @param pon_dir Directory containing PoN pileup files. Mutually exclusive
#'   with \code{pon_files} and \code{pon_pileups}.
#' @param pileup_pattern Regex used to filter files inside \code{pon_dir}.
#'   Default matches \code{.freq}, \code{.depth} and \code{.depth.gz}.
#' @param pon_files Character vector of explicit pileup paths (alternative to
#'   \code{pon_dir}).
#' @param pon_pileups List of in-memory pileup data frames (alternative to
#'   the path-based modes).
#' @param drop_zero_depth Logical; if \code{TRUE} (the default), positions
#'   whose averaged depth is exactly zero are dropped from the returned
#'   master profile. They cannot anchor a CBS segmentation and only inflate
#'   memory.
#' @param verbose Logical; print progress messages.
#'
#' @return A data frame with columns \code{CHR} (character), \code{POSITION}
#'   (integer) and \code{DEPTH} (numeric, the cross-sample mean), one row
#'   per genomic position present in the PoN.
#'
#' @examples
#' \dontrun{
#'   master <- build_master_depth(pon_dir = "path/to/pon",
#'                                pileup_pattern = "\\.freq$")
#'   head(master)
#' }
#' @export
build_master_depth <- function(pon_dir         = NULL,
                               pileup_pattern  = "\\.freq$|\\.depth(\\.gz)?$",
                               pon_files       = NULL,
                               pon_pileups     = NULL,
                               drop_zero_depth = TRUE,
                               verbose         = TRUE) {
  ## ---- 1. Resolve to a list of pileup data frames -------------------------
  modes <- c(!is.null(pon_dir), !is.null(pon_files), !is.null(pon_pileups))
  if (sum(modes) != 1L) {
    stop("Provide exactly one of `pon_dir`, `pon_files`, or `pon_pileups`.",
         call. = FALSE)
  }

  if (!is.null(pon_dir)) {
    if (!dir.exists(pon_dir)) {
      stop("Panel-of-normals directory not found: ", pon_dir, call. = FALSE)
    }
    pon_files <- list.files(pon_dir, pattern = pileup_pattern,
                            full.names = TRUE)
    if (length(pon_files) == 0L) {
      stop("No pileup files matched pattern '", pileup_pattern,
           "' in ", pon_dir, call. = FALSE)
    }
  }

  if (!is.null(pon_files)) {
    if (length(pon_files) < 2L) {
      stop("Need at least 2 PoN samples to average. Got ",
           length(pon_files), ".", call. = FALSE)
    }
    if (verbose) {
      message("CesaR: reading ", length(pon_files), " PoN pileups ...")
    }
    pon_pileups <- lapply(pon_files, function(p) {
      if (verbose) message("  ", basename(p))
      read_pileup(p)
    })
  } else {
    if (length(pon_pileups) < 2L) {
      stop("Need at least 2 PoN samples to average. Got ",
           length(pon_pileups), ".", call. = FALSE)
    }
    pon_pileups <- lapply(pon_pileups, function(d) {
      if (!is.data.frame(d) || ncol(d) < 3L) {
        stop("Each entry of `pon_pileups` must be a data frame with at ",
             "least CHR/POSITION/DEPTH columns.", call. = FALSE)
      }
      d <- as.data.frame(d, stringsAsFactors = FALSE)
      colnames(d)[1:3] <- c("CHR", "POSITION", "DEPTH")
      d
    })
  }

  ## ---- 2. Try the fast path: identical (CHR, POSITION) rows ---------------
  same_layout <- identical_position_layout(pon_pileups)

  if (same_layout) {
    if (verbose) {
      message("CesaR: position layout matches across PoN -> rowMeans path")
    }
    chr <- as.character(pon_pileups[[1L]]$CHR)
    pos <- as.integer(pon_pileups[[1L]]$POSITION)
    ## do.call(cbind, ...) keeps the matrix shape even when nrow == 1,
    ## which vapply() with numeric(1) does not.
    depth_mat <- do.call(cbind, lapply(pon_pileups,
                                       function(d) as.numeric(d$DEPTH)))
    master_depth <- rowMeans(depth_mat)
  } else {
    if (verbose) {
      message("CesaR: position layouts differ -> inner-join path")
    }
    merged <- Reduce(function(acc, d) {
      cur <- d[, c("CHR", "POSITION", "DEPTH")]
      colnames(cur)[3] <- paste0("DEPTH.", ncol(acc) - 1L)
      merge(acc, cur, by = c("CHR", "POSITION"), sort = FALSE)
    }, pon_pileups[-1L],
       init = pon_pileups[[1L]][, c("CHR", "POSITION", "DEPTH")])
    chr <- as.character(merged$CHR)
    pos <- as.integer(merged$POSITION)
    depth_mat <- as.matrix(merged[, -(1:2), drop = FALSE])
    master_depth <- rowMeans(depth_mat)
  }

  out <- data.frame(CHR      = chr,
                    POSITION = pos,
                    DEPTH    = master_depth,
                    stringsAsFactors = FALSE)

  ## Sort by chr (numeric ranking via chr_to_numeric) then position so that
  ## downstream CBS receives a tidy genome-ordered profile.
  ord <- order(chr_to_numeric(out$CHR), out$POSITION)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  if (isTRUE(drop_zero_depth)) {
    n_before <- nrow(out)
    out <- out[out$DEPTH > 0, , drop = FALSE]
    if (verbose && nrow(out) < n_before) {
      message("CesaR: dropped ", n_before - nrow(out),
              " zero-depth positions (", nrow(out), " kept).")
    }
    rownames(out) <- NULL
  }

  out
}

#' Test whether all pileups share identical (CHR, POSITION) rows in order.
#' @keywords internal
#' @noRd
identical_position_layout <- function(pileups) {
  if (length(pileups) < 2L) return(TRUE)
  chr0 <- as.character(pileups[[1L]]$CHR)
  pos0 <- as.integer(pileups[[1L]]$POSITION)
  for (k in 2:length(pileups)) {
    if (length(pileups[[k]]$POSITION) != length(pos0)) return(FALSE)
    if (!identical(as.integer(pileups[[k]]$POSITION), pos0)) return(FALSE)
    if (!identical(as.character(pileups[[k]]$CHR), chr0)) return(FALSE)
  }
  TRUE
}
