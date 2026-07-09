#' Read a pileup file
#'
#' Reads a tab-separated depth/pileup file. Required columns: chromosome (1st),
#' 1-based position (2nd), depth (3rd). Additional columns (e.g. allele
#' counts in a \code{.freq} file) are preserved but not used by CesaR.
#' Plain-text and gzip-compressed (\code{.gz}) files are both accepted.
#'
#' @param path Path to the pileup file (plain or gzip-compressed).
#' @param header Logical. Whether the file has a header row. Default \code{TRUE}.
#'
#' @return A data frame whose first three columns are renamed to
#'   \code{CHR}, \code{POSITION}, \code{DEPTH}.
#'
#' @export
read_pileup <- function(path, header = TRUE) {
  if (!file.exists(path)) {
    stop("Pileup file not found: ", path, call. = FALSE)
  }
  df <- utils::read.table(path, sep = "\t", header = header,
                          stringsAsFactors = FALSE)
  if (ncol(df) < 3) {
    stop("Pileup file '", path, "' has fewer than 3 columns; expected ",
         "chr/position/depth.", call. = FALSE)
  }
  colnames(df)[1:3] <- c("CHR", "POSITION", "DEPTH")
  df
}

#' Read and validate a 4-column BED file
#'
#' CesaR uses a 4-column BED throughout: chromosome, start, end, gene-name.
#' This helper reads the file, checks ordering and column count, and
#' canonicalizes column names to \code{c("chr","start","end","gene")}.
#' Whitespace (mixed tabs and spaces) is accepted as the column separator.
#'
#' @param path Path to the BED file.
#'
#' @return A data frame with columns \code{chr}, \code{start}, \code{end},
#'   \code{gene}.
#'
#' @export
read_panel_bed <- function(path) {
  if (!file.exists(path)) {
    stop("BED file not found: ", path, call. = FALSE)
  }
  bed <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE,
                           comment.char = "#")
  if (ncol(bed) < 3) {
    stop("BED file must have at least 3 columns (chr/start/end); ",
         "CesaR additionally requires a 4th gene-name column.",
         call. = FALSE)
  }
  if (ncol(bed) == 3) {
    bed[, 4] <- "Unknown"
  }
  bed <- bed[, 1:4]
  colnames(bed) <- c("chr", "start", "end", "gene")
  bed$start <- as.integer(bed$start)
  bed$end   <- as.integer(bed$end)
  bed$gene  <- as.character(bed$gene)
  if (anyNA(bed$start) || anyNA(bed$end)) {
    stop("BED file has non-integer start/end coordinates.", call. = FALSE)
  }
  if (any(bed$end < bed$start)) {
    stop("BED file has rows with end < start.", call. = FALSE)
  }
  bed
}

#' Compute mean depth in each BED region from a pileup
#'
#' For every row in \code{bed}, returns the mean read depth across positions
#' from \code{pileup} that fall within the region. CesaR's core API works on
#' precomputed depth matrices, but most upstream pipelines emit base-level
#' pileups (\code{.depth}, \code{.freq}, ...). This function bridges the two
#' representations: build the depth matrix once with this helper, then pass
#' the matrix to \code{\link{cesar_train}} and \code{\link{cesar_detect}}.
#'
#' @param pileup Either a path to a pileup file or a data frame already loaded
#'   via [read_pileup()]. The first three columns must be CHR / POSITION /
#'   DEPTH.
#' @param bed A 4-column BED data frame (see [read_panel_bed()]) or a path
#'   to a BED file.
#'
#' @return Numeric vector of length \code{nrow(bed)}; entries are 0 for
#'   regions with no overlapping pileup positions.
#'
#' @export
bed_depth_in_pileup <- function(pileup, bed) {
  if (is.character(pileup)) pileup <- read_pileup(pileup)
  if (is.character(bed))    bed    <- read_panel_bed(bed)
  if (!is.data.frame(pileup)) {
    pileup <- as.data.frame(pileup, stringsAsFactors = FALSE)
  }
  colnames(pileup)[1:3] <- c("CHR", "POSITION", "DEPTH")

  pos_key   <- chr_to_numeric(pileup$CHR) * 1e10 + pileup$POSITION
  bed_start <- chr_to_numeric(bed[, 1])    * 1e10 + bed[, 2]
  bed_end   <- chr_to_numeric(bed[, 1])    * 1e10 + bed[, 3]

  results <- numeric(nrow(bed))
  for (i in seq_len(nrow(bed))) {
    idx <- which(pos_key >= bed_start[i] & pos_key <= bed_end[i])
    if (length(idx) > 0) {
      results[i] <- mean(pileup$DEPTH[idx])
    }
  }
  results
}

#' Convert chromosome strings to numeric for sortable position keys
#' @keywords internal
#' @noRd
chr_to_numeric <- function(chr_vec) {
  v <- as.character(chr_vec)
  v <- gsub("^chr", "", v, ignore.case = TRUE)
  v <- gsub("X", "23", v, ignore.case = TRUE)
  v <- gsub("Y", "24", v, ignore.case = TRUE)
  v <- gsub("M(T)?$", "25", v, ignore.case = TRUE)
  suppressWarnings(as.numeric(v))
}
