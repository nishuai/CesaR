#' Re-segment a panel BED using CBS on a master depth profile
#'
#' Takes a coarse panel BED and the master depth from [build_master_depth()],
#' runs circular binary segmentation (CBS via \code{PSCBS::segmentByCBS}) on
#' every region longer than \code{min_region_length}, and returns a
#' fine-grained BED whose segment breaks reflect true coverage structure
#' rather than user-supplied exon boundaries.
#'
#' @section Algorithm:
#' For each row in \code{panel_bed}:
#' \enumerate{
#'   \item If \code{end - start < min_region_length}, keep it as one segment.
#'   \item Otherwise, extract depth values from \code{master_depth} that fall
#'     within \code{[start, end]}, run \code{PSCBS::segmentByCBS} with
#'     \code{min.width = cbs_min_width}, and apply four merge rules to the
#'     returned changepoints:
#'     \itemize{
#'       \item Adjacent segments whose mean ratio differs by less than
#'         \code{min_ratio_diff} are merged.
#'       \item Adjacent segments whose absolute mean difference is less than
#'         \code{total_mean * min_abs_diff_frac} are merged.
#'       \item Segments with fewer than \code{min_loci_per_segment} positions
#'         are merged with neighbors.
#'       \item Changepoints closer than \code{min_changepoint_gap} bp are
#'         collapsed.
#'     }
#'   \item Each output sub-segment inherits the \code{gene} label from its
#'     parent panel region.
#' }
#'
#' @param master_depth A data frame with columns \code{CHR}, \code{POSITION},
#'   \code{DEPTH}, typically from [build_master_depth()].
#' @param panel_bed A 4-column BED data frame (\code{chr, start, end, gene})
#'   or a path to a BED file.
#' @param min_region_length Integer. Regions shorter than this (bp) are kept
#'   as-is without calling CBS. Default 600.
#' @param min_ratio_diff Numeric in (0,1). Adjacent CBS segments whose mean
#'   depths differ by less than this fraction (e.g. 0.10 = 10%) are merged.
#'   Default 0.10.
#' @param min_abs_diff_frac Numeric. Adjacent segments whose absolute mean
#'   difference is smaller than \code{total_mean * min_abs_diff_frac} are
#'   merged. Default 1/20.
#' @param min_loci_per_segment Integer. Segments with fewer positions are
#'   merged. Default 5.
#' @param min_changepoint_gap Integer. Changepoints closer than this (bp)
#'   are collapsed. Default 50.
#' @param cbs_min_width Integer. Passed to \code{PSCBS::segmentByCBS} as
#'   \code{min.width}; the minimum number of consecutive positions a segment
#'   may span. Default 5.
#' @param verbose Logical. Print progress messages.
#'
#' @return A 4-column data frame \code{(chr, start, end, gene)} with one row
#'   per output segment. The number of rows is \eqn{\ge} \code{nrow(panel_bed)}.
#'
#' @examples
#' \dontrun{
#'   master <- build_master_depth(pon_dir = "pon/")
#'   panel  <- read_panel_bed("panel_raw.bed")
#'   seg_bed <- cesar_segment(master, panel)
#'   nrow(seg_bed)  # typically >> nrow(panel)
#' }
#' @export
cesar_segment <- function(master_depth,
                          panel_bed,
                          min_region_length    = 600,
                          min_ratio_diff       = 0.10,
                          min_abs_diff_frac    = 1/20,
                          min_loci_per_segment = 5,
                          min_changepoint_gap  = 50,
                          cbs_min_width        = 5,
                          verbose              = TRUE) {
  if (!requireNamespace("PSCBS", quietly = TRUE)) {
    stop("Package 'PSCBS' is required for cesar_segment().\n",
         "  Install it with: install.packages('PSCBS')",
         call. = FALSE)
  }
  if (is.character(panel_bed)) panel_bed <- read_panel_bed(panel_bed)
  if (!is.data.frame(master_depth) || !all(c("CHR","POSITION","DEPTH") %in%
                                            colnames(master_depth))) {
    stop("master_depth must have CHR/POSITION/DEPTH columns.", call. = FALSE)
  }

  ## Pre-compute total_mean for the abs-diff merge rule.
  total_mean <- mean(master_depth$DEPTH, na.rm = TRUE)

  ## Build a fast lookup: CHR -> sorted (POSITION, DEPTH) pairs.
  master_split <- split(master_depth[, c("POSITION","DEPTH")],
                        master_depth$CHR)
  master_split <- lapply(master_split, function(d) {
    ord <- order(d$POSITION)
    list(pos = d$POSITION[ord], depth = d$DEPTH[ord])
  })

  out_segments <- vector("list", nrow(panel_bed))

  for (i in seq_len(nrow(panel_bed))) {
    chr   <- as.character(panel_bed$chr[i])
    start <- as.integer(panel_bed$start[i])
    end   <- as.integer(panel_bed$end[i])
    gene  <- as.character(panel_bed$gene[i])

    region_len <- end - start
    if (region_len < min_region_length) {
      ## Too short -> keep as one segment
      out_segments[[i]] <- data.frame(chr = chr, start = start,
                                      end = end, gene = gene,
                                      stringsAsFactors = FALSE)
      next
    }

    ## Extract master depth positions within [start, end]
    if (!chr %in% names(master_split)) {
      ## No master depth for this chr -> keep as one segment
      if (verbose) {
        message("CesaR: no master depth for ", chr,
                " (region ", i, "); keeping unsegmented.")
      }
      out_segments[[i]] <- data.frame(chr = chr, start = start,
                                      end = end, gene = gene,
                                      stringsAsFactors = FALSE)
      next
    }
    m <- master_split[[chr]]
    idx <- which(m$pos >= start & m$pos <= end)
    if (length(idx) < cbs_min_width) {
      ## Not enough positions -> keep as one
      out_segments[[i]] <- data.frame(chr = chr, start = start,
                                      end = end, gene = gene,
                                      stringsAsFactors = FALSE)
      next
    }

    depth_vec <- m$depth[idx]
    pos_vec   <- m$pos[idx]

    ## Run CBS
    segs_raw <- tryCatch(
      PSCBS::segmentByCBS(depth_vec, min.width = cbs_min_width, verbose = 0),
      error = function(e) NULL
    )
    if (is.null(segs_raw) || nrow(segs_raw$output) == 0) {
      out_segments[[i]] <- data.frame(chr = chr, start = start,
                                      end = end, gene = gene,
                                      stringsAsFactors = FALSE)
      next
    }

    ## Apply merge rules to decide which changepoints to keep.
    changepoints <- find_changepoints_from_cbs(
      segs_raw, depth_vec, pos_vec,
      total_mean, min_ratio_diff, min_abs_diff_frac,
      min_loci_per_segment, min_changepoint_gap
    )

    ## Force the parent region's true start and end as anchors so that
    ## sub-segments tile the parent without gaps or overshoot. Internal
    ## CBS-derived changepoints sit between them.
    internal_breaks <- changepoints[changepoints > start &
                                    changepoints < end]
    boundaries <- sort(unique(c(start, internal_breaks, end)))

    ## Drop internal breaks that are too close to either parent boundary
    ## (rule 4 again, but enforced against the *true* start/end, not the
    ## first/last master-depth position which may be inside the region).
    if (length(boundaries) > 2L) {
      d_left  <- boundaries[-1L]            - boundaries[-length(boundaries)]
      keep    <- c(TRUE, d_left[-length(d_left)] >= min_changepoint_gap, TRUE)
      ## Also require the gap to the right boundary be >= min_changepoint_gap
      d_right <- c(d_left[-1L], 0)
      keep    <- keep & c(TRUE,
                          d_right[-length(d_right)] >= min_changepoint_gap,
                          TRUE)
      boundaries <- boundaries[keep]
    }

    ## Build sub-segments tiling [start, end].
    n_b <- length(boundaries)
    out_segments[[i]] <- data.frame(
      chr   = chr,
      start = boundaries[-n_b],
      end   = boundaries[-1L],
      gene  = gene,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, out_segments)
  rownames(result) <- NULL
  if (verbose) {
    message("CesaR: segmented ", nrow(panel_bed), " panel regions into ",
            nrow(result), " fine-grained segments.")
  }
  result
}

#' Apply the four merge rules from find_changepoints.R to CBS output.
#' Returns a vector of genomic coordinates (changepoints), including the
#' region's start and end.
#' @keywords internal
#' @noRd
find_changepoints_from_cbs <- function(segs_raw, depth_vec, pos_vec,
                                       total_mean, min_ratio_diff,
                                       min_abs_diff_frac,
                                       min_loci_per_segment,
                                       min_changepoint_gap) {
  loci <- segs_raw$segRows
  seg_output <- segs_raw$output

  ## Rule 1 & 2: adjacent segments differ by ratio > min_ratio_diff
  ##             AND abs diff > total_mean * min_abs_diff_frac
  if (nrow(seg_output) > 1L) {
    ratio_diff <- abs(seg_output$mean[2:nrow(seg_output)] /
                      seg_output$mean[1:(nrow(seg_output)-1L)] - 1)
    abs_diff   <- abs(seg_output$mean[2:nrow(seg_output)] -
                      seg_output$mean[1:(nrow(seg_output)-1L)])
    mergeable  <- (ratio_diff > min_ratio_diff) &
                  (abs_diff > total_mean * min_abs_diff_frac)

    ## Rule 3: segment must have at least min_loci_per_segment positions
    mergeable  <- mergeable & (seg_output$nbrOfLoci[-1L] >= min_loci_per_segment)

    ## Always keep the last segment's end
    mergeable  <- c(mergeable, TRUE)
    loci <- loci[mergeable, , drop = FALSE]
  }

  ## Build changepoint vector: start of first segment, end of each kept segment.
  changepoints <- c(loci$startRow[1L], loci$endRow)

  ## Map row indices back to genomic coordinates
  changepoints_pos <- pos_vec[changepoints]

  ## Rule 4: adjacent changepoints must be >= min_changepoint_gap bp apart
  ## (except the very first and very last, which are region boundaries)
  keep <- c(TRUE,
            diff(changepoints_pos[-length(changepoints_pos)]) >= min_changepoint_gap,
            TRUE)
  changepoints_pos <- changepoints_pos[keep]

  changepoints_pos
}
