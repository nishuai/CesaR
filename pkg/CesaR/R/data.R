#' CesaR spike-in dilution series — precomputed per-segment depths
#'
#' A compact, self-contained dataset built from every mpileup in the
#' original CESAR development cohort. Each sample is reduced to its
#' length-196 vector of mean depths against the bundled 30k panel BED,
#' so the entire 28-sample dataset (10 normals + 18 spike-ins across
#' three dilution levels) fits in ~28 KB.
#'
#' @format A named list with five elements:
#' \describe{
#'   \item{\code{bed}}{4-column BED data frame (\code{chr}, \code{start},
#'     \code{end}, \code{gene}) describing the 196 panel segments. Genes
#'     covered: EGFR (4 segs), MET (5 segs), ERBB2 (3 segs), \code{Other}
#'     (184 control segs).}
#'   \item{\code{depth_matrix}}{Numeric matrix, 28 samples x 196 segments,
#'     of segment-mean depths. Rows are named by \code{sample_id}, columns
#'     by \code{seg_NNN}. Suitable for direct use as the \code{depth_matrix}
#'     argument of \code{\link{cesar_train}} or as the \code{depth} argument
#'     of \code{\link{cesar_detect}} (one row at a time).}
#'   \item{\code{sample_meta}}{Per-sample metadata: \code{sample_id},
#'     \code{condition}, and the relative path to the source pileup file.}
#'   \item{\code{condition_meta}}{Per-condition truth: empirical CNV
#'     fold-change for MET and ERBB2 implied by each spike-in label.
#'     The \code{normal} row carries \code{1.0 / 1.0}.}
#'   \item{\code{description}, \code{built_on}}{Provenance strings.}
#' }
#'
#' @details
#' The four conditions encode the spike-in mixture's nominal copy-number
#' state for the two oncogenes:
#' \tabular{lrr}{
#'   condition              \tab MET   \tab ERBB2 \cr
#'   \code{normal}              \tab 1.000 \tab 1.000  \cr
#'   \code{met2.1875ERBB2.625}  \tab 2.188 \tab 0.625  \cr
#'   \code{met2.375ERBB3.25}    \tab 2.375 \tab 3.250  \cr
#'   \code{met2.75ERBB4.75}     \tab 2.750 \tab 4.750  \cr
#' }
#' The numbers are dilution-series targets, not absolute copy numbers in
#' a diploid genome — what the package is expected to recover is the
#' monotonic dose-response, not the literal multipliers (the ctDNA
#' fraction also enters the relationship).
#'
#' Because the depth matrix is precomputed, \code{cesar_demo} bypasses the
#' first two pipeline steps (\code{\link{build_master_depth}} and
#' \code{\link{cesar_segment}}). A full end-to-end walkthrough that
#' exercises those steps on real, full-size pileup data lives outside
#' the R package, under \code{examples/A269_MET/} in the project
#' repository.
#'
#' Build the dataset from raw mpileups with
#' \code{Rscript data-raw/build_cesar_demo.R} (the script lives in the
#' package source, not in the installed copy).
#'
#' @examples
#' data("cesar_demo", package = "CesaR")
#' table(cesar_demo$sample_meta$condition)
#'
#' # Train on the 10 normals, detect on a known-positive spike-in sample.
#' # The bundled depth matrix is trimmed (rows summed only to the panel),
#' # so we relax `min_mean_depth` from its default 200 to fit all segments.
#' normals <- cesar_demo$sample_meta$condition == "normal"
#' model <- cesar_train(
#'   depth_matrix      = cesar_demo$depth_matrix[normals, ],
#'   bed               = cesar_demo$bed,
#'   min_mean_depth    = 50,
#'   exclude_same_gene = TRUE,
#'   verbose           = FALSE
#' )
#' pos_idx <- which(cesar_demo$sample_meta$sample_id == "703-1")
#' res <- cesar_detect(model,
#'                     cesar_demo$depth_matrix[pos_idx, ],
#'                     sample_id = "703-1 (MET~2.75x, ERBB2~4.75x)")
#' summary(res)
#'
#' @source
#' Pileup files generated in-house at the BGI Tumor R&D Department
#' (now Huazhong University of Science and Technology) from a
#' commercial spike-in dilution panel run on the 30k targeted ctDNA assay.
"cesar_demo"
