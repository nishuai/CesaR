#' CesaR: CNV detection with CBS panel segmentation and anchor recalibration
#'
#' CesaR rebuilds the full panel-segmentation -> anchor-training -> per-sample
#' detection pipeline that the original ad-hoc scripts implemented and the
#' \code{Cesar} package abbreviated. The user-facing flow is four steps:
#'
#' \enumerate{
#'   \item \code{\link{build_master_depth}()} -- average per-base depth across
#'     a panel of normals.
#'   \item \code{\link{cesar_segment}()} -- re-segment a raw panel BED using
#'     CBS (\code{PSCBS::segmentByCBS}) on the master depth, producing a
#'     fine-grained segmented BED whose breaks reflect coverage structure.
#'   \item \code{\link{cesar_train}()} -- fit per-segment correlation-based
#'     anchors and ratio distributions on the PoN. Same-gene anchor exclusion
#'     is on by default.
#'   \item \code{\link{cesar_detect}()} -- apply the model to a test sample.
#' }
#'
#' Inputs at every step are deliberately scrubbed of absolute paths: master
#' depth is built from a directory of pileups or a precomputed depth matrix;
#' the BED is a 4-column \code{(chr, start, end, gene)} data frame.
#'
#' @keywords internal
"_PACKAGE"
