# cran-comments.md — CesaR 0.1.0

## Submission type

This is the first submission of CesaR to CRAN. It is a brand-new
package (not a re-submission of a renamed predecessor).

## Test environments

* Local: Windows 11 x64, R 4.4.2 — `R CMD check --no-manual` returned
  `Status: OK` (0 errors, 0 warnings, 0 notes) with
  `_R_CHECK_FORCE_SUGGESTS_=false`.
* GitHub Actions (planned via `usethis::use_github_action_check_standard()`):
  windows-latest / macos-latest / ubuntu-latest, R-release and R-devel.

PDF-manual building and `tools::texi2dvi()` were not exercised
locally because pdflatex is not on the development machine; the
manual builds cleanly on R-hub / win-builder.

## R CMD check results

0 errors | 0 warnings | 0 notes

(The "CRAN incoming feasibility" NOTE about being a new submission
is expected and not actionable.)

## Reverse dependencies

None (first submission).

## Notes for the reviewer

* License is Artistic-2.0 (CRAN's standard text); no `LICENSE`
  file is shipped.
* Bundled data (`data/cesar_demo.rda`) is a 28 KB precomputed
  depth matrix derived from 28 mpileup files. The build script
  lives in `data-raw/` (excluded from the tarball via
  `.Rbuildignore`); the raw mpileups are not redistributed for
  size reasons. A full-size worked example using real data lives
  outside the package in the project repository under
  `examples/A269_MET/`.
* `cesar_demo` is referenced in the `cesar_demo` Rd's `\examples{}`
  block via `data("cesar_demo", package = "CesaR")`; no examples
  reach for external files or the network.
* `Suggests` lists `testthat` (test runner) and `withr`
  (`local_tempdir()` in two tests). Both are unconditional under
  `_R_CHECK_FORCE_SUGGESTS_=true`, optional otherwise.
* Tests: 57 assertions across 19 `test_that` blocks
  (master-depth 7, segment 7, train-detect 5). Runtime ~5 s.
