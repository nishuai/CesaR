# CesaR 0.1.0

Initial release. CesaR rebuilds the full panel-segmentation -> anchor-training
-> per-sample CNV detection pipeline that ad-hoc scripts implemented (see
`R/cnv30k.R`, `R/find_changepoints.R`, `R/get_master_depth_from_normals.R` in
the project root) and the older `Cesar` package abbreviated.

## What's new versus `Cesar` 1.1.0

* **`build_master_depth()`** -- average per-base depth across a panel of
  normals (PoN). Three input modes: directory of pileups, explicit file
  paths, or a list of in-memory data frames. Fast `rowMeans()` path when
  all PoN samples share identical position layout, with an inner-join
  fallback for misaligned inputs.
* **`cesar_segment()`** -- re-segment a raw panel BED using CBS
  (`PSCBS::segmentByCBS`) on the master depth. Five tunable thresholds
  control the merge rules (ratio diff, abs diff, min loci per segment,
  changepoint gap, CBS min width). Region start/end are forced as anchors
  so sub-segments tile parent regions without gaps.
* **API discipline** -- both `cesar_train()` and `cesar_detect()` accept
  only precomputed depth matrices/vectors (no path inputs). Use
  `bed_depth_in_pileup()` to convert a pileup file to the depth
  representation.
* **`exclude_same_gene = TRUE`** by default in `cesar_train()` (was
  introduced in `Cesar` 1.1.0 with default `FALSE`). Restricts anchor
  candidates to segments whose gene label differs from the target's,
  which is the right default for whole-gene CNV studies and was always
  applied in the original hand-rolled `detect_a269_*` scripts.

## Bundled fixtures

The single bundled dataset is `data/cesar_demo.rda` -- a 28-sample
spike-in dilution series (10 normals + 18 spike-ins across 3 dilution
levels) reduced to per-segment mean depths against the 196-region
30k panel. Tests and the in-package demo (`inst/scripts/run_demo_full.R`)
operate entirely on this precomputed matrix. A full end-to-end run on
real, full-size pileup data lives outside the R package, under
`examples/A269_MET/` in the project repository.

## Verification

* 57 testthat assertions across 19 `test_that` blocks (master-depth 7,
  segment 7, train-detect 5).
* `R CMD check` 0 errors / 0 warnings / 3 environmental NOTEs.
* End-to-end demo on the bundled fixtures: positive sample
  `copy_ratio[ERBB2] = 2.32, conf = 100`; negative sample max conf < 1.4.
