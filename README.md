# CesaR

> **C**NV **E**stimation with CBS panel **S**egmentation and **A**nchor **R**ecalibration

<!-- badges: start -->
![R](https://img.shields.io/badge/R-%3E%3D3.5.0-276DC3?logo=r&logoColor=white)
![Version](https://img.shields.io/badge/version-0.1.0-blue)
![License](https://img.shields.io/badge/license-Artistic--2.0-green)
![R CMD check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)
<!-- badges: end -->

CesaR detects **copy-number variations (CNVs)** in targeted sequencing panels —
built for low-amplitude CNV calling in genes like *EGFR*, *ERBB2*, *MET* from
circulating tumor DNA or tissue panels.

It covers the full pipeline end to end: a **panel-of-normals master-depth
estimate**, **CBS-based re-segmentation** of the raw panel BED, and
correlation-based **anchor recalibration**, so training operates on segments
that truly reflect coverage structure rather than user-supplied bins.

---

## Why CesaR

Panel CNV callers lean on a *panel of normals* (PoN) to model expected depth.
Two problems recur, and CesaR addresses both:

1. **Bins don't match biology.** A raw panel BED chops targets into arbitrary
   windows. CBS re-segmentation on the PoN master depth produces breaks that
   follow real coverage structure.
2. **Same-gene anchors mask amplifications.** When a whole gene is amplified,
   its segments become each other's most-correlated PoN partners — the anchor
   mean tracks the target and the signal collapses. CesaR excludes same-gene
   anchors **by default**.

The payoff is monotonic, well-separated dose-response even at low spike-in
levels — see the [worked demo](#quick-demo).

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("nishuai/CesaR", subdir = "pkg/CesaR")
```

Or from a local clone:

```r
remotes::install_local("pkg/CesaR")
# or, from the shell:  R CMD INSTALL pkg/CesaR
```

**Dependencies:** `PSCBS`, `MASS` (plus base `stats`, `utils`).

---

## The four-step pipeline

Each step is independent — you can swap in your own depth matrix or BED at any
point. Inputs are deliberately free of absolute paths.

```
 PoN pileups          raw panel BED
      │                     │
      ▼                     │
┌───────────────────┐       │
│ build_master_depth│───────┤   average per-base depth across normals
└───────────────────┘       │
      │ master_depth        │
      ▼                     ▼
┌───────────────────────────────┐
│ cesar_segment                 │   CBS re-segmentation → fine-grained BED
└───────────────────────────────┘
      │ segmented BED
      ▼
┌───────────────────┐
│ cesar_train        │   per-segment correlation anchors + ratio model (PoN)
└───────────────────┘
      │ cesar_model
      ▼
┌───────────────────┐
│ cesar_detect       │   per-sample CNV calls (copy_ratio + confidence)
└───────────────────┘
```

| Step | Function | What it does |
|------|----------|--------------|
| 1 | `build_master_depth(pon_dir, ...)` | Average per-base depth across a PoN. Accepts a directory of pileups, explicit file paths, or in-memory data frames. Returns `data.frame(CHR, POSITION, DEPTH)`. |
| 2 | `cesar_segment(master_depth, panel_bed)` | CBS re-segmentation (`PSCBS::segmentByCBS`). Parent region start/end are forced as anchors so sub-segments tile without gaps. Five tunable merge thresholds; regions < 600 bp pass through unchanged. |
| 3 | `cesar_train(depth_matrix, bed)` | Fits per-segment correlation-based anchors and a ratio distribution on the PoN. `exclude_same_gene = TRUE` by default. |
| 4 | `cesar_detect(model, depth)` | Applies the model to one test sample. Takes a numeric depth vector of length `nrow(bed)` — use `bed_depth_in_pileup()` to convert a pileup. |

**Helpers:** `read_pileup()`, `read_panel_bed()`, `bed_depth_in_pileup()`.

---

## Quick demo

CesaR ships with `cesar_demo`, a spike-in dilution series (28 samples ×
196 segments across 4 genes) so you can run the whole pipeline with no external
data.

```r
library(CesaR)
data("cesar_demo", package = "CesaR")

# Train on the 10 normal samples (same-gene exclusion is on by default)
normals <- cesar_demo$sample_meta$condition == "normal"
model <- cesar_train(
  depth_matrix   = cesar_demo$depth_matrix[normals, ],
  bed            = cesar_demo$bed,
  min_mean_depth = 50        # matrix is trimmed; use 50, not the 200 default
)

# Detect on any sample
res <- cesar_detect(model, cesar_demo$depth_matrix[1, ])
res$genes
```

Run the full cohort dose-response in one line:

```r
source(system.file("scripts/run_demo_full.R", package = "CesaR"))
```

**Expected output** (median copy_ratio / confidence across replicates):

```
          condition   n  MET_copy  MET_conf  ERBB2_copy  ERBB2_conf
             normal  10     1.000       0.4       1.003         0.2
 met2.1875ERBB2.625   6     1.091       5.0       1.243        26.0
   met2.375ERBB3.25    6     1.238      22.2       1.533        72.7
    met2.75ERBB4.75    6     1.383      33.8       2.245        98.0
```

Copy ratio and confidence rise monotonically with spike-in level, while normals
stay flat at ~1.0 with confidence near 0 — exactly what a well-calibrated caller
should produce.

---

## Real-world example: A269 NSCLC panel

`examples/A269_MET/` walks through *MET* amplification detection on a real
94-gene NSCLC panel:

```r
source("examples/A269_MET/run_A269_MET.R")
```

> **Note:** the raw pileup data under `examples/A269_MET/data/pileups/`
> (~135 MB of `.freq` files) is **not tracked in git**. See
> `examples/A269_MET/README.md` for how to supply your own pileups, or use the
> bundled `cesar_demo` dataset above to run end-to-end without external data.

More runnable scripts live in `demo/` (MET, ERBB2, and EGFR-19del scenarios,
both hand-rolled and package-based).

---

## Package layout

```
pkg/CesaR/            # the R package (install this)
├── R/                # build_master_depth, cesar_segment, cesar_train, cesar_detect, io
├── data/             # cesar_demo.rda  (lazy-loaded dataset)
├── inst/scripts/     # run_demo_full.R
├── man/              # generated docs
└── tests/            # testthat suite (57 assertions)
demo/                 # runnable demo scripts + figures/results
examples/A269_MET/    # real-panel worked example (raw data not tracked)
```

---

## Testing

```r
setwd("pkg/CesaR")
testthat::test_local(".")
```

57 `testthat` assertions across master-depth, segmentation, and
train/detect blocks. `R CMD check` passes with 0 errors / 0 warnings.

---

## Citation

Citation details to be added.

## License

[Artistic-2.0](https://opensource.org/licenses/Artistic-2.0)

