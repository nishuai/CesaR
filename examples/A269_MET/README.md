# A269 MET amplification — CesaR worked example

A full-size, real-world example of the **CesaR** 4-step CNV pipeline
on a 9-sample circulating-tumor-DNA cohort. The original CNVkit run
missed all 3 MET amplifications; CesaR recovers them with confidence
scores between 14 and 100.

## Cohort

| Sample | Role | MET status |
|---|---|---|
| A269010002 | test | abnormality |
| A269010003 | test | abnormality |
| A269010004 | test | abnormality |
| A269010005 | PoN | negative |
| A269010006 | PoN | negative |
| A269010007 | PoN | negative |
| A269010009 | PoN | negative |
| A269010010 | PoN | negative |
| A269010011 | PoN | negative |

Per-investigator annotation. The 6 MET-negative samples form a
panel-of-normals; the 3 MET-positive samples are scored against it.

## Files

```
examples/A269_MET/
├── README.md
├── run_A269_MET.R                                 4-step pipeline
└── data/
    ├── segmented_bed_A269_hg38_annotated.bed      16-gene hg38 panel
    └── pileups/
        ├── A269010002F01.sort.bam.mpileup.freq    ~15 MB each
        ├── ...                                    9 files, ~135 MB total
        └── A269010011F01.sort.bam.mpileup.freq
```

Each `.freq` file has columns
`CHR  POSITION  DEPTH  REF  R+  R-  A+  A-  C+  C-  T+  T-  G+  g-`.
CesaR only uses the first three (CHR/POSITION/DEPTH) via
`read_pileup()`; the allele-count columns are preserved but ignored.

The BED is a 1645-row hg38-coordinate panel covering 16 genes
(MET, EGFR, ERBB2, ...). The reference `.bed` shipped here is the
already-segmented version; `run_A269_MET.R` collapses it back to a
"raw" 1355-region panel and re-segments it with CesaR's CBS path so
the full 4-step flow is exercised.

> Note: the vendor also distributes an hg19-coordinate panel BED for
> this assay, but its coordinates do not match these hg38 pileups, so
> this example uses the hg38 BED shipped above.

## How to run

From the repository root:

```bash
Rscript examples/A269_MET/run_A269_MET.R
```

Prerequisites:

```r
devtools::install("pkg/CesaR/")
```

(Or `install.packages("PSCBS")` + `R CMD INSTALL pkg/CesaR/` if you
prefer not to use `devtools`.)

## Pipeline walkthrough

The script applies CesaR's four public functions in sequence:

| Step | Function | Input | Output |
|---|---|---|---|
| 1 | `build_master_depth()` | 6 MET-negative pileups | 345,060-position depth profile |
| 2 | `cesar_segment()` | master depth + raw panel | 1363 fine-grained segments |
| 3 | `cesar_train()` | 6 × 1363 coverage matrix + segmented BED | `cesar_model` (anchors + ratio params) |
| 4 | `cesar_detect()` | model + per-sample depth vector | per-segment CNV calls aggregated to genes |

`exclude_same_gene = TRUE` in step 3 — anchors for any MET segment
are drawn from non-MET segments only, preventing whole-gene MET
amplification from collapsing the anchor-ratio back to 1.

## Expected output

```
=== MET CesaR summary (positive samples) ===
     sample is_MET_pos n_seg cn_abs_median cn_abs_mean cn_abs_max cn_abs_min
 A269010002       TRUE    19      3.228299    3.299398   3.928884   2.813364
 A269010003       TRUE    19      2.115418    2.120455   2.409768   1.756411
 A269010004       TRUE    19      1.850347    1.817479   2.074783   1.483069
     pval_min      conf
 1.00000e-100 100.00000
  1.32081e-14  13.87916
  3.23922e-33  32.48956
```

`cn_abs_*` are diploid-baseline absolute copy numbers (median /
mean / max / min across the 19 MET segments). `conf` is
`-log10(pval_min)`. All three MET-positive samples are called with
confidence ≥ 14; the strongest (002) saturates at 100 (the script's
`min_p = 1e-100` floor).

## Why this is a good test case

* **Real ctDNA data, not simulated.** Coverage is uneven and
  GC-correlated — the conditions that broke CNVkit and motivated
  CesaR's anchor recalibration.
* **Whole-gene amplification.** MET amplification spans all 19
  segments of the gene, which is exactly the case where
  `exclude_same_gene = TRUE` matters. With it off, the recalibrated
  ratio collapses toward 1 and the call disappears.
* **Modest-size PoN (n = 6).** Demonstrates that the pipeline runs
  with a small, clean PoN; the package's hard floor is 5.

## Versus the in-package demo

`pkg/CesaR/inst/scripts/run_demo.R` is the minimal in-package demo:
it ships with 9 trimmed PoN pileups, a pre-segmented 196-region
panel, and one positive + one negative test sample. It's the right
fixture for unit-style validation.

This `examples/A269_MET/` directory is the full-size companion: real
ctDNA data, a real diagnostic panel, and the precise cohort that
motivated the package. Use it when you want to see CesaR work on
inputs of production scale.
