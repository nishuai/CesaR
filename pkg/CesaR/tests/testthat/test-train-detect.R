test_that("cesar_train: rejects malformed inputs", {
  bed <- data.frame(chr = "chr1", start = 1L, end = 100L, gene = "A",
                    stringsAsFactors = FALSE)
  X <- matrix(100, nrow = 5, ncol = 1)
  expect_error(cesar_train(X, bed[, 1:3], verbose = FALSE),
               "chr/start/end/gene")
  expect_error(cesar_train(X, bed[rep(1, 2), ], verbose = FALSE),
               "must equal")  # ncol(X) != nrow(bed)
  expect_error(cesar_train(X[1:3, , drop = FALSE], bed, verbose = FALSE),
               "at least 5 PoN")
  expect_error(cesar_train(X, bed, min_anchors = 5, max_anchors = 5,
                           verbose = FALSE),
               "max_anchors")
})

test_that("cesar_train + cesar_detect: clean PoN and amplified test sample", {
  ## Build a synthetic 8-sample x 30-segment matrix:
  ##   - 25 'Other' segments draw from N(1000, 50)
  ##   - 5 'MET' segments draw from N(500, 25)
  ## Then for the test sample, MET segments are amplified 3x.
  set.seed(42)
  n_other <- 25; n_met <- 5
  bed <- data.frame(
    chr   = rep("chr7", n_other + n_met),
    start = seq_len(n_other + n_met) * 1000L,
    end   = seq_len(n_other + n_met) * 1000L + 500L,
    gene  = c(rep("Other", n_other), rep("MET", n_met)),
    stringsAsFactors = FALSE
  )
  ## PoN: 8 samples, each with a sample-specific scaling factor (1.0..1.5)
  pon <- matrix(0, nrow = 8L, ncol = nrow(bed))
  for (s in seq_len(8L)) {
    scale_s <- runif(1, 0.8, 1.4)
    pon[s, 1:n_other]      <- rnorm(n_other, 1000 * scale_s, 30)
    pon[s, (n_other+1):30] <- rnorm(n_met,    500 * scale_s, 15)
  }
  pon[pon < 0] <- 0

  model <- cesar_train(pon, bed, min_mean_depth = 50,
                       min_anchors = 3, max_anchors = 10,
                       exclude_same_gene = TRUE, verbose = FALSE)
  expect_s3_class(model, "cesar_model")
  expect_true(model$settings$exclude_same_gene)
  expect_equal(nrow(model$bed), 30L)

  ## Confirm same-gene anchor exclusion: no MET segment should anchor on MET
  for (j in (n_other + 1):30) {
    anc <- model$anchors[[j]]
    if (all(is.na(anc))) next
    expect_true(all(bed$gene[anc] != "MET"),
                info = sprintf("seg %d (MET) has same-gene anchor", j))
  }

  ## Build a "test" sample at PoN scale 1.0 with MET amplified 3x
  test_depth <- numeric(30)
  test_depth[1:n_other]      <- rnorm(n_other, 1000, 30)
  test_depth[(n_other+1):30] <- rnorm(n_met, 500 * 3, 30)

  res <- cesar_detect(model, test_depth, sample_id = "test_amp")
  expect_s3_class(res, "cesar_result")
  expect_equal(nrow(res$segments), 30L)

  ## MET should be the top hit, with copy_ratio > 2 and high confidence
  top <- res$genes[1L, ]
  expect_equal(top$gene, "MET")
  expect_gt(top$copy_ratio, 2)
  expect_gt(top$confidence, 5)
})

test_that("cesar_detect: held-out negative gives low confidence", {
  set.seed(7)
  bed <- data.frame(chr = "chr1", start = (1:20) * 1000L,
                    end = (1:20) * 1000L + 500L,
                    gene = c(rep("Other", 15), rep("MET", 5)),
                    stringsAsFactors = FALSE)
  pon <- matrix(0, nrow = 8L, ncol = 20)
  for (s in 1:8) {
    scale_s <- runif(1, 0.9, 1.1)
    pon[s, 1:15]  <- rnorm(15, 1000 * scale_s, 30)
    pon[s, 16:20] <- rnorm(5,   800 * scale_s, 25)
  }
  model <- cesar_train(pon, bed, min_mean_depth = 50,
                       min_anchors = 3, max_anchors = 8, verbose = FALSE)

  ## Generate a *new* normal at scale 1.0
  test_depth <- numeric(20)
  test_depth[1:15]  <- rnorm(15, 1000, 30)
  test_depth[16:20] <- rnorm(5,   800, 25)

  res <- cesar_detect(model, test_depth, sample_id = "neg")
  expect_lt(max(res$genes$confidence, na.rm = TRUE), 5)
})

test_that("cesar_detect: rejects bad inputs", {
  bed <- data.frame(chr = "chr1", start = 1:6 * 100L, end = 1:6 * 100L + 50L,
                    gene = "A", stringsAsFactors = FALSE)
  pon <- matrix(rnorm(48, 1000, 30), nrow = 8, ncol = 6)
  model <- cesar_train(pon, bed, min_mean_depth = 50,
                       min_anchors = 2, max_anchors = 4, verbose = FALSE)
  expect_error(cesar_detect("not a model", rep(100, 6)),
               "cesar_model")
  expect_error(cesar_detect(model, rep(100, 5)),
               "must equal nrow")
  expect_error(cesar_detect(model, matrix(100, 6, 1)),
               "numeric vector")
})

test_that("cesar_train: exclude_same_gene FALSE retains backward-compatible behavior", {
  set.seed(1)
  bed <- data.frame(chr = "chr1", start = (1:20) * 1000L,
                    end = (1:20) * 1000L + 500L,
                    gene = c(rep("Other", 15), rep("MET", 5)),
                    stringsAsFactors = FALSE)
  pon <- matrix(rnorm(8 * 20, 1000, 50), nrow = 8, ncol = 20)
  model_off <- cesar_train(pon, bed, min_mean_depth = 50,
                           min_anchors = 3, max_anchors = 8,
                           exclude_same_gene = FALSE, verbose = FALSE)
  expect_false(model_off$settings$exclude_same_gene)

  ## When exclude is OFF, MET targets MAY have MET anchors (no constraint).
  ## Just check the field is recorded; do not assert specific anchor identity.
  expect_equal(length(model_off$anchors), 20L)
})
