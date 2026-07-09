test_that("build_master_depth: in-memory pileups (fast path)", {
  ## Three samples with identical layout → fast rowMeans path.
  d1 <- data.frame(CHR = c("chr1","chr1","chr1"),
                   POSITION = c(100L, 101L, 102L),
                   DEPTH = c(100, 200, 300))
  d2 <- data.frame(CHR = c("chr1","chr1","chr1"),
                   POSITION = c(100L, 101L, 102L),
                   DEPTH = c(120, 180, 360))
  d3 <- data.frame(CHR = c("chr1","chr1","chr1"),
                   POSITION = c(100L, 101L, 102L),
                   DEPTH = c(80, 220, 240))
  m <- build_master_depth(pon_pileups = list(d1, d2, d3), verbose = FALSE)

  expect_s3_class(m, "data.frame")
  expect_equal(colnames(m), c("CHR", "POSITION", "DEPTH"))
  expect_equal(m$DEPTH, c(100, 200, 300))   # rowMeans of the three samples
  expect_equal(m$POSITION, c(100L, 101L, 102L))
})

test_that("build_master_depth: drops zero-depth positions when asked", {
  d1 <- data.frame(CHR = "chr1", POSITION = 1:4, DEPTH = c(0, 100, 0, 200))
  d2 <- data.frame(CHR = "chr1", POSITION = 1:4, DEPTH = c(0, 120, 0, 180))
  m <- build_master_depth(pon_pileups = list(d1, d2),
                          drop_zero_depth = TRUE, verbose = FALSE)
  expect_equal(m$POSITION, c(2L, 4L))
  expect_equal(m$DEPTH,    c(110, 190))

  m_keep <- build_master_depth(pon_pileups = list(d1, d2),
                               drop_zero_depth = FALSE, verbose = FALSE)
  expect_equal(nrow(m_keep), 4L)
})

test_that("build_master_depth: inner-join path when layouts differ", {
  ## d2 is missing position 101; only positions present in BOTH should survive.
  d1 <- data.frame(CHR = "chr1", POSITION = c(100L, 101L, 102L),
                   DEPTH = c(100, 200, 300))
  d2 <- data.frame(CHR = "chr1", POSITION = c(100L, 102L, 103L),
                   DEPTH = c(80, 240, 400))
  m <- build_master_depth(pon_pileups = list(d1, d2), verbose = FALSE)
  expect_equal(m$POSITION, c(100L, 102L))
  expect_equal(m$DEPTH,    c(90, 270))
})

test_that("build_master_depth: rejects ambiguous or empty inputs", {
  d1 <- data.frame(CHR = "chr1", POSITION = 1L, DEPTH = 100)
  expect_error(
    build_master_depth(pon_pileups = list(d1), verbose = FALSE),
    "at least 2"
  )
  expect_error(
    build_master_depth(pon_dir = "no_such_dir", verbose = FALSE),
    "not found"
  )
  expect_error(
    build_master_depth(pon_dir = ".",
                       pon_files = "x",
                       verbose = FALSE),
    "exactly one"
  )
})

test_that("build_master_depth: pon_files mode reads pileups from disk", {
  tmp <- withr::local_tempdir()
  d1 <- data.frame(CHR = c("chr1","chr1"),
                   POSITION = c(100L, 101L),
                   DEPTH = c(100, 200))
  d2 <- data.frame(CHR = c("chr1","chr1"),
                   POSITION = c(100L, 101L),
                   DEPTH = c(120, 220))
  utils::write.table(d1, file.path(tmp, "s1.depth"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(d2, file.path(tmp, "s2.depth"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  m <- build_master_depth(
    pon_files = c(file.path(tmp, "s1.depth"), file.path(tmp, "s2.depth")),
    verbose = FALSE
  )
  expect_equal(m$DEPTH, c(110, 210))
})

test_that("build_master_depth: pon_dir mode picks up files via pattern", {
  tmp <- withr::local_tempdir()
  d <- data.frame(CHR = "chr1", POSITION = 1L, DEPTH = 100)
  utils::write.table(d, file.path(tmp, "a.depth"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(d, file.path(tmp, "b.depth"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(d, file.path(tmp, "ignore.txt"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  m <- build_master_depth(pon_dir = tmp, pileup_pattern = "\\.depth$",
                          verbose = FALSE)
  expect_equal(nrow(m), 1L)
  expect_equal(m$DEPTH, 100)
})

test_that("build_master_depth: sorts output by chr then position", {
  d1 <- data.frame(CHR = c("chr2", "chr1", "chr1"),
                   POSITION = c(50L, 200L, 100L),
                   DEPTH = c(100, 200, 300))
  d2 <- d1; d2$DEPTH <- d2$DEPTH + 20
  m <- build_master_depth(pon_pileups = list(d1, d2), verbose = FALSE)
  ## chr1 positions first (100, 200), then chr2 position 50
  expect_equal(m$CHR, c("chr1", "chr1", "chr2"))
  expect_equal(m$POSITION, c(100L, 200L, 50L))
})
