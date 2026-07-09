test_that("cesar_segment: short regions are kept unsegmented", {
  ## A 100-bp region < min_region_length=600 should pass through as 1 segment.
  master <- data.frame(CHR = "chr1", POSITION = 100:200,
                       DEPTH = rep(500, 101), stringsAsFactors = FALSE)
  panel <- data.frame(chr = "chr1", start = 100L, end = 200L,
                      gene = "EGFR", stringsAsFactors = FALSE)
  seg <- cesar_segment(master, panel, verbose = FALSE)
  expect_equal(nrow(seg), 1L)
  expect_equal(seg$gene, "EGFR")
  expect_equal(seg$start, 100L)
  expect_equal(seg$end, 200L)
})

test_that("cesar_segment: long uniform-depth regions yield one segment", {
  ## A 1000-bp region with constant depth should not produce extra breakpoints.
  master <- data.frame(CHR = "chr1", POSITION = 1000:2000,
                       DEPTH = rep(800, 1001), stringsAsFactors = FALSE)
  panel <- data.frame(chr = "chr1", start = 1000L, end = 2000L,
                      gene = "MET", stringsAsFactors = FALSE)
  seg <- cesar_segment(master, panel, verbose = FALSE)
  expect_equal(nrow(seg), 1L)
  expect_equal(seg$gene, "MET")
})

test_that("cesar_segment: a clear depth step is detected as a changepoint", {
  ## 1000-bp region, half low (200x) and half high (800x).
  ## CBS should find the change at position ~1500.
  pos   <- 1000:2000
  depth <- ifelse(pos < 1500, 200, 800)
  master <- data.frame(CHR = "chr1", POSITION = pos, DEPTH = depth,
                       stringsAsFactors = FALSE)
  panel <- data.frame(chr = "chr1", start = 1000L, end = 2000L,
                      gene = "ERBB2", stringsAsFactors = FALSE)
  seg <- cesar_segment(master, panel, verbose = FALSE)

  expect_gte(nrow(seg), 2L)  # the step should produce at least 2 sub-segments
  expect_true(all(seg$gene == "ERBB2"))
  ## The break should land near 1500, within 50 bp tolerance
  break_pos <- seg$end[seg$end != 2000L]
  expect_true(any(abs(break_pos - 1500) < 50))
})

test_that("cesar_segment: gene labels propagate to all sub-segments", {
  pos   <- 1000:2000
  depth <- ifelse(pos < 1500, 200, 800)
  master <- data.frame(CHR = "chr1", POSITION = pos, DEPTH = depth,
                       stringsAsFactors = FALSE)
  panel <- data.frame(chr = "chr1", start = 1000L, end = 2000L,
                      gene = "MyGene", stringsAsFactors = FALSE)
  seg <- cesar_segment(master, panel, verbose = FALSE)
  expect_true(all(seg$gene == "MyGene"))
})

test_that("cesar_segment: regions with no master depth coverage are kept", {
  master <- data.frame(CHR = "chr1", POSITION = 100:200, DEPTH = 500,
                       stringsAsFactors = FALSE)
  ## Panel region on chr2, no master coverage → should still appear in output.
  panel <- data.frame(chr = c("chr1", "chr2"),
                      start = c(100L, 5000L),
                      end   = c(200L, 6000L),
                      gene  = c("A", "B"),
                      stringsAsFactors = FALSE)
  seg <- cesar_segment(master, panel, verbose = FALSE)
  ## chr1 region: kept as 1 segment (short region < 600 bp)
  ## chr2 region: kept as 1 segment (no master depth)
  expect_equal(nrow(seg), 2L)
  expect_setequal(seg$gene, c("A", "B"))
})

test_that("cesar_segment: rejects malformed master_depth", {
  bad <- data.frame(POS = 1:10, DEPTH = rep(100, 10))
  panel <- data.frame(chr = "chr1", start = 1L, end = 100L, gene = "A",
                      stringsAsFactors = FALSE)
  expect_error(
    cesar_segment(bad, panel, verbose = FALSE),
    "CHR/POSITION/DEPTH"
  )
})

test_that("cesar_segment: total segment count is at least nrow(panel)", {
  ## Mix short and long, with a step in one of the long ones.
  pos   <- 1:3000
  depth <- ifelse(pos < 1500, 200, 800)
  master <- data.frame(CHR = "chr1", POSITION = pos, DEPTH = depth,
                       stringsAsFactors = FALSE)
  panel <- data.frame(
    chr   = c("chr1", "chr1", "chr1"),
    start = c(1L,    1000L,  2500L),
    end   = c(500L,  2000L,  2900L),  # row 1 short, row 2 long with step, row 3 short
    gene  = c("A",   "B",    "C"),
    stringsAsFactors = FALSE
  )
  seg <- cesar_segment(master, panel, verbose = FALSE)
  expect_gte(nrow(seg), 3L)
  ## Genes preserved per parent region
  expect_true(all(seg$gene[seg$start < 500]   == "A"))
  expect_true(all(seg$gene[seg$start >= 1000 & seg$end <= 2000] == "B"))
  expect_true(all(seg$gene[seg$start >= 2500] == "C"))
})
