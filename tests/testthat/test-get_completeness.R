test_that("completeness is the unambiguous-base fraction", {
  ss <- Biostrings::DNAStringSet(c(
    full = "ACGT",  # 4/4  = 1.0
    half = "ACNN",  # 2/4  = 0.5  (N counts as missing)
    gappy = "AC--", # 2/4  = 0.5  (gaps count as missing)
    none = "NNNN"   # 0/4  = 0.0
  ))
  out <- get_completeness(ss)
  expect_equal(out$completeness[out$seq == "full"], 1.0)
  expect_equal(out$completeness[out$seq == "half"], 0.5)
  expect_equal(out$completeness[out$seq == "gappy"], 0.5)
  expect_equal(out$completeness[out$seq == "none"], 0.0)
})

test_that("zero-length sequences return 0, not NaN", {
  ss <- Biostrings::DNAStringSet(c(empty = ""))
  out <- get_completeness(ss)
  expect_equal(out$completeness, 0)
})

test_that("output is a data.frame with seq and completeness columns", {
  ss <- Biostrings::DNAStringSet(c(a = "ACGT", b = "ACGN"))
  out <- get_completeness(ss)
  expect_s3_class(out, "data.frame")
  expect_equal(names(out), c("seq", "completeness"))
  expect_equal(out$seq, c("a", "b"))
})

test_that("errors on a non-DNAStringSet input", {
  expect_error(get_completeness("ACGT"), "DNAStringSet")
})
