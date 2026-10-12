aln <- Biostrings::DNAStringSet(c(s1 = "ACGT", s2 = "ACNT", s3 = "A-GA"))

test_that("counts and percents are correct at requested positions", {
  out <- nt_freq(aln, positions = c(2, 4))
  expect_equal(out$position, c(2, 4))
  expect_equal(out$C, c(2, 0))
  expect_equal(out$gap, c(1, 0))
  expect_equal(out$ambiguous, c(0, 0))
  expect_equal(out$pct_T, c(0, 66.67))
})

test_that("ambiguous codes are counted and excluded from percent denominator", {
  out <- nt_freq(aln, positions = 3)
  expect_equal(out$ambiguous, 1)
  expect_equal(out$n_acgt, 2)
  expect_equal(out$pct_G, 100)
})

test_that("DNAMultipleAlignment input matches DNAStringSet input", {
  msa <- Biostrings::DNAMultipleAlignment(as.character(aln))
  expect_equal(nt_freq(msa), nt_freq(aln))
})

test_that("default positions cover all columns; bad input errors", {
  expect_equal(nrow(nt_freq(aln)), 4)
  expect_error(nt_freq(aln, positions = 5), "within")
  expect_error(nt_freq("x"), "must be")
  expect_error(nt_freq(Biostrings::DNAStringSet(c("AC", "ACG"))), "same width")
})
