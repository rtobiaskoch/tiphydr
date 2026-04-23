test_that("is_aa_pattern detects amino acid patterns", {
  expect_true(is_aa_pattern("K"))
  expect_true(is_aa_pattern("M"))
  expect_true(is_aa_pattern("KM"))
  expect_true(is_aa_pattern("ACDEFGHIKLMNPQRSTVWY"))
})

test_that("is_aa_pattern detects nucleotide patterns", {
  expect_false(is_aa_pattern("ATG"))
  expect_false(is_aa_pattern("TAA"))
  expect_false(is_aa_pattern("ACGT"))
})

test_that("is_aa_pattern prioritizes amino acids for ambiguous patterns", {
  # A is both AA and nucleotide; should return TRUE (prioritize AA)
  expect_true(is_aa_pattern("A"))
})

test_that("is_aa_pattern is case insensitive", {
  expect_true(is_aa_pattern("k"))
  expect_true(is_aa_pattern("atg") == FALSE)
})

test_that("is_aa_pattern rejects invalid patterns", {
  expect_false(is_aa_pattern("XYZ123"))
  expect_false(is_aa_pattern(""))
})
