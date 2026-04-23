test_that("search_pattern_in_xstring finds nucleotide matches in order", {
  dna <- Biostrings::DNAStringSet(c("Seq1" = "ATGATGATG"))

  result <- search_pattern_in_xstring(dna[[1]], "ATG")

  expect_type(result, "integer")
  expect_equal(result, c(1L, 4L, 7L))
})

test_that("search_pattern_in_xstring returns empty for no matches", {
  dna <- Biostrings::DNAStringSet(c("Seq1" = "CGCGCGCG"))

  result <- search_pattern_in_xstring(dna[[1]], "ATG")

  expect_type(result, "integer")
  expect_equal(length(result), 0L)
})

test_that("search_pattern_in_xstring finds multiple matches", {
  dna <- Biostrings::DNAStringSet(c("Seq1" = "GAATTATTAG"))

  result <- search_pattern_in_xstring(dna[[1]], "AT")

  expect_equal(result, c(3L, 6L))
})

test_that("search_pattern_in_xstring works with amino acids", {
  aa <- Biostrings::AAStringSet(c("Seq1" = "MKVLKFL"))

  result <- search_pattern_in_xstring(aa[[1]], "K")

  expect_equal(result, c(2L, 5L))
})

test_that("search_pattern_in_xstring is case insensitive for amino acids", {
  aa <- Biostrings::AAStringSet(c("Seq1" = "MKVLKFL"))

  result <- search_pattern_in_xstring(aa[[1]], "k")

  expect_equal(result, c(2L, 5L))
})

test_that("search_pattern_in_xstring handles lowercase DNA patterns", {
  dna <- Biostrings::DNAStringSet(c("Seq1" = "ATGATGATG"))

  result <- search_pattern_in_xstring(dna[[1]], "atg")

  expect_equal(result, c(1L, 4L, 7L))
})

test_that("search_pattern_in_xstring single match", {
  dna <- Biostrings::DNAStringSet(c("Seq1" = "GATATATG"))

  result <- search_pattern_in_xstring(dna[[1]], "ATG")

  expect_equal(result, 6L)
})
