test_that("replaces IUPAC ambiguity codes and gaps with N, preserving class and names", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "ACGTRYN", s2 = "AC-GT.N"))
  result <- mask_ambiguous(seqs)

  expect_s4_class(result, "DNAStringSet")
  expect_equal(names(result), c("s1", "s2"))
  expect_equal(as.character(result), c(s1 = "ACGTNNN", s2 = "ACNGTNN"))
})

test_that("keep-set preserves gaps while still masking ambiguity codes", {
  seqs <- Biostrings::DNAStringSet(c(s1 = "AC-GTR"))
  # Keep gaps: the gap survives, the ambiguous R becomes N
  result <- mask_ambiguous(seqs, keep = "ACGT-")
  expect_equal(as.character(result), c(s1 = "AC-GTN"))
})

test_that("non-DNAStringSet input errors", {
  expect_error(mask_ambiguous("ACGT"), "DNAStringSet")
})
