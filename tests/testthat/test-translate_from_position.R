test_that("translate_from_position translates DNA from a start position", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA")
  )

  result <- translate_from_position(dna_seq, start_position = 1)

  expect_s4_class(result, "AAStringSet")
  expect_equal(length(result), 2)
  expect_equal(names(result), names(dna_seq))
  expect_equal(as.character(result)[[1]], "MK*")
  expect_equal(as.character(result)[[2]], "MP*")
})

test_that("translate_from_position respects start position offset", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG")
  )

  result <- translate_from_position(dna_seq, start_position = 2)

  expect_s4_class(result, "AAStringSet")
  # Position 2: "TGAAATAG" -> "TGA" (stop), then "AAT" -> N, then "AG" (incomplete codon ignored)
  # Result should be "*N"
  expect_equal(as.character(result)[[1]], "*N")
})

test_that("translate_from_position handles fuzzy codons", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGNNATAG")
  )

  result <- translate_from_position(dna_seq, start_position = 1)

  expect_match(as.character(result)[[1]], "X")
})

test_that("translate_from_position preserves sequence names", {
  dna_seq <- Biostrings::DNAStringSet(
    c("seq_A" = "ATGAAATAG", "seq_B" = "ATGCCCTAA")
  )

  result <- translate_from_position(dna_seq, start_position = 1)

  expect_equal(names(result), c("seq_A", "seq_B"))
})
