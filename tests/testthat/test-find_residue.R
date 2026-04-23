test_that("find_residue finds DNA pattern in DNA sequences", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGATGATG",
    "Seq2" = "ATGCCCTAG"
  ))

  result <- find_residue(dna, pattern = "ATG", verbose = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("sequence_name", "location", "pattern"))
  # Seq1: ATG at 1, 4, 7
  # Seq2: ATG at 1
  expect_equal(nrow(result), 4L)
  expect_equal(result$sequence_name, c("Seq1", "Seq1", "Seq1", "Seq2"))
  expect_equal(result$location, c(1L, 4L, 7L, 1L))
  expect_equal(result$pattern, c("ATG", "ATG", "ATG", "ATG"))
})

test_that("find_residue finds protein pattern in protein sequences", {
  aa <- Biostrings::AAStringSet(c(
    "Seq1" = "MKVLKFL",
    "Seq2" = "MKSLKFL"
  ))

  result <- find_residue(aa, pattern = "K", verbose = FALSE)

  expect_s3_class(result, "tbl_df")
  # Seq1: K at positions 2 and 5
  # Seq2: K at positions 2 and 5
  expect_equal(nrow(result), 4L)
  expect_equal(result$sequence_name, c("Seq1", "Seq1", "Seq2", "Seq2"))
})

test_that("find_residue auto-translates DNA with AA pattern", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAATAG"  # ATG (M) AAA (K) TAG (*)
  ))

  result <- find_residue(dna, pattern = "K", verbose = FALSE)

  expect_s3_class(result, "tbl_df")
  # Should translate and find K at position 2
  expect_equal(nrow(result), 1L)
  expect_equal(result$sequence_name, "Seq1")
  expect_equal(result$location, 2L)
})

test_that("find_residue handles multiple patterns", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAATAG",  # ATG, TAA, TAG
    "Seq2" = "ATGCCCTAA"   # ATG, TAA
  ))

  result <- find_residue(dna, pattern = c("ATG", "TAA"), verbose = FALSE)

  expect_s3_class(result, "tbl_df")
  # Long format: one row per match per pattern
  expect_true(all(result$pattern %in% c("ATG", "TAA")))
})

test_that("find_residue allows start_position override", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAATAG"
  ))

  result <- find_residue(
    dna,
    pattern = "K",
    start_position = 1,
    verbose = FALSE
  )

  # Should translate from position 1
  expect_s3_class(result, "tbl_df")
})

test_that("find_residue preserves sequence names", {
  dna <- Biostrings::DNAStringSet(c(
    "SeqA" = "ATGATGATG",
    "SeqB" = "ATGCCCTAG"
  ))

  result <- find_residue(dna, pattern = "ATG", verbose = FALSE)

  expect_true(all(result$sequence_name %in% c("SeqA", "SeqB")))
})

test_that("find_residue returns empty tibble when no matches", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "GCGCGCGC"
  ))

  result <- find_residue(dna, pattern = "ATG", verbose = FALSE)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0L)
  expect_named(result, c("sequence_name", "location", "pattern"))
})

test_that("find_residue returns long format (one row per match)", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGATGATG"
  ))

  result <- find_residue(dna, pattern = "ATG", verbose = FALSE)

  # Should have 3 rows (one per ATG match)
  expect_equal(nrow(result), 3L)
  # Each row represents one match
  expect_equal(result$pattern, c("ATG", "ATG", "ATG"))
})
