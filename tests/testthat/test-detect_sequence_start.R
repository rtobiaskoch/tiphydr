test_that("detect_sequence_start finds modal position when all at same position", {
  # All sequences have ATG at position 1
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAACCC",
    "Seq2" = "ATGGGGTTT",
    "Seq3" = "ATGCCCCAA"
  ))

  result <- detect_sequence_start(dna, pattern = "ATG", verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("position", "table"))
  expect_equal(result$position, 1L)
  expect_s3_class(result$table, "table")
})

test_that("detect_sequence_start returns mode when positions vary", {
  # Sequences have ATG at different positions; position 1 appears most (mode)
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAACCC",  # position 1
    "Seq2" = "ATGGGGTTT",  # position 1
    "Seq3" = "GGGATGCAA",  # position 4
    "Seq4" = "ATGCCCCAA"   # position 1
  ))

  result <- detect_sequence_start(dna, pattern = "ATG", verbose = FALSE)

  expect_equal(result$position, 1L)
})

test_that("detect_sequence_start excludes sequences missing pattern", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAACCC",  # has ATG at 1
    "Seq2" = "GGGGGGGG",   # no ATG
    "Seq3" = "ATGCCCCAA"   # has ATG at 1
  ))

  result <- detect_sequence_start(dna, pattern = "ATG", verbose = FALSE)

  expect_equal(result$position, 1L)
  # table should contain only 2 sequences (at position 1)
  expect_equal(sum(result$table), 2L)
})

test_that("detect_sequence_start messages when verbose=TRUE", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAACCC",
    "Seq2" = "ATGGGGTTT"
  ))

  expect_message(
    detect_sequence_start(dna, pattern = "ATG", verbose = TRUE),
    "Start codon 'ATG' detected at position"
  )
})

test_that("detect_sequence_start stops with error if no sequences contain pattern", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "GGGGGGGG",
    "Seq2" = "CCCCCCCC"
  ))

  expect_error(
    detect_sequence_start(dna, pattern = "ATG", verbose = FALSE),
    "No sequences contain pattern"
  )
})

test_that("detect_sequence_start works with different patterns", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "TAACCCTAG",  # TAA at position 1
    "Seq2" = "TAAGGGGG"    # TAA at position 1
  ))

  result <- detect_sequence_start(dna, pattern = "TAA", verbose = FALSE)

  expect_equal(result$position, 1L)
})

test_that("detect_sequence_start default pattern is ATG", {
  dna <- Biostrings::DNAStringSet(c(
    "Seq1" = "ATGAAACCC",
    "Seq2" = "ATGGGGTTT"
  ))

  result <- detect_sequence_start(dna, verbose = FALSE)

  expect_equal(result$position, 1L)
})
