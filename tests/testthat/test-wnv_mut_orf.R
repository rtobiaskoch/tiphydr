# Exercises the shipped wnv_mut_orf lookup table's cross-family combination
# rows (data-raw/wnv_mutations.R) against define_lineage(). Builds a minimal
# synthetic AAStringSet with residues set only at the diagnostic ORF
# positions used by wnv_mut_orf (449, 2209, 2842, 1331, 2513); background
# positions are filled with a residue ("S") that never appears in the table,
# so it never coincidentally matches a diagnostic position.

wnv_test_seq <- function(env = "A", ns4a = FALSE, ns5 = FALSE, ns2a = FALSE, ns4b = FALSE) {
  seq_len <- 2842L
  chars <- rep("S", seq_len)
  chars[1] <- "M" # start codon so aa_start = 1, no fuzzy-start warning
  chars[449] <- env
  chars[2209] <- if (ns4a) "T" else "S"
  chars[2842] <- if (ns5) "R" else "S"
  chars[1331] <- if (ns2a) "K" else "S"
  chars[2513] <- if (ns4b) "M" else "S"
  Biostrings::AAStringSet(setNames(paste(chars, collapse = ""), "test_seq"))
}

test_that("cross-family combination: NY10-SW03_NS5", {
  # All of NY10's markers (NS2A, NS4B) plus SW03's NS5 marker only (no NS4A)
  seq <- wnv_test_seq(env = "A", ns5 = TRUE, ns2a = TRUE, ns4b = TRUE)
  result <- suppressMessages(define_lineage(seq, muts = wnv_mut_orf, mut_type = "aa"))
  expect_equal(result$lineage, "NY10-SW03_NS5")
})

test_that("cross-family combination: NY10-SW03 (all five markers)", {
  seq <- wnv_test_seq(env = "A", ns4a = TRUE, ns5 = TRUE, ns2a = TRUE, ns4b = TRUE)
  result <- suppressMessages(define_lineage(seq, muts = wnv_mut_orf, mut_type = "aa"))
  expect_equal(result$lineage, "NY10-SW03")
})

test_that("plain NY10 does not trip a cross-family row", {
  seq <- wnv_test_seq(env = "A", ns2a = TRUE, ns4b = TRUE)
  result <- suppressMessages(define_lineage(seq, muts = wnv_mut_orf, mut_type = "aa"))
  expect_equal(result$lineage, "NY10")
})

test_that("plain SW03 does not trip a cross-family row", {
  seq <- wnv_test_seq(env = "A", ns4a = TRUE, ns5 = TRUE)
  result <- suppressMessages(define_lineage(seq, muts = wnv_mut_orf, mut_type = "aa"))
  expect_equal(result$lineage, "SW03")
})
