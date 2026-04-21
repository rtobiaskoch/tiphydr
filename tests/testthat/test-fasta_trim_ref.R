skip_if_no_mafft <- function() {
  testthat::skip_if(nchar(Sys.which("mafft")) == 0, "MAFFT not on PATH")
}

test_that("returns DNAStringSet with ref excluded", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()
  result <- suppressMessages(fasta_trim_ref(seqs, ref))
  expect_s4_class(result, "DNAStringSet")
  expect_false(names(ref) %in% names(result))
})

test_that("all trimmed sequences equal reference length", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()
  result <- suppressMessages(fasta_trim_ref(seqs, ref))
  # Reference is 84nt; all trimmed seqs must also be 84nt (gaps allowed)
  expect_true(all(Biostrings::width(result) == Biostrings::width(ref)))
})

test_that("sequences with leading insertions trim to reference coordinates", {
  skip_if_no_mafft()
  ref <- make_test_ref()
  # Add 5 extra nt before the reference region in one sequence
  padded_seq <- Biostrings::DNAStringSet(
    c(test_padded = paste0("AAAAA", as.character(ref[[1]])))
  )
  result <- suppressMessages(fasta_trim_ref(padded_seq, ref))
  # After trimming, the padded sequence should be the same length as reference
  # Use single-bracket [1] to keep a DNAStringSet (width() requires XStringSet, not XString)
  expect_equal(Biostrings::width(result[1]), Biostrings::width(ref))
})

test_that("ref of length > 1 throws error", {
  two_refs <- c(make_test_ref(), make_test_ref())
  expect_error(fasta_trim_ref(make_test_biostring(), two_refs), "length 1")
})

test_that("absent MAFFT gives informative stop message", {
  skip_if_no_mafft()
  with_mocked_bindings(
    `Sys.which` = function(x) if (x == "mafft") "" else base::Sys.which(x),
    expect_error(fasta_trim_ref(make_test_biostring(), make_test_ref()),
                 "MAFFT not found"),
    .package = "base"
  )
})

test_that("drop=TRUE removes extra sequences from alignment", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()

  # First align all 5 sequences
  full_aln <- suppressMessages(fasta_trim_ref(seqs, ref))

  # Provide subset of 2 sequences, use full alignment as starting point
  # With drop=TRUE, only the 2 sequences should appear in output
  subset_seqs <- seqs[1:2]
  result <- suppressMessages(
    fasta_trim_ref(subset_seqs, ref, alignment = full_aln, drop = TRUE)
  )
  expect_equal(length(result), 2L)
  expect_true(all(names(result) %in% names(subset_seqs)))
})
