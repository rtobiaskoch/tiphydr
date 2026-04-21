test_that("combines a list of DNAStringSets into one", {
  ref  <- make_test_ref()
  seqs <- make_test_biostring()
  result <- suppressMessages(fasta_combine(list(ref, seqs)))
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), length(ref) + length(seqs))
})

test_that("warns on duplicate sequence names", {
  seqs <- make_test_biostring()
  expect_warning(
    suppressMessages(fasta_combine(list(seqs, seqs))),
    "Duplicate sequence names"
  )
})

test_that("stops on empty list", {
  expect_error(fasta_combine(list()), "non-empty list")
})

test_that("verbose=FALSE suppresses messages", {
  expect_silent(fasta_combine(list(make_test_biostring()), verbose = FALSE))
})

test_that("verbose=FALSE still emits warning for duplicate names", {
  seqs <- make_test_biostring()
  expect_warning(
    fasta_combine(list(seqs, seqs), verbose = FALSE),
    "Duplicate sequence names"
  )
})

test_that("single-element list combines without warning", {
  seqs <- make_test_biostring()
  result <- suppressMessages(fasta_combine(list(seqs)))
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), length(seqs))
})
