test_that("single path returns DNAStringSet directly (not a list)", {
  path   <- system.file("extdata", "wnv_ref.fasta", package = "tiphydr")
  result <- suppressMessages(fasta_read(path))
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), 1L)
})

test_that("multiple paths return named list of DNAStringSets", {
  ref_path <- system.file("extdata", "wnv_ref.fasta",  package = "tiphydr")
  seq_path <- system.file("extdata", "wnv_seqs.fasta", package = "tiphydr")
  result   <- suppressMessages(fasta_read(c(ref_path, seq_path)))
  expect_type(result, "list")
  expect_length(result, 2L)
  expect_named(result, c("wnv_ref", "wnv_seqs"))
  expect_s4_class(result[["wnv_ref"]],  "DNAStringSet")
  expect_s4_class(result[["wnv_seqs"]], "DNAStringSet")
  expect_equal(length(result[["wnv_seqs"]]), 5L)
})

test_that("non-existent file throws informative error", {
  expect_error(fasta_read("ghost.fasta"), "Files not found")
})

test_that("verbose=FALSE suppresses all messages", {
  path <- system.file("extdata", "wnv_ref.fasta", package = "tiphydr")
  expect_silent(fasta_read(path, verbose = FALSE))
})

test_that("duplicate paths warn about duplicate list names", {
  path   <- system.file("extdata", "wnv_ref.fasta", package = "tiphydr")
  expect_warning(
    suppressMessages(fasta_read(c(path, path))),
    "Duplicate list names"
  )
})
