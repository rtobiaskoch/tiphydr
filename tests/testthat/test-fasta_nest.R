test_that("adds sequence list-column to metadata", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|")
  ))
  expect_true("sequence" %in% names(result))
  expect_equal(nrow(result), nrow(meta))
  expect_s4_class(result$sequence[[1]], "DNAStringSet")
  expect_equal(length(result$sequence[[1]]), 1L)
})

test_that("matched sequence name matches id_col value via delimiter split", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|")
  ))
  # NY10_001 should match WNV|2021|CO|NY10_001
  expect_equal(names(result$sequence[[1]]), "WNV|2021|CO|NY10_001")
})

test_that("unmatched row stores NULL and raises warning", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  meta_extra <- tibble::add_row(meta, strain_id = "GHOST_999",
                                 date = as.Date("2021-01-01"),
                                 location = "Colorado", host = "Bird")
  expect_warning(
    suppressMessages(fasta_nest(meta_extra, seqs, id_col = "strain_id", delim = "|")),
    regexp = "no matching sequence"
  )
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta_extra, seqs, id_col = "strain_id", delim = "|")
  ))
  expect_null(result$sequence[[nrow(result)]])
})

test_that("wrong id_col name throws informative error", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  expect_error(
    fasta_nest(meta, seqs, id_col = "bad_col", delim = "|"),
    "not found in df columns"
  )
})

test_that("verbose=FALSE suppresses messages", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  expect_no_message(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|", verbose = FALSE)
  ))
})
