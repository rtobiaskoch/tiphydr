test_that("returns list with metadata and biostring elements", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  result <- suppressMessages(fasta_unnest(nested))
  expect_type(result, "list")
  expect_named(result, c("metadata", "biostring"))
  expect_s3_class(result$metadata, "data.frame")
  expect_s4_class(result$biostring, "DNAStringSet")
})

test_that("metadata does not contain sequence column", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  result <- suppressMessages(fasta_unnest(nested))
  expect_false("sequence" %in% names(result$metadata))
})

test_that("biostring length equals number of non-NULL rows", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  # All 4 metadata rows match — biostring should have 4 sequences
  result <- suppressMessages(fasta_unnest(nested))
  expect_equal(length(result$biostring), 4L)
})

test_that("NULL rows are dropped and messaged", {
  meta_extra <- tibble::add_row(make_test_metadata(),
                                 strain_id = "GHOST_999",
                                 date = as.Date("2021-01-01"),
                                 location = "Colorado", host = "Bird")
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(meta_extra, make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  expect_message(fasta_unnest(nested), "Dropping 1 rows")
  result <- suppressMessages(fasta_unnest(nested))
  expect_equal(nrow(result$metadata), 4L)
})

test_that("stops if sequence column is missing", {
  meta <- make_test_metadata()
  expect_error(fasta_unnest(meta), "must have a 'sequence' list-column")
})

test_that("verbose=FALSE suppresses drop message", {
  meta_extra <- tibble::add_row(make_test_metadata(),
                                 strain_id = "GHOST_999",
                                 date = as.Date("2021-01-01"),
                                 location = "Colorado", host = "Bird")
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(meta_extra, make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  expect_no_message(fasta_unnest(nested, verbose = FALSE))
})
