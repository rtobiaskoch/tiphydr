test_that("empty accession throws error", {
  expect_error(fetch_fasta(character(0), "out.fasta"), "non-empty character vector")
})

test_that("non-character accession throws error", {
  expect_error(fetch_fasta(123, "out.fasta"), "non-empty character vector")
})

test_that("bad filename throws error", {
  expect_error(fetch_fasta("AF481864", ""), "non-empty string")
})

test_that("single accession fetches and writes fasta", {
  skip_if_offline()
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp))

  result <- suppressMessages(fetch_fasta("AF481864", tmp))

  expect_s4_class(result, "DNAStringSet")
  expect_gte(length(result), 1L)
  expect_true(file.exists(tmp))
})

test_that("multiple accessions write to single file", {
  skip_if_offline()
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp))

  result <- suppressMessages(fetch_fasta(c("AF481864", "NC_009942"), tmp))

  expect_s4_class(result, "DNAStringSet")
  expect_gte(length(result), 2L)
})
