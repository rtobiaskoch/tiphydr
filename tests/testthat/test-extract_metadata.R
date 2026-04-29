# test-extract_metadata.R — unit tests for extract_metadata()
# Shared fixture: make_test_biostring() from helper-fixtures.R
# Names follow schema: VIRUS|YEAR|STATE|STRAIN_ID (e.g. "WNV|2021|CO|NY10_001")

test_that("errors when biostring is not an XStringSet", {
  expect_error(
    extract_metadata("not_a_biostring", names = c("virus", "year")),
    "biostring must be an XStringSet"
  )
})

test_that("errors when names is not a character vector", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = 1:3),
    "names must be a non-empty character vector"
  )
})

test_that("errors when names is empty", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = character(0)),
    "names must be a non-empty character vector"
  )
})

test_that("errors when names contains duplicates", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = c("virus", "virus", "year")),
    "names contains duplicate column names"
  )
})

test_that("errors when delim is not a single non-empty string", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = c("virus", "year"), delim = ""),
    "delim must be a single non-empty string"
  )
  expect_error(
    extract_metadata(bs, names = c("virus", "year"), delim = c("|", "/")),
    "delim must be a single non-empty string"
  )
})
