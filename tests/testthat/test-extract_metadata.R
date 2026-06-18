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

test_that("returns a tibble with correct dimensions for exact field count", {
  # make_test_biostring() names: "WNV|YEAR|STATE|STRAIN_ID" (4 pipe-fields)
  bs     <- make_test_biostring()
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 6L)   # 6 sequences in fixture
  expect_equal(ncol(result), 5L)   # taxa + 4 fields
  expect_named(result, c("taxa", "virus", "year", "state", "strain_id"))
})

test_that("extracts correct values from pipe-delimited headers", {
  bs     <- make_test_biostring()
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  # First sequence: "WNV|2021|CO|NY10_001"
  expect_equal(result$taxa[1],      "WNV|2021|CO|NY10_001")
  expect_equal(result$virus[1],     "WNV")
  expect_equal(result$year[1],      "2021")
  expect_equal(result$state[1],     "CO")
  expect_equal(result$strain_id[1], "NY10_001")
})

test_that("extra header fields beyond length(names) are silently dropped", {
  bs <- make_test_biostring()
  # Request only 2 columns — 3rd and 4th fields are dropped with no warning
  expect_no_warning(
    result <- extract_metadata(bs, names = c("virus", "year"), delim = "|")
  )

  expect_equal(ncol(result), 3L)   # taxa + 2 fields
  expect_named(result, c("taxa", "virus", "year"))
  expect_equal(result$virus[1], "WNV")
  expect_equal(result$year[1],  "2021")
})

test_that("works with a single sequence", {
  bs     <- make_test_biostring()[1]   # subset to 1 sequence
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  expect_equal(nrow(result), 1L)
  expect_equal(result$strain_id, "NY10_001")
})

test_that("works with a custom delimiter", {
  # Build a DNAStringSet with slash-delimited names
  bs <- Biostrings::DNAStringSet(c(
    "WNV/2021/CO" = "ATCG",
    "WNV/2020/WY" = "GCTA"
  ))
  result <- extract_metadata(bs, names = c("virus", "year", "state"), delim = "/")

  expect_equal(result$state, c("CO", "WY"))
})

test_that("pads with NA and warns when a header has fewer fields than names", {
  # Build a set where one sequence is missing the 4th pipe-field
  bs <- Biostrings::DNAStringSet(c(
    "WNV|2021|CO|NY10_001" = "ATCG",   # 4 fields — OK
    "WNV|2020|WY"          = "GCTA"    # 3 fields — short by 1
  ))

  expect_warning(
    result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|"),
    "fewer fields than `names`"
  )

  # Full row: all fields correct
  expect_equal(result$strain_id[1], "NY10_001")
  # Short row: present fields not corrupted, missing field is NA
  expect_equal(result$state[2], "WY")
  expect_true(is.na(result$strain_id[2]))
})

test_that("warning message lists the affected sequence names", {
  bs <- Biostrings::DNAStringSet(c(
    "WNV|2021|CO|NY10_001" = "ATCG",
    "WNV|2020|WY"          = "GCTA"
  ))

  expect_warning(
    extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|"),
    "WNV\\|2020\\|WY"   # the short header should appear in the warning text
  )
})

test_that("no warning when all headers have exactly the right number of fields", {
  bs <- make_test_biostring()   # all 6 sequences have 4 pipe-fields
  expect_no_warning(
    extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")
  )
})

test_that("no warning when headers have more fields than names", {
  bs <- make_test_biostring()   # 4 pipe-fields; request only 2
  expect_no_warning(
    extract_metadata(bs, names = c("virus", "year"), delim = "|")
  )
})

test_that("names = NULL extracts all fields and names them V1, V2, ...", {
  bs     <- make_test_biostring()   # headers: "WNV|YEAR|STATE|STRAIN_ID" (4 fields)
  result <- extract_metadata(bs, delim = "|")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 6L)
  expect_equal(ncol(result), 5L)   # taxa + 4 fields
  expect_named(result, c("taxa", "V1", "V2", "V3", "V4"))
  expect_equal(result$taxa[1], "WNV|2021|CO|NY10_001")
  expect_equal(result$V1[1],   "WNV")
  expect_equal(result$V4[1],   "NY10_001")
})

test_that("names = NULL produces no warning even when headers differ in field count", {
  bs <- Biostrings::DNAStringSet(c(
    "WNV|2021|CO|NY10_001" = "ATCG",
    "WNV|2020|WY"          = "GCTA"   # shorter header — pads to V4 = NA, no warning
  ))
  expect_no_warning(result <- extract_metadata(bs, delim = "|"))

  expect_equal(ncol(result), 5L)   # taxa + 4 fields (widest header)
  expect_equal(result$taxa[2], "WNV|2020|WY")
  expect_true(is.na(result$V4[2]))  # short header padded with NA
})
