test_that("fasta_name_switch renames all sequences when every name matches", {
  fasta <- make_test_biostring()
  old_names <- names(fasta)
  new_names <- paste0("Sample_", seq_along(old_names))

  result <- fasta_name_switch(fasta, match_col = old_names, replace_col = new_names)

  expect_setequal(names(result), new_names)
  expect_equal(length(result), length(fasta))
})

test_that("fasta_name_switch matches on a substring of the fasta name", {
  fasta <- make_test_biostring()
  old_names <- names(fasta)

  # match by accession (last field) only, replace the whole header
  accessions <- vapply(strsplit(old_names, "\\|"), function(x) x[length(x)], character(1))
  new_names <- paste0("renamed|", accessions)

  result <- fasta_name_switch(fasta, match_col = accessions, replace_col = new_names)

  expect_setequal(names(result), new_names)
})

test_that("fasta_name_switch leaves unmatched names unchanged with a warning", {
  fasta <- make_test_biostring()
  old_names <- names(fasta)

  # only provide a lookup for the first sequence
  expect_warning(
    result <- fasta_name_switch(
      fasta,
      match_col   = old_names[1],
      replace_col = "Sample_1"
    ),
    regexp = "no match_col entry"
  )

  expect_equal(names(result)[1], "Sample_1")
  expect_equal(names(result)[-1], old_names[-1])
})

test_that("fasta_name_switch errors on duplicate match_col values", {
  fasta <- make_test_biostring()
  old_names <- names(fasta)

  expect_error(
    fasta_name_switch(
      fasta,
      match_col   = c(old_names[1], old_names[1]),
      replace_col = c("Sample_1", "Sample_1_dupe")
    ),
    regexp = "duplicate"
  )
})

test_that("fasta_name_switch errors when match_col and replace_col differ in length", {
  fasta <- make_test_biostring()

  expect_error(
    fasta_name_switch(fasta, match_col = names(fasta), replace_col = "Sample_1"),
    regexp = "length"
  )
})

test_that("fasta_name_switch errors on non-XStringSet input", {
  expect_error(
    fasta_name_switch(list(), match_col = "a", replace_col = "b"),
    regexp = "XStringSet"
  )
})
