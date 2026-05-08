# Regenerate tests/testthat/testdata/ from fixture definitions in helper-fixtures.R.
# Run once from the project root whenever fixture data changes:
#
#   source("tests/write_fixtures.R")
#
# Do NOT call this inside test files — write_test_fixtures() must never run
# automatically during devtools::test().

library(Biostrings)
library(tibble)
source("tests/testthat/helper-fixtures.R")
write_test_fixtures()
