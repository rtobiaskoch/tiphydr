test_that("returns dataframe with strain and lineage columns", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_s3_class(result, "data.frame")
  expect_named(result, c("strain", "lineage"))
  expect_equal(nrow(result), length(make_test_biostring()))
})

test_that("strict match: NY10 requires both mutations", {
  # NY10_001 has pos10=A and pos25=C → matches NY10 (more specific than PARTIAL_A)
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("partial match assigns less-specific lineage, not unknown", {
  # NS2A_002 has pos10=A but pos25=T → matches PARTIAL_A (1 mutation) but not NY10 (2 mutations)
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "PARTIAL_A")
})

test_that("no match assigns unknown", {
  # OTHER_004 has no mutations
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|OTHER_004"], "unknown")
  # NS4B_003 has pos25=C but not pos10=A → neither lineage fully matched
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "unknown")
})

test_that("multiple-match: assigns most-specific (most mutations)", {
  # NY10_001 matches both NY10 (2 muts) and PARTIAL_A (1 mut) → NY10 wins
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("muts missing required columns throws error", {
  bad_muts <- tibble::tibble(lineage = "NY10", position = 10L, base = "A")
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), bad_muts),
    "lineage.*pos.*residue"
  )
})

test_that("invalid type throws error", {
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts(),
                   type = "protein"),
    "type must be"
  )
})
