test_that("most_likely_deme picks the highest-probability deme per node", {
  np <- data.frame(
    node        = c(1L, 2L, 3L),
    regionA     = c(0.1, 0.8, 0.0),
    regionB     = c(0.7, 0.2, 0.4),
    regionC     = c(0.2, 0.0, 0.6),
    check.names = FALSE
  )

  res <- most_likely_deme(np)

  expect_equal(res$node, c(1L, 2L, 3L))
  expect_equal(res$inferred_state, c("regionB", "regionA", "regionC"))
  expect_equal(res$confidence_state, c(0.7, 0.8, 0.6))
})

test_that("most_likely_deme de-duplicates identical rows (e.g. doubled tip rows)", {
  np <- data.frame(
    node        = c(1L, 1L, 2L),
    regionA     = c(1, 1, 0),
    regionB     = c(0, 0, 1),
    check.names = FALSE
  )

  res <- most_likely_deme(np)

  expect_equal(nrow(res), 2L)
  expect_equal(res$confidence_state, c(1, 1))
})

test_that("most_likely_deme breaks ties by column order", {
  np <- data.frame(node = 1L, regionA = 0.5, regionB = 0.5, check.names = FALSE)
  expect_equal(most_likely_deme(np)$inferred_state, "regionA")
})

test_that("most_likely_deme errors without a node column or deme columns", {
  expect_error(most_likely_deme(data.frame(regionA = 1)), "node")
  expect_error(most_likely_deme(data.frame(node = 1L)), "deme")
})
