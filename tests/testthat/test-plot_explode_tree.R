# test-plot_explode_tree.R

test_that("plot_explode_tree returns a ggplot for valid input", {
  tree   <- make_intro_tree()
  result <- explode_tree(tree, make_intro_node_probs(tree))
  p      <- plot_explode_tree(result$introductions)
  expect_s3_class(p, "ggplot")
})

test_that("plot_explode_tree errors on missing required columns", {
  bad_df <- tibble::tibble(intro_clade_id = 1L, x = "a")
  expect_error(
    plot_explode_tree(bad_df),
    regexp = "missing required columns"
  )
})

test_that("plot_explode_tree produces one facet per destination deme", {
  tree   <- make_intro_tree()
  result <- explode_tree(tree, make_intro_node_probs(tree))
  p      <- plot_explode_tree(result$introductions)
  built  <- ggplot2::ggplot_build(p)
  n_demes <- dplyr::n_distinct(result$introductions$deme)
  expect_equal(length(built$layout$panel_params), n_demes)
})
