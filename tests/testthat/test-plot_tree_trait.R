test_that("plot_tree_trait returns a ggplot object for valid input", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_tip_dates(tree))
  result  <- plot_tree_trait(tree_df, "inferred_state")
  expect_s3_class(result, "ggplot")
})

test_that("plot_tree_trait errors if tree_df is not from build_tree_df()", {
  bad_df <- tibble::tibble(x = 1:3, y = letters[1:3])
  expect_error(
    plot_tree_trait(bad_df, "y"),
    regexp = "build_tree_df"
  )
})

test_that("plot_tree_trait errors if trait column is absent", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_tip_dates(tree))
  expect_error(
    plot_tree_trait(tree_df, "nonexistent_column"),
    regexp = "not found in tree_df"
  )
})

test_that("plot_tree_trait errors if trait column is numeric", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_tip_dates(tree))
  # confidence_state is numeric — the discrete-only guard should fire
  expect_error(
    plot_tree_trait(tree_df, "confidence_state"),
    regexp = "numeric"
  )
})
