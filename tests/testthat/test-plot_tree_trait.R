test_that("plot_tree_trait returns a ggplot object for a discrete trait", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  result  <- plot_tree_trait(tree_df, "inferred_state")
  expect_s3_class(result, "ggplot")
})

test_that("plot_tree_trait handles a continuous (numeric) trait", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  # confidence_state is numeric — should now render as a continuous gradient,
  # not error. Pass a non-existent confidence column so the size/alpha
  # encoding is skipped (and does not collide with the trait itself).
  result <- plot_tree_trait(tree_df, "confidence_state", confidence = "none")
  expect_s3_class(result, "ggplot")
})

test_that("plot_tree_trait drops confidence encoding when column is absent", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  tree_df$confidence_state <- NULL
  result  <- plot_tree_trait(tree_df, "inferred_state")
  expect_s3_class(result, "ggplot")
})

test_that("plot_tree_trait errors if tree_df lacks node/parent columns", {
  bad_df <- tibble::tibble(x = 1:3, y = letters[1:3])
  expect_error(
    plot_tree_trait(bad_df, "y"),
    regexp = "node, parent"
  )
})

test_that("plot_tree_trait errors if trait column is absent", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  expect_error(
    plot_tree_trait(tree_df, "nonexistent_column"),
    regexp = "not found in tree_df"
  )
})
