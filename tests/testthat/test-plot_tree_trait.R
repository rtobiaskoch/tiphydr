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

test_that("branch_state switches branch colouring between parent and node", {
  tree    <- make_intro_tree()
  tree_df <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))

  # Read the colour ggplot2 actually resolved for each branch out of the BUILT
  # layer data, keyed by node -- not the aes() spec. A branch_state that was
  # accepted but never wired through would still pass a spec-level check.
  branch_colours <- function(p) {
    d <- ggplot2::ggplot_build(p)$data[[1]]
    stats::setNames(as.character(d$colour), as.character(d$node))
  }

  p_parent <- plot_tree_trait(tree_df, "inferred_state", branch_state = "parent")
  p_node   <- plot_tree_trait(tree_df, "inferred_state", branch_state = "node")

  expect_s3_class(p_parent, "ggplot")
  expect_s3_class(p_node, "ggplot")

  # Find a branch whose parent state differs from its own -- that is what an
  # introduction IS, so the fixture must contain at least one. Without such a
  # node the two modes would agree trivially and the test would prove nothing.
  pd <- p_parent$data
  differing <- pd$node[!is.na(pd$parent_trait) & pd$parent_trait != pd$inferred_state]
  expect_gt(length(differing), 0)

  cp <- branch_colours(p_parent)
  cn <- branch_colours(p_node)
  n  <- as.character(differing[[1]])

  # The two modes must disagree on exactly the branches whose endpoints differ.
  expect_false(identical(cp[[n]], cn[[n]]))

  # And "node" mode must paint that branch the same colour as the branch whose
  # PARENT is that node -- i.e. the colour the state itself maps to.
  child <- pd$node[!is.na(pd$parent) & pd$parent == differing[[1]]]
  child <- child[pd$inferred_state[match(child, pd$node)] ==
                   pd$inferred_state[match(differing[[1]], pd$node)]]
  skip_if(length(child) == 0, "fixture has no same-state child to compare against")
  expect_identical(cn[[n]], cp[[as.character(child[[1]])]])
})
