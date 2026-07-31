test_that("build_tree_df annotates every node with state, confidence and date", {
  tree <- make_intro_tree()
  res  <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))

  n_nodes <- ape::Ntip(tree) + tree$Nnode
  expect_equal(nrow(res), n_nodes)
  expect_true(all(c("inferred_state", "confidence_state",
                    "inferred_date", "is_tip") %in% names(res)))
  expect_false(anyNA(res$inferred_state))
  expect_equal(sum(res$is_tip), ape::Ntip(tree))
})

test_that("build_tree_df keeps the tbl_tree class for tidytree traversal", {
  tree <- make_intro_tree()
  res  <- build_tree_df(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))

  expect_s3_class(res, "tbl_tree")
  # offspring must work on the annotated table
  root <- ape::Ntip(tree) + 1L
  expect_equal(nrow(tidytree::offspring(res, root, tiponly = TRUE)),
               ape::Ntip(tree))
})
