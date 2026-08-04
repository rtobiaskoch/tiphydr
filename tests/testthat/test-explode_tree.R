test_that("explode_tree returns introductions and one subtree per multi-tip clade", {
  tree <- make_intro_tree()
  res  <- explode_tree(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))

  expect_named(res, c("introductions", "trees"))

  # 3 kept tips (regionB clade of 2 + regionC singleton)
  expect_equal(nrow(res$introductions), 3L)
  # intro_node is exposed (scripts/dta_explode.R in the calling pipeline
  # depends on this column being present)
  expect_true("intro_node" %in% names(res$introductions))

  # exactly one subtree (the established regionB clade), with 2 tips
  expect_length(res$trees, 1L)
  expect_equal(ape::Ntip(res$trees[[1]]), 2L)
  expect_setequal(
    res$trees[[1]]$tip.label,
    grep("^tB", tree$tip.label, value = TRUE)
  )
})

test_that("explode_tree errors when tree is not a phylo object", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree("not_a_tree", make_intro_node_probs(tree), make_intro_annotated_tree_path(tree)),
    "phylo"
  )
})

test_that("explode_tree errors when node_probs lacks a node column", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree, data.frame(regionA = 1, regionB = 0), make_intro_annotated_tree_path(tree)),
    "node"
  )
})

test_that("explode_tree errors when annotated_tree_path topology does not match tree", {
  # Regression coverage for the same topology-mismatch guard tested directly
  # in test-node_dates_from_timetree.R -- confirms explode_tree() propagates
  # it rather than swallowing/masking the error.
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(tree)
  pruned <- ape::drop.tip(tree, "tC1|2022-01-01|regionC")

  expect_error(
    explode_tree(pruned, make_intro_node_probs(pruned), path),
    "do not share the same topology"
  )
})

test_that("explode_tree errors when neither source is supplied", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree),
    "Supply exactly one source"
  )
})

test_that("explode_tree errors when both sources are supplied", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(
      tree,
      node_probs = make_intro_node_probs(tree),
      annotated_tree_path = make_intro_annotated_tree_path(tree),
      clade_dwell = make_intro_clade_dwell(tree),
      tip_membership = make_intro_tip_membership(tree)
    ),
    "Supply exactly one source"
  )
})

test_that("explode_tree errors when only clade_dwell is supplied without tip_membership", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree, clade_dwell = make_intro_clade_dwell(tree)),
    "clade_dwell and tip_membership must both be supplied together"
  )
})

test_that("explode_tree errors when only tip_membership is supplied without clade_dwell", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree, tip_membership = make_intro_tip_membership(tree)),
    "clade_dwell and tip_membership must both be supplied together"
  )
})

test_that("explode_tree errors when node_probs is supplied without annotated_tree_path", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree, node_probs = make_intro_node_probs(tree)),
    "annotated_tree_path is required"
  )
})

test_that("explode_tree (simmap source) matches the node-marginal path's clade structure for the baseline fixture", {
  tree <- make_intro_tree()
  res <- explode_tree(
    tree,
    clade_dwell = make_intro_clade_dwell(tree),
    tip_membership = make_intro_tip_membership(tree)
  )

  expect_named(res, c("introductions", "trees"))
  expect_true("intro_node" %in% names(res$introductions))

  # Same structure as the node-marginal test: 3 kept tips (regionB clade of
  # 2 + regionC singleton), one multi-tip subtree.
  expect_equal(nrow(res$introductions), 3L)
  expect_length(res$trees, 1L)
  expect_equal(ape::Ntip(res$trees[[1]]), 2L)
  expect_setequal(
    res$trees[[1]]$tip.label,
    grep("^tB", tree$tip.label, value = TRUE)
  )

  # intro_confidence_state carries posterior_support for the simmap path.
  b_row <- dplyr::filter(res$introductions, .data$deme == "regionB")
  expect_equal(unique(b_row$intro_confidence_state), 0.9)
})
