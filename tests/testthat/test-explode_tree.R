test_that("explode_tree returns introductions and one subtree per multi-tip clade", {
  tree <- make_intro_tree()
  res  <- explode_tree(tree, make_intro_node_probs(tree))

  expect_named(res, c("introductions", "trees"))

  # 3 kept tips (regionB clade of 2 + regionC singleton)
  expect_equal(nrow(res$introductions), 3L)
  # internal intro_node column is not exposed
  expect_false("intro_node" %in% names(res$introductions))

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
  expect_error(explode_tree("not_a_tree", make_intro_node_probs(tree)), "phylo")
})

test_that("explode_tree errors when node_probs lacks a node column", {
  tree <- make_intro_tree()
  expect_error(
    explode_tree(tree, data.frame(regionA = 1, regionB = 0)),
    "node"
  )
})

test_that("explode_tree errors when the date field cannot be parsed", {
  tree <- ape::read.tree(
    text = "((tA|notadate|regionA:1,tB|2020-01-01|regionB:1):1,tC|2021-01-01|regionC:2);"
  )
  expect_error(
    explode_tree(tree, make_intro_node_probs(make_intro_tree())),
    "date"
  )
})
