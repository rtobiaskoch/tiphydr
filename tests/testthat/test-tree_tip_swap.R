test_that("tree_tip_swap renames all tips when every id matches a tip", {
  tree <- make_test_tree()
  id   <- tree$tip.label
  new  <- c("NY10_001", "NS2A_002", "NS4B_003", "OTHER_004")

  result <- tree_tip_swap(tree, id = id, newnames = new)

  expect_setequal(result$tip.label, new)
})

test_that("tree_tip_swap preserves Ntip and Nnode", {
  tree <- make_test_tree()

  result <- tree_tip_swap(
    tree,
    id       = tree$tip.label,
    newnames = paste0("Sample_", seq_along(tree$tip.label))
  )

  expect_equal(ape::Ntip(result),  ape::Ntip(tree))
  expect_equal(ape::Nnode(result), ape::Nnode(tree))
})

test_that("tree_tip_swap errors on unmatched id when match = TRUE", {
  tree <- make_test_tree()
  id   <- c(tree$tip.label[1:3], "WNV|2099|XX|GHOST_999")  # GHOST not in tree
  new  <- c("A", "B", "C", "D")

  expect_error(
    tree_tip_swap(tree, id = id, newnames = new, match = TRUE),
    regexp = "GHOST_999"
  )
})

test_that("tree_tip_swap skips unmatched id silently when match = FALSE", {
  tree <- make_test_tree()
  id   <- c(tree$tip.label[1:2], "WNV|2099|XX|GHOST_999")
  new  <- c("Sample1", "Sample2", "SampleGhost")

  result <- tree_tip_swap(tree, id = id, newnames = new, match = FALSE)

  # renamed
  expect_true("Sample1" %in% result$tip.label)
  expect_true("Sample2" %in% result$tip.label)
  # unreferenced tips kept as-is
  expect_true("WNV|2022|WY|NS4B_003"  %in% result$tip.label)
  expect_true("WNV|2019|CO|OTHER_004" %in% result$tip.label)
  # ghost id was ignored — not injected into the tree
  expect_false("SampleGhost" %in% result$tip.label)
})

test_that("tree_tip_swap errors on duplicate id values", {
  tree <- make_test_tree()
  id   <- c(tree$tip.label[1], tree$tip.label[1], tree$tip.label[2])
  new  <- c("Sample1", "SampleDupe", "Sample2")

  expect_error(
    tree_tip_swap(tree, id = id, newnames = new),
    regexp = "duplicate"
  )
})

test_that("tree_tip_swap errors when id and newnames have different lengths", {
  tree <- make_test_tree()

  expect_error(
    tree_tip_swap(tree, id = tree$tip.label, newnames = c("A", "B")),
    regexp = "length"
  )
})
