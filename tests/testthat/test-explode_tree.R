# Write the in-memory fixtures to disk so explode_tree can read them as it would
# real TreeTime output.
write_intro_fixture <- function() {
  tree <- make_intro_tree()
  dir  <- withr::local_tempdir(.local_envir = parent.frame())

  tree_path  <- file.path(dir, "tree.nwk")
  probs_path <- file.path(dir, "node_probs.tsv")

  ape::write.tree(tree, tree_path)
  utils::write.table(make_intro_node_probs(tree), probs_path,
                     sep = "\t", row.names = FALSE, quote = FALSE)

  list(tree_path = tree_path, probs_path = probs_path)
}

test_that("explode_tree returns introductions and one subtree per multi-tip clade", {
  f   <- write_intro_fixture()
  res <- explode_tree(f$tree_path, f$probs_path)

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
    grep("^tB", make_intro_tree()$tip.label, value = TRUE)
  )
})

test_that("explode_tree errors on a missing input file", {
  f <- write_intro_fixture()
  expect_error(explode_tree("nope.nwk", f$probs_path), "does not exist")
  expect_error(explode_tree(f$tree_path, "nope.tsv"), "does not exist")
})

test_that("explode_tree errors when the date field cannot be parsed", {
  dir  <- withr::local_tempdir()
  tree <- ape::read.tree(text = "((tA|notadate|regionA:1,tB|2020-01-01|regionB:1):1,tC|2021-01-01|regionC:2);")
  tp   <- file.path(dir, "t.nwk"); ape::write.tree(tree, tp)
  pp   <- file.path(dir, "p.tsv")
  utils::write.table(make_intro_node_probs(make_intro_tree()), pp,
                     sep = "\t", row.names = FALSE, quote = FALSE)

  expect_error(explode_tree(tp, pp), "date")
})
