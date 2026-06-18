test_that("node_dates_from_timetree returns observed tip dates unchanged", {
  tree      <- make_intro_tree()
  tip_dates <- make_intro_tip_dates(tree)

  res <- node_dates_from_timetree(tree, tip_dates)

  tip_res <- dplyr::filter(res, node <= ape::Ntip(tree))
  # join back by tip order
  expect_equal(
    tip_res$inferred_date[order(tip_res$node)],
    tip_dates$date
  )
})

test_that("node_dates_from_timetree dates internal nodes before their tips", {
  tree      <- make_intro_tree()
  tip_dates <- make_intro_tip_dates(tree)

  res  <- node_dates_from_timetree(tree, tip_dates)
  ntip <- ape::Ntip(tree)
  root_date <- res$inferred_date[res$node == ntip + 1L]

  # Root must be no later than the earliest tip (it is ancestral).
  expect_lte(root_date, min(tip_dates$date))
})

test_that("node_dates_from_timetree is unit-agnostic (scaling branch lengths is invariant)", {
  tree      <- make_intro_tree()
  tip_dates <- make_intro_tip_dates(tree)

  tree_scaled <- tree
  tree_scaled$edge.length <- tree$edge.length * 365  # e.g. years -> days

  d1 <- node_dates_from_timetree(tree, tip_dates)
  d2 <- node_dates_from_timetree(tree_scaled, tip_dates)

  # Predicted dates should match within rounding regardless of branch unit.
  expect_true(all(abs(as.numeric(d1$inferred_date - d2$inferred_date)) <= 1))
})

test_that("node_dates_from_timetree errors on a tree without branch lengths", {
  tree <- ape::read.tree(text = "((a,b),c);")
  td   <- tibble::tibble(tip_label = c("a", "b", "c"),
                         date = as.Date(c("2020-01-01", "2021-01-01", "2022-01-01")))
  expect_error(node_dates_from_timetree(tree, td), "branch length")
})

test_that("node_dates_from_timetree does not warn on a clocklike time tree", {
  tree      <- make_intro_tree()
  tip_dates <- make_intro_tip_dates(tree)
  expect_no_warning(node_dates_from_timetree(tree, tip_dates))
})

test_that("node_dates_from_timetree warns when temporal signal is weak", {
  tree <- make_intro_tree()
  # dates that zig-zag against depth order -> near-zero R2 (no clock signal)
  d <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
  scrambled <- tibble::tibble(
    tip_label = tree$tip.label,
    date = as.Date(c("2019-01-01", "2022-01-01",
                     "2019-01-01", "2022-01-01"))[rank(d, ties.method = "first")]
  )

  expect_warning(node_dates_from_timetree(tree, scrambled), "temporal signal")
})

test_that("node_dates_from_timetree errors when all tips are equidistant from root", {
  # ultrametric tree: every tip at depth 2 -> no spread in the predictor
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  td   <- tibble::tibble(
    tip_label = c("a", "b", "c", "d"),
    date = as.Date(c("2019-01-01", "2020-01-01", "2021-01-01", "2022-01-01"))
  )
  expect_error(node_dates_from_timetree(tree, td), "equidistant")
})
