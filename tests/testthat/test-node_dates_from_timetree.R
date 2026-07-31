test_that("node_dates_from_timetree returns annotated tip dates unchanged", {
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(tree)
  tip_dates <- make_intro_tip_dates(tree)

  res <- node_dates_from_timetree(tree, path)

  tip_res <- dplyr::filter(res, node <= ape::Ntip(tree))
  expect_equal(
    tip_res$inferred_date[order(tip_res$node)],
    tip_dates$date
  )
})

test_that("node_dates_from_timetree dates internal nodes before their tips", {
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(tree)
  tip_dates <- make_intro_tip_dates(tree)

  res  <- node_dates_from_timetree(tree, path)
  ntip <- ape::Ntip(tree)
  root_date <- res$inferred_date[res$node == ntip + 1L]

  # Root must be no later than the earliest tip (it is ancestral).
  expect_lte(root_date, min(tip_dates$date))
})

test_that("node_dates_from_timetree parses a bare-year annotation as that year's midpoint", {
  tree <- make_intro_tree()
  # tB1's own annotation reduced to a bare year -- LSD2's shape when the
  # fitted date has no finer resolution.
  path <- make_intro_annotated_tree_path(
    tree,
    tip_dates_override = c("tB1|2020-01-01|regionB" = "2020")
  )

  res <- node_dates_from_timetree(tree, path)
  tB1_date <- res$inferred_date[res$node == which(tree$tip.label == "tB1|2020-01-01|regionB")]

  expect_equal(tB1_date, as.Date("2020-07-02"))
})

test_that("node_dates_from_timetree parses a TreeTime decimal-year annotation", {
  tree <- make_intro_tree()
  # tA1's own annotation as a TreeTime-style decimal year instead of ISO.
  path <- make_intro_annotated_tree_path(
    tree,
    tip_dates_override = c("tA1|2019-01-01|regionA" = "2019.42")
  )

  res <- node_dates_from_timetree(tree, path)
  tA1_date <- res$inferred_date[res$node == which(tree$tip.label == "tA1|2019-01-01|regionA")]

  # 2019 is not a leap year (365 days); day 0.42*365 = 153.3 -> round to 153.
  expect_equal(tA1_date, as.Date("2019-01-01") + round(0.42 * 365))
})

test_that("node_dates_from_timetree returns a non-clocklike annotated date exactly, not a model prediction", {
  # This is the actual regression test for the bug this design fixes: the old
  # median-anchor/regression approach MODELED every node's date from a single
  # uniform-clock assumption (depth * time_unit_days + one shared anchor),
  # which can (and on real data, did) place an ancestor's modeled date later
  # than one of its own descendant tips' true observed dates whenever a tip's
  # real date didn't fit that uniform model. The new design never models
  # anything -- it only reads whatever the dating tool itself already wrote.
  # Prove that by annotating a tip with a date that is NOT what the
  # depth-based clock model of make_intro_tree() would predict (tC1's "true"
  # clock-implied date is 2022-01-01; give it a wildly different annotated
  # date instead) and confirm the function returns exactly that value,
  # unmodified -- there is no regression/anchor step left that could
  # "correct" it.
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(
    tree,
    tip_dates_override = c("tC1|2022-01-01|regionC" = "1999-03-15")
  )

  res <- node_dates_from_timetree(tree, path)
  tC1_date <- res$inferred_date[res$node == which(tree$tip.label == "tC1|2022-01-01|regionC")]

  expect_equal(tC1_date, as.Date("1999-03-15"))
})

test_that("node_dates_from_timetree errors when tree and annotated_tree_path topology differ", {
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(tree)

  pruned <- ape::drop.tip(tree, "tC1|2022-01-01|regionC")
  expect_error(
    node_dates_from_timetree(pruned, path),
    "do not share the same topology"
  )
})

test_that("node_dates_from_timetree errors when a node has no annotated date", {
  tree <- make_intro_tree()
  path <- make_intro_annotated_tree_path(tree)
  # Corrupt the file so one node's [&date=...] comment is gone entirely.
  lines <- readLines(path)
  lines <- sub("\\[&date=2020-01-01\\]", "", lines, fixed = FALSE)
  writeLines(lines, path)

  expect_error(node_dates_from_timetree(tree, path), "no annotated date found")
})
