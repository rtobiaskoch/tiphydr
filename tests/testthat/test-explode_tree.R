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

test_that("explode_tree (simmap source) excludes clades below the posterior_support floor", {
  tree <- make_intro_tree()
  clade_dwell <- make_intro_clade_dwell(tree)
  tip_membership <- make_intro_tip_membership(tree)

  # Add a low-support phantom clade at a real tree node (node_AC) -- grounded
  # in the tree so tree-extraction wouldn't choke if the filter failed to
  # exclude it, but posterior_support = 0.2 should keep it out entirely. Give
  # it a REAL candidate tip (tA1, which has zero rows in the baseline
  # tip_membership -- it never transitions) with membership_prob = 1: if the
  # establishment filter were broken or absent, tA1 would wrongly appear as
  # a regionD introduction. Without this candidate row the test would pass
  # trivially even with a broken filter, since there'd be nothing to
  # wrongly include.
  node_AC <- ape::getMRCA(
    tree,
    grep("^tA1|^tC1", tree$tip.label, value = TRUE)
  )
  clade_dwell <- dplyr::bind_rows(
    clade_dwell,
    tibble::tibble(
      intro_clade_id = 3L,
      intro_node = node_AC,
      deme = "regionD",
      posterior_support = 0.2,
      inferred_intro_source = "regionA",
      inferred_intro_source_probability = 1,
      clade_size = 1L,
      inferred_intro_date = as.Date("2017-06-01"),
      last_sample_date = as.Date("2019-01-01")
    )
  )
  tip_membership <- dplyr::bind_rows(
    tip_membership,
    tibble::tibble(
      tipname = grep("^tA1", tree$tip.label, value = TRUE),
      intro_node = node_AC,
      deme = "regionD",
      membership_prob = 1,
      is_modal = TRUE
    )
  )

  res <- explode_tree(
    tree,
    clade_dwell = clade_dwell,
    tip_membership = tip_membership
  )

  expect_false("regionD" %in% res$introductions$deme)
  expect_false(grep("^tA1", tree$tip.label, value = TRUE) %in% res$introductions$tipname)
  # baseline clades unaffected
  expect_equal(nrow(res$introductions), 3L)
})

test_that("explode_tree (simmap source) warns and returns zero rows when nothing clears confidence", {
  tree <- make_intro_tree()
  clade_dwell <- make_intro_clade_dwell(tree) |>
    dplyr::mutate(posterior_support = 0.1)

  expect_warning(
    res <- explode_tree(
      tree,
      clade_dwell = clade_dwell,
      tip_membership = make_intro_tip_membership(tree)
    ),
    "none reach"
  )
  expect_equal(nrow(res$introductions), 0L)
  expect_length(res$trees, 0L)
})

test_that("explode_tree (simmap source) drops a tip whose membership_prob is below confidence", {
  tree <- make_intro_tree()
  tip_membership <- make_intro_tip_membership(tree) |>
    dplyr::mutate(
      membership_prob = dplyr::if_else(
        .data$tipname == grep("^tB2", tree$tip.label, value = TRUE),
        0.3, # below default confidence = 0.5
        .data$membership_prob
      )
    )

  res <- explode_tree(
    tree,
    clade_dwell = make_intro_clade_dwell(tree),
    tip_membership = tip_membership
  )

  # tB2 dropped, tB1 kept, tC1 (singleton) unaffected -- 2 rows, not 3.
  expect_equal(nrow(res$introductions), 2L)
  expect_false(grep("^tB2", tree$tip.label, value = TRUE) %in% res$introductions$tipname)
  expect_true(grep("^tB1", tree$tip.label, value = TRUE) %in% res$introductions$tipname)

  # clade_size still reflects the simmap's aggregate estimate (2), even
  # though only 1 tip actually survives the per-tip floor in this output --
  # intentional: clade_size is the model's median across replicates, not a
  # recount of this particular filtered view. See the spec's "documented
  # limitation" note.
  b_row <- dplyr::filter(res$introductions, .data$deme == "regionB")
  expect_equal(nrow(b_row), 1L)
  expect_equal(b_row$clade_size, 2L)
})

test_that("explode_tree (simmap source) recomputes modal pick among established candidates, not the raw global argmax", {
  tree <- make_intro_tree()
  clade_dwell <- make_intro_clade_dwell(tree)
  node_AC <- ape::getMRCA(
    tree,
    grep("^tA1|^tC1", tree$tip.label, value = TRUE)
  )
  # A low-support phantom clade (won't establish) that tB2 has HIGHER
  # membership_prob in than its true (established) regionB clade -- if the
  # implementation used the raw precomputed is_modal flag instead of
  # recomputing among survivors, tB2 would be wrongly dropped entirely.
  clade_dwell <- dplyr::bind_rows(
    clade_dwell,
    tibble::tibble(
      intro_clade_id = 3L,
      intro_node = node_AC,
      deme = "regionD",
      posterior_support = 0.2,
      inferred_intro_source = "regionA",
      inferred_intro_source_probability = 1,
      clade_size = 1L,
      inferred_intro_date = as.Date("2017-06-01"),
      last_sample_date = as.Date("2021-01-01")
    )
  )
  tip_membership <- make_intro_tip_membership(tree) |>
    dplyr::mutate(
      membership_prob = dplyr::if_else(
        .data$tipname == grep("^tB2", tree$tip.label, value = TRUE),
        0.6, # still clears confidence, but no longer 1.0
        .data$membership_prob
      ),
      # tB2's regionB row is no longer its raw global modal once the regionD
      # row below (0.7) is added -- is_modal reflects that honestly, even
      # though .introductions_from_simmap() never reads this column itself
      # (it recomputes modal among established survivors instead; see the
      # docstring on that helper).
      is_modal = dplyr::if_else(
        .data$tipname == grep("^tB2", tree$tip.label, value = TRUE),
        FALSE,
        .data$is_modal
      )
    ) |>
    dplyr::bind_rows(
      tibble::tibble(
        tipname = grep("^tB2", tree$tip.label, value = TRUE),
        intro_node = node_AC,
        deme = "regionD",
        membership_prob = 0.7, # higher than its regionB row -- global argmax
        is_modal = TRUE
      )
    )

  res <- explode_tree(
    tree,
    clade_dwell = clade_dwell,
    tip_membership = tip_membership
  )

  tb2 <- dplyr::filter(res$introductions, .data$tipname == grep("^tB2", tree$tip.label, value = TRUE))
  expect_equal(nrow(tb2), 1L)
  expect_equal(tb2$deme, "regionB")
})
