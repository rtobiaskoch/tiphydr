# test-plot_clade_persistence.R

test_that("plot_clade_persistence returns a ggplot for valid input", {
  tree <- make_intro_tree()
  p <- plot_clade_persistence(
    make_intro_raw_clade_dwell(tree),
    make_intro_clade_dwell(tree)
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_clade_persistence errors on missing raw_clade_dwell columns", {
  tree <- make_intro_tree()
  bad_df <- tibble::tibble(intro_node = 1L, x = "a")
  expect_error(
    plot_clade_persistence(bad_df, make_intro_clade_dwell(tree)),
    "raw_clade_dwell is missing required columns"
  )
})

test_that("plot_clade_persistence errors on missing clade_dwell columns", {
  tree <- make_intro_tree()
  bad_df <- tibble::tibble(intro_clade_id = 1L, x = "a")
  expect_error(
    plot_clade_persistence(make_intro_raw_clade_dwell(tree), bad_df),
    "clade_dwell is missing required columns"
  )
})

test_that("plot_clade_persistence produces one facet per destination deme", {
  tree <- make_intro_tree()
  p <- plot_clade_persistence(
    make_intro_raw_clade_dwell(tree),
    make_intro_clade_dwell(tree)
  )
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(built$layout$panel_params), 2L) # regionB, regionC
})

test_that("plot_clade_persistence renders a degenerate box for a singleton-support clade without error", {
  tree <- make_intro_tree()
  # tip_C has exactly 1 replicate draw in make_intro_raw_clade_dwell() --
  # confirms geom_boxplot() doesn't error on n = 1.
  expect_no_error(
    p <- plot_clade_persistence(
      make_intro_raw_clade_dwell(tree),
      make_intro_clade_dwell(tree)
    )
  )
  built <- ggplot2::ggplot_build(p)
  expect_true(nrow(built$data[[1]]) > 0L)
})

test_that("plot_clade_persistence's row ordering matches plot_explode_tree()'s for the same clade_dwell", {
  tree <- make_intro_tree()
  base_clade <- make_intro_clade_dwell(tree) |> dplyr::filter(deme == "regionB")

  # 3 clades, all in ONE deme (regionB), with distinct out-of-insertion-order
  # inferred_intro_date values, so .order_clades()'s arrange-by-date produces a
  # non-trivial ordering. Single-deme fixture avoids any facet/PANEL cross-deme
  # y-value collision when reading raw y positions off the built plot objects below.
  # intro_node values MUST be distinct across the three rows -- plot_clade_persistence()
  # inner-joins raw_clade_dwell to the ordered clade table by (intro_node, deme);
  # reusing one intro_node across multiple intro_clade_id rows would fan the join out.
  clade_dwell <- dplyr::bind_rows(
    base_clade,
    tibble::tibble(
      intro_clade_id = 3L, intro_node = 9001L, deme = "regionB",
      posterior_support = 0.9, inferred_intro_source = "regionA",
      inferred_intro_source_probability = 1, clade_size = 1L,
      inferred_intro_date = as.Date("2015-01-01"), last_sample_date = as.Date("2015-06-01")
    ),
    tibble::tibble(
      intro_clade_id = 4L, intro_node = 9002L, deme = "regionB",
      posterior_support = 0.9, inferred_intro_source = "regionA",
      inferred_intro_source_probability = 1, clade_size = 1L,
      inferred_intro_date = as.Date("2023-01-01"), last_sample_date = as.Date("2023-06-01")
    )
  )
  raw_clade_dwell <- dplyr::bind_rows(
    dplyr::filter(make_intro_raw_clade_dwell(tree), deme == "regionB"),
    tibble::tibble(sim_id = 1L, intro_node = 9001L, deme = "regionB", persistence_days = 150),
    tibble::tibble(sim_id = 1L, intro_node = 9002L, deme = "regionB", persistence_days = 150)
  )

  # Canonical, single source of truth -- called exactly once.
  expected_order <- .order_clades(clade_dwell)
  # Named vector: intro_clade_id (as name) -> expected local_clade_num, for direct
  # per-clade lookup rather than a set-membership check that can't detect misordering.
  expected_y_by_clade <- stats::setNames(
    expected_order$local_clade_num, expected_order$intro_clade_id
  )

  expl_built <- ggplot2::ggplot_build(plot_explode_tree(clade_dwell))
  pers_built <- ggplot2::ggplot_build(plot_clade_persistence(raw_clade_dwell, clade_dwell))

  # ggplot2 assigns each layer's `group` column as the rank of each row's discrete
  # grouping value among the sorted unique values present in that layer's data
  # (ggplot2:::add_group() / vctrs::vec_group_id()). Since group = intro_clade_id was
  # mapped explicitly in both functions, `group` here is rank(intro_clade_id), not
  # intro_clade_id itself -- recover the real intro_clade_id via the same sorted-unique
  # lookup, so the comparison below is keyed by clade identity, not position alone.
  sorted_ids <- sort(unique(clade_dwell$intro_clade_id))

  extract_y_by_clade <- function(built) {
    d <- built$data[[1]]
    stats::setNames(d$y, sorted_ids[d$group])
  }

  expl_y_by_clade <- extract_y_by_clade(expl_built)
  pers_y_by_clade <- extract_y_by_clade(pers_built)

  # Index both by the same name order so expect_equal compares position-for-position
  # by clade identity, not just as sets.
  ord_ids <- as.character(sorted_ids)
  expect_equal(expl_y_by_clade[ord_ids], expected_y_by_clade[ord_ids])
  expect_equal(pers_y_by_clade[ord_ids], expected_y_by_clade[ord_ids])
})
