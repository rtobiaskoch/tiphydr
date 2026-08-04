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
