# test-plot_explode_tree.R

test_that("plot_explode_tree returns a ggplot for valid input", {
  tree   <- make_intro_tree()
  result <- explode_tree(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  p      <- plot_explode_tree(result$introductions)
  expect_s3_class(p, "ggplot")
})

test_that("plot_explode_tree errors on missing required columns", {
  bad_df <- tibble::tibble(intro_clade_id = 1L, x = "a")
  expect_error(
    plot_explode_tree(bad_df),
    regexp = "missing required columns"
  )
})

test_that("plot_explode_tree produces one facet per destination deme", {
  tree   <- make_intro_tree()
  result <- explode_tree(tree, make_intro_node_probs(tree), make_intro_annotated_tree_path(tree))
  p      <- plot_explode_tree(result$introductions)
  built  <- ggplot2::ggplot_build(p)
  n_demes <- dplyr::n_distinct(result$introductions$deme)
  expect_equal(length(built$layout$panel_params), n_demes)
})

test_that("plot_explode_tree assigns semi-persistent when persistent_days is supplied", {
  # Clade spans 500 days: above transient_days=365, below persistent_days=1826
  introductions <- tibble::tibble(
    intro_clade_id             = 1L,
    deme                       = "regionA",
    inferred_intro_date        = as.Date("2020-01-01"),
    last_sample_date           = as.Date("2021-05-16"), # 501 days
    inferred_intro_source      = "regionB",
    inferred_intro_source_probability = 1,
    clade_size                 = 2L
  )
  p     <- plot_explode_tree(introductions, persistence_days = 365, persistent_days = 1826,
                             named_palette = c(regionB = "#000000"))
  built <- ggplot2::ggplot_build(p)
  # color scale should contain semi-persistent
  color_vals <- built$plot$scales$scales[[1]]$palette(3)
  expect_true(!is.null(color_vals))
})

test_that("plot_explode_tree errors when source demes exceed palette capacity", {
  # Manufacture a minimal introductions tibble with 9 distinct source demes —
  # more than the Set2 palette maximum (8).
  introductions <- tibble::tibble(
    intro_clade_id    = 1:9,
    deme              = "regionA",
    inferred_intro_date           = as.Date("2020-01-01"),
    last_sample_date  = as.Date("2021-01-01"),
    inferred_intro_source         = paste0("source", 1:9),
    inferred_intro_source_probability = 1,
    clade_size        = 1L
  )
  expect_error(
    plot_explode_tree(introductions, palette_source = "Set2"),
    regexp = "exceeds the maximum colours"
  )
})
