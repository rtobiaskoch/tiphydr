make_intro_df <- function() {
  tree <- make_intro_tree()
  build_tree_df(tree, make_intro_node_probs(tree), make_intro_tip_dates(tree))
}

test_that("detect_introductions finds the established clade and the singleton", {
  res <- detect_introductions(make_intro_df())

  # regionB established clade (2 tips) + regionC singleton (1 tip) = 3 tip rows
  expect_equal(nrow(res), 3L)

  clades <- dplyr::distinct(res, intro_clade_id, deme, clade_size,
                            inferred_intro_source)

  # established clade: regionB, size 2, sourced from regionA
  b <- dplyr::filter(clades, deme == "regionB")
  expect_equal(b$clade_size, 2L)
  expect_equal(b$inferred_intro_source, "regionA")

  # singleton: regionC, size 1, sourced from regionA
  cc <- dplyr::filter(clades, deme == "regionC")
  expect_equal(cc$clade_size, 1L)
  expect_equal(cc$inferred_intro_source, "regionA")
})

test_that("detect_introductions reports the clade's most recent tip date", {
  res <- detect_introductions(make_intro_df())
  clades <- dplyr::distinct(res, deme, last_sample_date)

  # regionB clade holds tB1 (2020) and tB2 (2021) -> last = 2021
  expect_equal(clades$last_sample_date[clades$deme == "regionB"],
               as.Date("2021-01-01"))
  # regionC singleton: last sample date is its own date (2022)
  expect_equal(clades$last_sample_date[clades$deme == "regionC"],
               as.Date("2022-01-01"))
})

test_that("detect_introductions numbers clades chronologically", {
  res <- detect_introductions(make_intro_df())
  ord <- dplyr::distinct(res, intro_clade_id, inferred_intro_date) |>
    dplyr::arrange(intro_clade_id)
  # clade ids increase with introduction date
  expect_false(is.unsorted(ord$inferred_intro_date))
})

test_that("detect_introductions confidence filter drops low-probability clades", {
  df <- make_intro_df()
  # Halve the regionB node's confidence so it falls below 0.95.
  df$confidence_state[df$inferred_state == "regionB" & !df$is_tip] <- 0.5

  # all internal transitions now fall below threshold -> diagnostic warning
  expect_warning(
    res <- detect_introductions(df, confidence = 0.95),
    "none reach confidence"
  )
  # established regionB clade is gone; only the regionC singleton remains
  expect_false("regionB" %in% res$deme)
  expect_true("regionC" %in% res$deme)
})

test_that("detect_introductions returns expected columns", {
  res <- detect_introductions(make_intro_df())
  expect_named(res, c("tipname", "deme", "intro_clade_id", "inferred_intro_date",
                      "last_sample_date", "inferred_intro_source",
                      "inferred_intro_source_probability",
                      "intro_confidence_state",
                      "clade_size", "intro_node"))
})
