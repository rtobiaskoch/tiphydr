# Synthetic mutation table exercising the base/partial/hybrid naming grammar
# independently of the real wnv_mut_orf data, to confirm the function is
# generic to any table following the same convention.
synthetic_muts <- function() {
  tibble::tribble(
    ~lineage,   ~pos, ~residue,
    "FOO",         1,      "A",
    "FOO",         2,      "A",
    "FOO",         3,      "A",
    "FOO_X",       1,      "A",
    "FOO_X",       2,      "A",
    "BAR",         4,      "A",
    "BAR",         5,      "A",
    "FOO-BAR",     1,      "A",
    "FOO-BAR",     2,      "A",
    "FOO-BAR",     4,      "A",
    "FOO-BAR",     5,      "A"
  )
}

synthetic_base_colors <- function() {
  tibble::tribble(
    ~lineage, ~color,
    "FOO",    "#FF0000",
    "BAR",    "#0000FF"
  )
}

test_that("a bare base family gets its exact base color", {
  result <- assign_lineage_colors(synthetic_muts(), synthetic_base_colors())
  expect_equal(result$color[result$lineage == "FOO"], "#FF0000")
  expect_equal(result$color[result$lineage == "BAR"], "#0000FF")
})

test_that("a single-family partial gets a lighter tint of its family color", {
  result <- assign_lineage_colors(synthetic_muts(), synthetic_base_colors())
  foo_x <- result$color[result$lineage == "FOO_X"]
  expect_false(foo_x == "#FF0000")
  # Lighter means higher RGB channel values on average than the base red.
  base_rgb <- as.vector(grDevices::col2rgb("#FF0000"))
  tint_rgb <- as.vector(grDevices::col2rgb(foo_x))
  expect_true(mean(tint_rgb) > mean(base_rgb))
})

test_that("a hybrid gets a blend weighted toward the higher-mutation-count side", {
  # FOO has 3 rows, BAR has 2 rows -> hybrid should be closer to FOO's red
  # than to BAR's blue.
  result <- assign_lineage_colors(synthetic_muts(), synthetic_base_colors())
  hybrid_rgb <- as.vector(grDevices::col2rgb(result$color[result$lineage == "FOO-BAR"]))
  foo_rgb <- as.vector(grDevices::col2rgb("#FF0000"))
  bar_rgb <- as.vector(grDevices::col2rgb("#0000FF"))
  dist_to_foo <- sqrt(sum((hybrid_rgb - foo_rgb)^2))
  dist_to_bar <- sqrt(sum((hybrid_rgb - bar_rgb)^2))
  expect_true(dist_to_foo < dist_to_bar)
})

test_that("unknown and ambiguous always appear with their given colors", {
  result <- assign_lineage_colors(
    synthetic_muts(), synthetic_base_colors(),
    unknown_color = "#111111", ambiguous_color = "#222222"
  )
  expect_equal(result$color[result$lineage == "unknown"], "#111111")
  expect_equal(result$color[result$lineage == "ambiguous"], "#222222")
})

test_that("muts missing required columns throws an error", {
  bad_muts <- tibble::tibble(lineage = "FOO", position = 1L, base = "A")
  expect_error(
    assign_lineage_colors(bad_muts, synthetic_base_colors()),
    "lineage.*pos.*residue"
  )
})

test_that("base_colors missing required columns throws an error", {
  bad_base <- tibble::tibble(lineage = "FOO", hex = "#FF0000")
  expect_error(
    assign_lineage_colors(synthetic_muts(), bad_base),
    "lineage.*color"
  )
})

# ── real WNV data regression checks ──────────────────────────────────────────

test_that("real wnv_mut_orf data: base families get their default colors", {
  result <- assign_lineage_colors(wnv_mut_orf)
  expect_equal(result$color[result$lineage == "NY10"], "#1B9E77")
  expect_equal(result$color[result$lineage == "WN02"], "#D95F02")
  expect_equal(result$color[result$lineage == "SW03"], "#7570B3")
  expect_equal(result$color[result$lineage == "NY99"], "#E7298A")
})

test_that("real wnv_mut_orf data: equal-fraction partials get identical colors", {
  # SW03_NS4A and SW03_NS5 both carry 2 of SW03's 3 mutations.
  result <- assign_lineage_colors(wnv_mut_orf)
  expect_equal(
    result$color[result$lineage == "SW03_NS4A"],
    result$color[result$lineage == "SW03_NS5"]
  )
})

test_that("real wnv_mut_orf data: output covers every muts lineage plus sentinels", {
  result <- assign_lineage_colors(wnv_mut_orf)
  expect_setequal(result$lineage, c(unique(wnv_mut_orf$lineage), "unknown", "ambiguous"))
})
