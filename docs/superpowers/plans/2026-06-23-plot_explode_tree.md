# plot_explode_tree Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `plot_explode_tree()` — a function that takes `explode_tree()$introductions` and produces a horizontal timeline plot: one row per introduction clade, with the introduction date (coloured by source deme) and last sample date (coloured by persistence), faceted by destination deme.

**Architecture:** Pure function over a tibble; aggregates tip-level rows to clade-level with `dplyr::distinct`, constructs two separated `geom_point` layers (one per date type) so `color` (persistence) and `fill` (source deme) map to independent scales without any extra packages. No modification to upstream functions.

**Tech Stack:** `ggplot2`, `dplyr`, `RColorBrewer` — all already imported. No new dependencies.

## Global Constraints

- All functions: `verb_noun()` snake_case names, explicit namespacing (`dplyr::select`, `ggplot2::aes`, etc.)
- `ggplot2::theme_classic()` on all plots
- Use `.data$` pronoun inside `aes()` for tidy eval
- `@export` every public function; run `devtools::document()` after
- Tests use fixtures from `tests/testthat/helper-fixtures.R` (`make_intro_tree()`, `make_intro_node_probs()`)
- Run tests with: `devtools::test(filter = "plot_explode_tree")`
- No new package dependencies — use `dplyr::bind_rows` instead of `tidyr::pivot_longer`

---

### Data contract (read before touching any code)

`explode_tree()$introductions` — one row **per tip** (not per clade). Required columns:

| Column | Type | Notes |
|--------|------|-------|
| `intro_clade_id` | integer | unique per introduction event |
| `deme` | character | destination deme |
| `inferred_intro_date` | Date | date of the ancestral introduction node |
| `last_sample_date` | Date | latest tip date in the clade |
| `inferred_intro_source` | character | parent (source) deme |
| `clade_size` | integer | number of tips in clade |

`plot_explode_tree()` collapses these to one row per clade with `dplyr::distinct`, then computes `persistence`.

---

### Task 1: Write `R/plot_explode_tree.R`

**Files:**
- Create: `R/plot_explode_tree.R`

**Interfaces:**
- Consumes: `explode_tree()$introductions` tibble (see data contract above)
- Produces: `ggplot` object, exportable as `plot_explode_tree(introductions, persistence_days, palette_persistence, palette_source)`

- [ ] **Step 1: Create the file with roxygen header and function signature**

```r
#' Plot introduction clade timelines from explode_tree output
#'
#' Visualises each introduction clade as a horizontal segment: from
#' \code{inferred_intro_date} to \code{last_sample_date}, with the
#' introduction point (filled diamond) coloured by \code{inferred_intro_source}
#' and the last-sample point (open circle) and segment coloured by persistence
#' (\code{last_sample_date - inferred_intro_date > persistence_days}).
#' One facet per destination deme.
#'
#' @param introductions A tibble from \code{\link{explode_tree}()$introductions}.
#'   Required columns: \code{intro_clade_id}, \code{deme},
#'   \code{inferred_intro_date}, \code{last_sample_date},
#'   \code{inferred_intro_source}, \code{clade_size}.
#' @param persistence_days Numeric scalar. Days a clade must span (intro date to
#'   last sample) to be labelled \code{"persistent"}. Default 365 (one year).
#' @param palette_persistence RColorBrewer palette for the persistence colour
#'   scale (line + last-sample point). Default \code{"Dark2"}.
#' @param palette_source RColorBrewer palette for the source-deme fill scale
#'   (introduction-date point). Default \code{"Set2"}.
#'
#' @return A \code{ggplot} object. Additional \pkg{ggplot2} layers can be
#'   added with \code{+}.
#' @export
plot_explode_tree <- function(introductions,
                              persistence_days    = 365,
                              palette_persistence = "Dark2",
                              palette_source      = "Set2") {
```

- [ ] **Step 2: Add input validation**

```r
  # --- Validate required columns -------------------------------------------
  required_cols <- c(
    "intro_clade_id", "deme", "inferred_intro_date",
    "last_sample_date", "inferred_intro_source", "clade_size"
  )
  missing_cols <- setdiff(required_cols, names(introductions))
  if (length(missing_cols) > 0L) {
    stop(
      "introductions is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ".\nPass the $introductions element from explode_tree()."
    )
  }
```

- [ ] **Step 3: Aggregate to clade level and compute persistence**

```r
  # Introductions tibble has one row per tip; collapse to one row per clade.
  # persistence: did the clade persist for longer than persistence_days?
  clade_tbl <- introductions |>
    dplyr::distinct(
      .data$intro_clade_id,
      .data$deme,
      .data$inferred_intro_date,
      .data$last_sample_date,
      .data$inferred_intro_source,
      .data$clade_size
    ) |>
    dplyr::mutate(
      persistence = dplyr::if_else(
        as.numeric(.data$last_sample_date - .data$inferred_intro_date) > persistence_days,
        "persistent",
        "transient"
      )
    )
```

- [ ] **Step 4: Create long format and split into per-date-type data frames**

```r
  # Long format: two rows per clade (one per date type).
  # Using bind_rows avoids a tidyr dependency.
  clade_long <- dplyr::bind_rows(
    dplyr::mutate(clade_tbl,
                  date      = .data$inferred_intro_date,
                  date_type = "inferred_intro_date"),
    dplyr::mutate(clade_tbl,
                  date      = .data$last_sample_date,
                  date_type = "last_sample_date")
  )

  # Separate layers: intro date (fill = source) vs last sample (color = persistence).
  # color and fill are independent scales in ggplot2 — no extra package needed.
  intro_pts  <- dplyr::filter(clade_long, .data$date_type == "inferred_intro_date")
  sample_pts <- dplyr::filter(clade_long, .data$date_type == "last_sample_date")
```

- [ ] **Step 5: Build and return the ggplot**

```r
  ggplot2::ggplot(clade_long,
                  ggplot2::aes(y = .data$intro_clade_id)) +
    # Segment colored by persistence (was the clade transient or long-lived?)
    ggplot2::geom_line(
      ggplot2::aes(x     = .data$date,
                   group = .data$intro_clade_id,
                   color = .data$persistence),
      linewidth = 1,
      alpha     = 0.8
    ) +
    # Last-sample point: open circle, outline colored by persistence
    ggplot2::geom_point(
      data = sample_pts,
      ggplot2::aes(x     = .data$date,
                   color = .data$persistence,
                   size  = .data$clade_size),
      shape = 21,
      fill  = "white"
    ) +
    # Introduction-date point: filled diamond, fill colored by source deme
    ggplot2::geom_point(
      data = intro_pts,
      ggplot2::aes(x    = .data$date,
                   fill = .data$inferred_intro_source,
                   size = .data$clade_size),
      shape = 22,
      color = "white"
    ) +
    ggplot2::scale_color_brewer(palette = palette_persistence, name = "Persistence") +
    ggplot2::scale_fill_brewer(palette = palette_source,       name = "Source deme") +
    ggplot2::scale_size_continuous(range = c(2.5, 8),          guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(.data$deme), scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Introduction clade") +
    ggplot2::theme_classic()
}
```

- [ ] **Step 6: Run `devtools::document()` and verify `NAMESPACE` contains `export(plot_explode_tree)`**

```bash
Rscript -e "devtools::document()"
grep "plot_explode_tree" NAMESPACE
```

Expected: `export(plot_explode_tree)`

- [ ] **Step 7: Smoke-test interactively (in R)**

```r
devtools::load_all(".")
tree   <- make_intro_tree()
result <- explode_tree(tree, make_intro_node_probs(tree))
plot_explode_tree(result$introductions)
```

Expected: a ggplot with two facets (regionB, regionC), a line per clade, diamond + circle points.

---

### Task 2: Write `tests/testthat/test-plot_explode_tree.R`

**Files:**
- Create: `tests/testthat/test-plot_explode_tree.R`

**Interfaces:**
- Consumes: `plot_explode_tree()` from Task 1; fixtures `make_intro_tree()`, `make_intro_node_probs()` from `helper-fixtures.R`
- Produces: passing testthat suite

- [ ] **Step 1: Write the three failing tests**

```r
# test-plot_explode_tree.R

test_that("plot_explode_tree returns a ggplot for valid input", {
  tree   <- make_intro_tree()
  result <- explode_tree(tree, make_intro_node_probs(tree))
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
  result <- explode_tree(tree, make_intro_node_probs(tree))
  p      <- plot_explode_tree(result$introductions)
  built  <- ggplot2::ggplot_build(p)
  n_demes <- dplyr::n_distinct(result$introductions$deme)
  expect_equal(length(built$layout$panel_params), n_demes)
})
```

- [ ] **Step 2: Run tests to verify they fail (function does not exist yet)**

```bash
Rscript -e "devtools::test(filter = 'plot_explode_tree')"
```

Expected: 3 failures — `could not find function "plot_explode_tree"`.

- [ ] **Step 3: Go implement Task 1, then return here**

- [ ] **Step 4: Run tests to verify they pass**

```bash
Rscript -e "devtools::test(filter = 'plot_explode_tree')"
```

Expected: `[ PASS  3 ]`

- [ ] **Step 5: Run full test suite — no regressions**

```bash
Rscript -e "devtools::test()"
```

Expected: all prior tests still pass.

- [ ] **Step 6: Commit**

```bash
git add R/plot_explode_tree.R tests/testthat/test-plot_explode_tree.R NAMESPACE
git commit -m "feat: add plot_explode_tree — introduction timeline plot"
```

---

## Design decisions logged here for driver review

| Decision | Choice made | Alternative |
|----------|-------------|-------------|
| Date range | `inferred_intro_date` → `last_sample_date` | Add `fst_sample_date` to `detect_introductions` output (requires modifying that function + its tests) |
| Dual colour scales | `color` = persistence (line + sample point outline), `fill` = source (intro point interior) | `ggnewscale::new_scale_color()` — avoided new dependency |
| Persistence | Derived inside function: `last_sample - intro_date > persistence_days` | Pass as column in `introductions` |
| Facet scales | `scales = "free_y"` (shared x-axis for date comparison across demes) | `scales = "fixed"` (y-axis gaps between demes) |
| `fst_sample_date` | Not shown — not in `detect_introductions` output | Reference Rmd uses it; would require adding per-tip `sample_date` back to output |
