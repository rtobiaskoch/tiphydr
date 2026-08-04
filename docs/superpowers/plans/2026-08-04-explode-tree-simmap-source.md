# explode_tree() Dual-Source + plot_clade_persistence() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `tiphydr::explode_tree()` accept either node-marginal (`node_probs`) or simmap-branch (`clade_dwell` + `tip_membership`) input and return the identical `list(introductions, trees)` shape either way; add a new `plot_clade_persistence()` boxplot function alongside it.

**Architecture:** `explode_tree()` gains two new parameters and dispatches by which are supplied. The simmap branch is a new internal helper, `.introductions_from_simmap()`, that reduces `clade_dwell`/`tip_membership` down to a tibble with the exact same columns `detect_introductions()` produces, via a two-step confidence filter (clade-level `posterior_support`, then tip-level `membership_prob`, recomputed after the clade filter — not the precomputed `is_modal` flag). The existing tree-extraction code at the bottom of `explode_tree()` is untouched and consumes either branch's output identically. Row-ordering logic shared between `plot_explode_tree()` and the new `plot_clade_persistence()` is extracted into `.order_clades()`.

**Tech Stack:** R package (`tiphydr`), `dplyr`/`tibble` (native pipe `|>`, `.data$col` NSE), `ggplot2` (`ggplot2::geom_boxplot(orientation = "y")`), `testthat` edition 3, `ape`.

## Global Constraints

- Package: `/Users/user/Programming_Directory/tiphydr` (a separate git repo from the calling pipeline). Every commit in this plan happens there.
- Existing positional call sites (`explode_tree(tree, node_probs, annotated_tree_path)`, 3 positional args, used by `scripts/dta_explode.R` in the `wnv-ge_dta_stan` pipeline and by this package's own tests) must keep working unchanged — `node_probs`/`annotated_tree_path` stay parameters 2 and 3 in that order.
- `intro_node` **stays** in `explode_tree()`'s `$introductions` output for both source paths — `scripts/dta_explode.R` (a different repo) reads it directly and it's already a real column in production output (`runs/R028/R028_flyway_1_dta_explode.tsv`). Do not drop it. The existing docstring/test claiming it's absent is the actual bug being fixed (Task 1), not the code.
- No new package dependencies. `ggplot2` 4.0.2 is already a dependency and supports `orientation = "y"` natively (available since ggplot2 3.3).
- Style: native pipe `|>`, `.data$col` NSE per existing files (`R/explode_tree.R`, `R/detect_introductions.R`, `R/plot_explode_tree.R`) — follow, don't reintroduce `%>%`.
- `confidence` (default `0.5`) is the single shared threshold for both the clade-establishment filter (`posterior_support >= confidence`) and the tip-level filter (`membership_prob >= confidence`) — not two separate parameters.

---

### Task 1: Fix the `intro_node` docstring/test bug (baseline cleanup)

**Files:**
- Modify: `R/explode_tree.R:31` (Roxygen `@return` doc)
- Modify: `tests/testthat/test-explode_tree.R:10`

**Interfaces:**
- Consumes: nothing new.
- Produces: a green baseline (`devtools::test()` all-passing) before any new code is added, and confirms `intro_node`'s presence is the *intended* contract for every later task in this plan to build against.

- [ ] **Step 1: Confirm the current (buggy) state**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: 1 failure — `test-explode_tree.R:10` — `Expected "intro_node" %in% names(res$introductions) to be FALSE. actual: TRUE, expected: FALSE`.

- [ ] **Step 2: Fix the docstring**

In `R/explode_tree.R`, find this line inside the `@return` block:
```r
#'     \item{\code{introductions}}{A tibble, one row per kept tip (see
#'       \code{\link{detect_introductions}}), without the internal \code{intro_node}.}
```
Replace with:
```r
#'     \item{\code{introductions}}{A tibble, one row per kept tip (see
#'       \code{\link{detect_introductions}}), including \code{intro_node} (the
#'       ape node id of the introduction) so callers can join back to the
#'       source clade table.}
```

- [ ] **Step 3: Fix the test assertion**

In `tests/testthat/test-explode_tree.R`, change:
```r
  # internal intro_node column is not exposed
  expect_false("intro_node" %in% names(res$introductions))
```
to:
```r
  # intro_node is exposed (scripts/dta_explode.R in the calling pipeline
  # depends on this column being present)
  expect_true("intro_node" %in% names(res$introductions))
```

- [ ] **Step 4: Run tests to verify the whole suite is green**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(reporter = "summary")'`
Expected: 0 failures.

- [ ] **Step 5: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add R/explode_tree.R tests/testthat/test-explode_tree.R
git commit -m "Fix explode_tree() docstring/test: intro_node is intentionally exposed

scripts/dta_explode.R (wnv-ge_dta_stan) reads intro_node directly from
explode_tree()'s output and it's already a real column in production
output files — the code was correct, the docstring and test asserting
its absence were stale."
```

---

### Task 2: Add simmap-source fixtures to `helper-fixtures.R`

**Files:**
- Modify: `tests/testthat/helper-fixtures.R` (append new section)

**Interfaces:**
- Consumes: `make_intro_tree()` (existing, `tests/testthat/helper-fixtures.R:227`).
- Produces: `make_intro_clade_dwell(tree = make_intro_tree())` → tibble with columns `intro_clade_id, intro_node, deme, posterior_support, inferred_intro_source, inferred_intro_source_probability, clade_size, inferred_intro_date, last_sample_date`. `make_intro_tip_membership(tree = make_intro_tree())` → tibble with columns `tipname, intro_node, deme, membership_prob, is_modal`. Both encode the exact same clade structure as the existing node-marginal fixtures (`make_intro_node_probs()`): an established `regionB` clade of 2 (`tB1`, `tB2`) and a singleton `regionC` clade (`tC1`) — `tA1` never introduces, same as the node path.

- [ ] **Step 1: Write the fixtures**

Append to `tests/testthat/helper-fixtures.R`, after the existing `make_intro_annotated_tree_path()` function (before the "File writer" section):

```r
#' Simmap-shaped clade_dwell tibble matching make_intro_tree()'s structure.
#'
#' Encodes the same clade structure as make_intro_node_probs(): an established
#' regionB clade of size 2 (tB1, tB2) and a singleton regionC clade (tC1).
#' Dates match make_intro_annotated_tree_path()'s node-date scheme exactly
#' (node_B = 2018-01-01, tC1's own date = 2022-01-01) so cross-fixture
#' comparisons (e.g. against the node-marginal path's output) are meaningful.
#'
#' @param tree The phylo from make_intro_tree().
#' @return tibble, one row per clade.
make_intro_clade_dwell <- function(tree = make_intro_tree()) {
  lbl <- function(p) grep(p, tree$tip.label, value = TRUE)
  node_B <- ape::getMRCA(tree, c(lbl("^tB1"), lbl("^tB2")))
  tip_C <- match(lbl("^tC1"), tree$tip.label)

  tibble::tibble(
    intro_clade_id = 1:2,
    intro_node = c(node_B, tip_C),
    deme = c("regionB", "regionC"),
    posterior_support = c(0.9, 0.9),
    inferred_intro_source = c("regionA", "regionA"),
    inferred_intro_source_probability = c(1, 1),
    clade_size = c(2L, 1L),
    inferred_intro_date = as.Date(c("2018-01-01", "2022-01-01")),
    last_sample_date = as.Date(c("2021-01-01", "2022-01-01"))
  )
}

#' Simmap-shaped tip_membership tibble matching make_intro_clade_dwell().
#'
#' Every tip has membership_prob = 1 / is_modal = TRUE in its one true clade
#' (the "clean" baseline scenario) -- tA1 never appears, matching the node
#' path's fixture where tA1 never transitions.
#'
#' @param tree The phylo from make_intro_tree().
#' @return tibble, one row per (tipname, intro_node).
make_intro_tip_membership <- function(tree = make_intro_tree()) {
  lbl <- function(p) grep(p, tree$tip.label, value = TRUE)
  node_B <- ape::getMRCA(tree, c(lbl("^tB1"), lbl("^tB2")))
  tip_C <- match(lbl("^tC1"), tree$tip.label)

  tibble::tibble(
    tipname = c(lbl("^tB1"), lbl("^tB2"), lbl("^tC1")),
    intro_node = c(node_B, node_B, tip_C),
    deme = c("regionB", "regionB", "regionC"),
    membership_prob = c(1, 1, 1),
    is_modal = TRUE
  )
}

#' Per-replicate persistence_days draws matching make_intro_clade_dwell().
#'
#' node_B (the established clade) gets 5 replicate draws -- enough to form a
#' real distribution for a boxplot. tip_C (the singleton) gets exactly 1 draw,
#' deliberately, to exercise geom_boxplot()'s degenerate n=1 case.
#'
#' @param tree The phylo from make_intro_tree().
#' @return tibble, one row per (sim_id, intro_node).
make_intro_raw_clade_dwell <- function(tree = make_intro_tree()) {
  lbl <- function(p) grep(p, tree$tip.label, value = TRUE)
  node_B <- ape::getMRCA(tree, c(lbl("^tB1"), lbl("^tB2")))
  tip_C <- match(lbl("^tC1"), tree$tip.label)

  tibble::tibble(
    sim_id = c(1:5, 1L),
    intro_node = c(rep(node_B, 5), tip_C),
    deme = c(rep("regionB", 5), "regionC"),
    persistence_days = c(1080, 1050, 1100, 1090, 1070, 240)
  )
}
```

- [ ] **Step 2: Verify the fixtures load and produce the expected shape**

Run:
```bash
cd /Users/user/Programming_Directory/tiphydr && Rscript -e '
devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-fixtures.R")
tree <- make_intro_tree()
print(make_intro_clade_dwell(tree))
print(make_intro_tip_membership(tree))
print(make_intro_raw_clade_dwell(tree))
'
```
Expected: three tibbles print with no error; `make_intro_clade_dwell()` has 2 rows (`regionB`, `regionC`); `make_intro_tip_membership()` has 3 rows (`tB1`, `tB2`, `tC1` — no `tA1`); `make_intro_raw_clade_dwell()` has 6 rows (5 for the `regionB` node, 1 for the `regionC` tip).

- [ ] **Step 3: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add tests/testthat/helper-fixtures.R
git commit -m "Add simmap-source test fixtures matching make_intro_tree()'s structure"
```

---

### Task 3: `explode_tree()` — dual-source input validation

**Files:**
- Modify: `R/explode_tree.R`
- Modify: `tests/testthat/test-explode_tree.R`

**Interfaces:**
- Consumes: `make_intro_tree()`, `make_intro_node_probs()`, `make_intro_annotated_tree_path()`, `make_intro_clade_dwell()`, `make_intro_tip_membership()` (all from `helper-fixtures.R`).
- Produces: `explode_tree()`'s new signature `function(tree, node_probs = NULL, annotated_tree_path = NULL, clade_dwell = NULL, tip_membership = NULL, confidence = 0.5)` with validation only — the simmap branch itself is stubbed to `stop("not yet implemented")` in this task, filled in by Task 4-6.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-explode_tree.R`:

```r
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
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: 5 new failures (the validation doesn't exist yet — calls either error with base R's "argument is missing" or succeed/error for the wrong reason).

- [ ] **Step 3: Implement the validation**

In `R/explode_tree.R`, replace the function signature and the top of the body:

```r
explode_tree <- function(
  tree,
  node_probs = NULL,
  annotated_tree_path = NULL,
  clade_dwell = NULL,
  tip_membership = NULL,
  confidence = 0.5
) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be an ape::phylo object")
  }

  has_node <- !is.null(node_probs)
  has_simmap <- !is.null(clade_dwell) || !is.null(tip_membership)

  if (has_node && has_simmap) {
    stop(
      "Supply exactly one source: node_probs (+ annotated_tree_path), OR ",
      "both clade_dwell and tip_membership -- not both at once."
    )
  }
  if (!has_node && !has_simmap) {
    stop(
      "Supply exactly one source: node_probs (+ annotated_tree_path), OR ",
      "both clade_dwell and tip_membership."
    )
  }
  if (has_simmap && (is.null(clade_dwell) || is.null(tip_membership))) {
    stop("clade_dwell and tip_membership must both be supplied together.")
  }

  if (has_node) {
    if (is.null(annotated_tree_path)) {
      stop("annotated_tree_path is required when using node_probs.")
    }
    if (!is.data.frame(node_probs) || !"node" %in% names(node_probs)) {
      stop("node_probs must be a data frame with a 'node' column")
    }
    tree_df <- build_tree_df(tree, node_probs, annotated_tree_path)
    intros <- detect_introductions(tree_df, confidence = confidence)
  } else {
    stop("simmap source not yet implemented") # replaced in Task 4-6
  }

  # Build one subtree per multi-tip clade: descend to the introduction node, then
  # keep only the tips that passed the continued-transmission filter.
  multi <- intros |>
    dplyr::filter(.data$clade_size > 1) |>
    dplyr::arrange(.data$intro_clade_id)

  if (nrow(multi) == 0L) {
    trees <- list()
  } else {
    trees <- multi |>
      dplyr::group_by(.data$intro_clade_id) |>
      dplyr::group_map(function(rows, key) {
        clade <- ape::extract.clade(tree, rows$intro_node[[1]])
        ape::keep.tip(clade, rows$tipname)
      }) |>
      stats::setNames(unique(multi$intro_clade_id))
  }

  list(
    introductions = intros,
    trees = trees
  )
}
```

Also update the `@param` Roxygen block above the function to document the four new parameters (`clade_dwell`, `tip_membership`) and that `node_probs`/`annotated_tree_path` are now optional-but-required-together — match the existing doc style in the file.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: all pass, including the 5 new validation tests and the pre-existing node-path tests (still using the `has_node` branch, unchanged behavior).

- [ ] **Step 5: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add R/explode_tree.R tests/testthat/test-explode_tree.R
git commit -m "Add dual-source input validation to explode_tree()

Accepts either node_probs (+ annotated_tree_path) or clade_dwell +
tip_membership, exactly one. Simmap branch itself is stubbed pending
the next tasks."
```

---

### Task 4: `.introductions_from_simmap()` — happy-path equivalence to the node path

**Files:**
- Modify: `R/explode_tree.R`
- Modify: `tests/testthat/test-explode_tree.R`

**Interfaces:**
- Consumes: `make_intro_clade_dwell()`, `make_intro_tip_membership()` (Task 2).
- Produces: internal (unexported) `.introductions_from_simmap(clade_dwell, tip_membership, confidence)` → tibble with columns `tipname, deme, intro_clade_id, inferred_intro_date, last_sample_date, inferred_intro_source, inferred_intro_source_probability, intro_confidence_state, clade_size, intro_node` — identical column set to `detect_introductions()`'s output. Wired into `explode_tree()`'s `has_simmap` branch, replacing the Task 3 stub.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-explode_tree.R`:

```r
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: FAIL with "simmap source not yet implemented".

- [ ] **Step 3: Implement `.introductions_from_simmap()`**

In `R/explode_tree.R`, add this function above `explode_tree()` (it's a package-internal helper, no `@export`):

```r
#' Build a per-tip introductions tibble from simmap-aggregated clade tables
#'
#' Two-step confidence filter, mirroring detect_introductions()'s
#' `confidence_state >= confidence` gate but split across the two axes the
#' simmap data actually has: clade-level (posterior_support -- is this
#' introduction event real at all) and tip-level (membership_prob -- is this
#' specific tip's assignment to it reliable). See
#' docs/superpowers/specs/2026-08-04-explode-tree-simmap-source-design.md for
#' the full rationale.
#'
#' @param clade_dwell Tibble, one row per clade: intro_node, deme,
#'   posterior_support, inferred_intro_source,
#'   inferred_intro_source_probability, clade_size, inferred_intro_date,
#'   last_sample_date (the shape of dta_compute.R's _dta_clade_dwell.tsv).
#' @param tip_membership Tibble, one row per (tipname, intro_node): tipname,
#'   intro_node, deme, membership_prob (the shape of
#'   _dta_tip_membership.tsv) -- the FULL table, not pre-filtered to
#'   is_modal == TRUE rows, since the modal pick here is recomputed after
#'   establishment filtering, not read off the file's own is_modal column.
#' @param confidence Minimum value, shared by both filter steps: clade-level
#'   posterior_support and (post-establishment) tip-level membership_prob.
#' @return A tibble with the same columns as detect_introductions()'s output:
#'   tipname, deme, intro_clade_id, inferred_intro_date, last_sample_date,
#'   inferred_intro_source, inferred_intro_source_probability,
#'   intro_confidence_state, clade_size, intro_node.
.introductions_from_simmap <- function(clade_dwell, tip_membership, confidence) {
  established <- dplyr::filter(clade_dwell, .data$posterior_support >= confidence)

  if (nrow(clade_dwell) > 0L && nrow(established) == 0L) {
    warning(
      sprintf(
        paste0(
          "%d candidate clade(s) found in clade_dwell but none reach ",
          "confidence %.2f (max posterior_support available %.2f); ",
          "explode_tree() will return zero introductions. Lower ",
          "`confidence` to include them."
        ),
        nrow(clade_dwell),
        confidence,
        max(clade_dwell$posterior_support)
      ),
      call. = FALSE
    )
  }

  if (nrow(established) == 0L) {
    return(tibble::tibble(
      tipname = character(0),
      deme = character(0),
      intro_clade_id = integer(0),
      inferred_intro_date = as.Date(character(0)),
      last_sample_date = as.Date(character(0)),
      inferred_intro_source = character(0),
      inferred_intro_source_probability = numeric(0),
      intro_confidence_state = numeric(0),
      clade_size = integer(0),
      intro_node = integer(0)
    ))
  }

  # Restrict tip candidates to established clades BEFORE ranking -- a tip's
  # global argmax (raw is_modal) can point at a clade that didn't establish.
  candidates <- tip_membership |>
    dplyr::inner_join(
      dplyr::select(established, "intro_node", "deme"),
      by = c("intro_node", "deme")
    )

  # Recompute each tip's modal pick among survivors only.
  modal <- candidates |>
    dplyr::group_by(.data$tipname) |>
    dplyr::slice_max(.data$membership_prob, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$membership_prob >= confidence)

  modal |>
    dplyr::left_join(established, by = c("intro_node", "deme")) |>
    dplyr::rename(intro_confidence_state = "posterior_support") |>
    dplyr::select(
      "tipname",
      "deme",
      "intro_clade_id",
      "inferred_intro_date",
      "last_sample_date",
      "inferred_intro_source",
      "inferred_intro_source_probability",
      "intro_confidence_state",
      "clade_size",
      "intro_node"
    )
}
```

Then replace the Task 3 stub inside `explode_tree()`:
```r
  } else {
    stop("simmap source not yet implemented") # replaced in Task 4-6
  }
```
with:
```r
  } else {
    intros <- .introductions_from_simmap(clade_dwell, tip_membership, confidence)
  }
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add R/explode_tree.R tests/testthat/test-explode_tree.R
git commit -m "Implement .introductions_from_simmap() happy path

Produces the same $introductions column set as detect_introductions()
for the baseline fixture: established regionB clade (2 tips) + regionC
singleton, matching the node-marginal path's output exactly."
```

---

### Task 5: Establishment filter (`posterior_support`) + warning

**Files:**
- Modify: `tests/testthat/test-explode_tree.R`

**Interfaces:**
- Consumes: `.introductions_from_simmap()` (Task 4), `make_intro_clade_dwell()`, `make_intro_tip_membership()`.
- Produces: no production code change — this task is regression coverage proving the establishment filter already implemented in Task 4 behaves correctly on two edge cases the happy-path test didn't exercise.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-explode_tree.R`:

```r
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
```

- [ ] **Step 2: Run the tests to verify they fail or pass unexpectedly**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: both new tests already pass — Task 4's implementation already contains the establishment filter and warning logic. This task exists to lock that behavior in with explicit regression coverage, not to add new code.

- [ ] **Step 3: (No implementation step — Task 4 already covers this.) Confirm via Step 2's output that both tests are green.**

- [ ] **Step 4: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add tests/testthat/test-explode_tree.R
git commit -m "Add regression coverage for the posterior_support establishment filter

Confirms a low-support phantom clade is excluded even though it's
grounded in a real tree node, and that explode_tree() warns + returns
zero rows when nothing in clade_dwell clears confidence."
```

---

### Task 6: Tip-level floor + recompute-modal-among-survivors correctness

**Files:**
- Modify: `tests/testthat/test-explode_tree.R`

**Interfaces:**
- Consumes: `.introductions_from_simmap()` (Task 4).
- Produces: regression coverage for the two remaining edge cases the spec calls out explicitly: a tip dropped for low `membership_prob`, and a tip whose raw global `is_modal` points at a non-established clade but is correctly rescued via its established-clade candidate row.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-explode_tree.R`:

```r
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
```

- [ ] **Step 2: Run the tests to verify they pass**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "explode_tree", reporter = "summary")'`
Expected: both pass against Task 4's existing implementation — `dplyr::slice_max()` inside `.introductions_from_simmap()` already operates on `candidates` (post-establishment-filter), not the raw `tip_membership` table, so the recompute-among-survivors behavior is already correct. This task locks it in with explicit coverage for the specific failure mode described in the spec.

- [ ] **Step 3: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add tests/testthat/test-explode_tree.R
git commit -m "Add regression coverage for tip-level floor + modal recompute

Confirms a low-membership_prob tip is dropped even when its clade
establishes, and that a tip whose raw global is_modal points at a
non-established clade is still correctly assigned via its established
candidate row (not silently dropped)."
```

---

### Task 7: Full-suite regression check for `explode_tree()`

**Files:**
- None (verification-only task).

**Interfaces:**
- Consumes: everything from Tasks 1-6.
- Produces: confidence that the node-marginal path is untouched and the simmap path is fully wired end-to-end, before moving to the plotting half of this plan.

- [ ] **Step 1: Run the full test suite**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(reporter = "summary")'`
Expected: 0 failures across every test file, including `test-plot_explode_tree.R` (untouched by this task's changes — confirms the node-marginal path's downstream consumer still works).

- [ ] **Step 2: Regenerate docs**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::document()'`
Expected: `man/explode_tree.Rd` updates to reflect the new parameters; no errors.

- [ ] **Step 3: Commit the regenerated docs**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add man/explode_tree.Rd NAMESPACE
git commit -m "Regenerate docs for explode_tree()'s new dual-source parameters"
```

---

### Task 8: Extract `.order_clades()` from `plot_explode_tree()`

**Files:**
- Modify: `R/plot_explode_tree.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: internal `.order_clades(clade_tbl, min_clade = 1)` → tibble with `local_clade_num` added, filtered to `clade_size >= min_clade`, ordered by `(deme, inferred_intro_date)`. Used by `plot_explode_tree()` (refactored to call it) and, in Task 9, by `plot_clade_persistence()`.

- [ ] **Step 1: Confirm the current test suite is green before refactoring**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "plot_explode_tree", reporter = "summary")'`
Expected: all pass (this is the baseline the refactor must not change).

- [ ] **Step 2: Extract the helper**

In `R/plot_explode_tree.R`, find:
```r
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
      .persist_num = as.numeric(.data$last_sample_date - .data$inferred_intro_date),
      persistence = if (is.null(persistent_days)) {
        dplyr::if_else(.data$.persist_num > persistence_days, "persistent", "transient")
      } else {
        dplyr::case_when(
          .data$.persist_num <  persistence_days ~ "transient",
          .data$.persist_num >  persistent_days  ~ "persistent",
          .default                               = "semi-persistent"
        )
      }
    ) |>
    # Number clades 1..N within each deme, ordered chronologically by intro date
    dplyr::filter(.data$clade_size >= min_clade) |>
    dplyr::arrange(.data$deme, .data$inferred_intro_date) |>
    dplyr::mutate(
      local_clade_num = dplyr::row_number(),
      .by = "deme"
    )
```

Replace with:
```r
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
      .persist_num = as.numeric(.data$last_sample_date - .data$inferred_intro_date),
      persistence = if (is.null(persistent_days)) {
        dplyr::if_else(.data$.persist_num > persistence_days, "persistent", "transient")
      } else {
        dplyr::case_when(
          .data$.persist_num <  persistence_days ~ "transient",
          .data$.persist_num >  persistent_days  ~ "persistent",
          .default                               = "semi-persistent"
        )
      }
    ) |>
    .order_clades(min_clade)
```

Add the new helper function to the same file, above `plot_explode_tree()` (no `@export`):
```r
#' Filter, order, and number clades within each deme for plotting
#'
#' Shared row-ordering logic between plot_explode_tree() and
#' plot_clade_persistence() -- keeping both in one place guarantees their
#' rows line up when the two plots are viewed side by side.
#'
#' @param clade_tbl A tibble with (at least) clade_size, deme,
#'   inferred_intro_date columns, one row per clade.
#' @param min_clade Minimum clade_size to keep.
#' @return clade_tbl, filtered, arranged by (deme, inferred_intro_date), with
#'   a local_clade_num column added (1..N within each deme).
.order_clades <- function(clade_tbl, min_clade = 1) {
  clade_tbl |>
    dplyr::filter(.data$clade_size >= min_clade) |>
    dplyr::arrange(.data$deme, .data$inferred_intro_date) |>
    dplyr::mutate(
      local_clade_num = dplyr::row_number(),
      .by = "deme"
    )
}
```

- [ ] **Step 3: Run tests to verify no behavior change**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "plot_explode_tree", reporter = "summary")'`
Expected: identical pass/fail result as Step 1 — pure refactor, no behavior change.

- [ ] **Step 4: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add R/plot_explode_tree.R
git commit -m "Extract .order_clades() helper from plot_explode_tree()

Shared row-ordering (filter by min_clade, arrange by deme + intro date,
number local_clade_num) so plot_clade_persistence() (next task) can
guarantee its rows line up with plot_explode_tree()'s."
```

---

### Task 9: `plot_clade_persistence()` — implementation + basic tests

**Files:**
- Create: `R/plot_clade_persistence.R`
- Create: `tests/testthat/test-plot_clade_persistence.R`

**Interfaces:**
- Consumes: `.order_clades()` (Task 8), `make_intro_raw_clade_dwell()`, `make_intro_clade_dwell()` (Task 2).
- Produces: exported `plot_clade_persistence(raw_clade_dwell, clade_dwell, min_clade = 1)` → `ggplot` object.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-plot_clade_persistence.R`:

```r
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "plot_clade_persistence", reporter = "summary")'`
Expected: FAIL — `plot_clade_persistence` not found.

- [ ] **Step 3: Implement `plot_clade_persistence()`**

Create `R/plot_clade_persistence.R`:

```r
#' Plot per-clade persistence_days distributions from simmap replicates
#'
#' Companion to plot_explode_tree(): one horizontal box-and-whisker per
#' introduction clade, showing the distribution of persistence_days across
#' the simmap posterior (the n_sim make.simmap() replicates). persistence_days
#' is a duration, not a calendar date, so it can't share plot_explode_tree()'s
#' x-axis -- this is a separate plot, not a layer added to that one. Row
#' ordering matches plot_explode_tree() exactly (via the shared
#' .order_clades() helper) so the two plots align when viewed side by side.
#'
#' @param raw_clade_dwell Tibble, one row per (sim_id, intro_node): intro_node,
#'   deme, persistence_days (the shape of dta_compute.R's
#'   _dta_clade_dwell_raw.tsv). Values come from here only -- never from
#'   clade_dwell -- so the plot reflects the actual posterior, not a summary
#'   of a summary.
#' @param clade_dwell Tibble, one row per clade: intro_clade_id, intro_node,
#'   deme, inferred_intro_date, clade_size (the shape of
#'   _dta_clade_dwell.tsv). Used only for row ordering and the min_clade
#'   filter.
#' @param min_clade Minimum clade_size to plot. Default 1.
#'
#' @return A ggplot object.
#' @export
plot_clade_persistence <- function(raw_clade_dwell, clade_dwell, min_clade = 1) {
  required_raw <- c("intro_node", "deme", "persistence_days")
  missing_raw <- setdiff(required_raw, names(raw_clade_dwell))
  if (length(missing_raw) > 0L) {
    stop(
      "raw_clade_dwell is missing required columns: ",
      paste(missing_raw, collapse = ", ")
    )
  }

  required_clade <- c("intro_clade_id", "intro_node", "deme", "inferred_intro_date", "clade_size")
  missing_clade <- setdiff(required_clade, names(clade_dwell))
  if (length(missing_clade) > 0L) {
    stop(
      "clade_dwell is missing required columns: ",
      paste(missing_clade, collapse = ", ")
    )
  }

  ordered <- .order_clades(clade_dwell, min_clade)

  plot_df <- raw_clade_dwell |>
    dplyr::inner_join(
      dplyr::select(ordered, "intro_node", "deme", "intro_clade_id", "local_clade_num"),
      by = c("intro_node", "deme")
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$persistence_days,
      y = .data$local_clade_num,
      group = .data$intro_clade_id
    )
  ) +
    ggplot2::geom_boxplot(
      orientation = "y",
      fill = "grey70",
      color = "grey40"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$deme), scales = "free_y") +
    ggplot2::scale_y_continuous(breaks = function(x) {
      seq(ceiling(x[1]), floor(x[2]))
    }) +
    ggplot2::labs(x = "Persistence (days)", y = "Introduction") +
    ggplot2::theme_classic()
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "plot_clade_persistence", reporter = "summary")'`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add R/plot_clade_persistence.R tests/testthat/test-plot_clade_persistence.R
git commit -m "Add plot_clade_persistence(): boxplot of simmap persistence_days per clade

Companion to plot_explode_tree() -- persistence_days is a duration and
can't share that function's calendar-date x-axis, so this is a
separate function sharing row ordering via .order_clades()."
```

---

### Task 10: Row-alignment test vs `plot_explode_tree()`

**Files:**
- Modify: `tests/testthat/test-plot_clade_persistence.R`

**Interfaces:**
- Consumes: `plot_explode_tree()`, `plot_clade_persistence()`, `.order_clades()`.
- Produces: regression coverage guarding the shared `.order_clades()` helper against future drift between the two functions.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-plot_clade_persistence.R`:

```r
test_that("plot_clade_persistence's row ordering matches plot_explode_tree()'s for the same clade_dwell", {
  tree <- make_intro_tree()
  clade_dwell <- make_intro_clade_dwell(tree)

  # plot_explode_tree() accepts clade_dwell directly (its required columns
  # are a subset of clade_dwell's, per the tiphydr-compatibility design in
  # docs/superpowers/plans/isnt-the-dwell-time-eager-cocoa.md in the calling
  # pipeline repo).
  expl_built <- ggplot2::ggplot_build(plot_explode_tree(clade_dwell))
  pers_built <- ggplot2::ggplot_build(
    plot_clade_persistence(make_intro_raw_clade_dwell(tree), clade_dwell)
  )

  expl_order <- .order_clades(
    dplyr::distinct(
      clade_dwell,
      .data$intro_clade_id, .data$deme, .data$inferred_intro_date, .data$clade_size
    )
  )
  pers_order <- .order_clades(clade_dwell)

  expect_equal(
    dplyr::arrange(expl_order, .data$intro_clade_id)$local_clade_num,
    dplyr::arrange(pers_order, .data$intro_clade_id)$local_clade_num
  )
})
```

- [ ] **Step 2: Run the test to verify it fails or passes**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(filter = "plot_clade_persistence", reporter = "summary")'`
Expected: passes immediately — both functions call the same `.order_clades()` on equivalent inputs by construction (Task 8/9). This test exists as a permanent regression guard, not to drive new implementation.

- [ ] **Step 3: Commit**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add tests/testthat/test-plot_clade_persistence.R
git commit -m "Add regression test guarding plot_clade_persistence()/plot_explode_tree() row alignment"
```

---

### Task 11: Final full-suite check, docs, and manual verification

**Files:**
- None (verification-only task).

**Interfaces:**
- Consumes: everything from Tasks 1-10.
- Produces: a fully green package, regenerated `NAMESPACE`/`man/`, and one manual check against real pipeline output confirming the new code handles production-scale data.

- [ ] **Step 1: Run the full test suite**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::test(reporter = "summary")'`
Expected: 0 failures across all test files.

- [ ] **Step 2: Regenerate docs**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::document()'`
Expected: `man/plot_clade_persistence.Rd` created, `NAMESPACE` updated to export `plot_clade_persistence`; no errors.

- [ ] **Step 3: `R CMD check` style validation**

Run: `cd /Users/user/Programming_Directory/tiphydr && Rscript -e 'devtools::check(document = FALSE, quiet = TRUE)' 2>&1 | tail -60`
Expected: no new `ERROR`/`WARNING` introduced by this plan's changes (pre-existing NOTEs from unrelated in-progress work in the repo, per the uncommitted `CLAUDE.md`/`README.md`/`man/*.Rd` changes noted earlier in this session, are not this plan's concern).

- [ ] **Step 4: Manual check against real pipeline output**

Run:
```bash
cd /Users/user/Programming_Directory/tiphydr && Rscript -e '
devtools::load_all(".", quiet = TRUE)
clade_dwell <- readr::read_tsv(
  "/Users/user/Programming_Directory/Ebel_Lab/wnv-ge_dta_stan/runs/R028/R028_flyway_1_dta_clade_dwell.tsv",
  show_col_types = FALSE
)
tip_membership <- readr::read_tsv(
  "/Users/user/Programming_Directory/Ebel_Lab/wnv-ge_dta_stan/runs/R028/R028_flyway_1_dta_tip_membership.tsv",
  show_col_types = FALSE
)
tree <- ape::read.nexus(
  "/Users/user/Programming_Directory/Ebel_Lab/wnv-ge_dta_stan/runs/R028/R028_flyway_1_dated_dayadj.timetree.nex"
)
res <- explode_tree(tree, clade_dwell = clade_dwell, tip_membership = tip_membership)
cat("introductions:", nrow(res$introductions), "rows\n")
cat("trees:", length(res$trees), "subtrees\n")
cat("columns match detect_introductions() shape:",
    setequal(
      names(res$introductions),
      c("tipname", "deme", "intro_clade_id", "inferred_intro_date", "last_sample_date",
        "inferred_intro_source", "inferred_intro_source_probability",
        "intro_confidence_state", "clade_size", "intro_node")
    ), "\n")

raw_clade_dwell <- readr::read_tsv(
  "/Users/user/Programming_Directory/Ebel_Lab/wnv-ge_dta_stan/runs/R028/R028_flyway_1_dta_clade_dwell_raw.tsv",
  show_col_types = FALSE
)
p <- plot_clade_persistence(raw_clade_dwell, clade_dwell)
cat("plot_clade_persistence() built ok:", inherits(p, "ggplot"), "\n")
'
```
Expected: no errors; `introductions` and `trees` both non-empty; `columns match detect_introductions() shape: TRUE`; `plot_clade_persistence() built ok: TRUE`. This is a real 300+ tip run — confirms the plan's implementation holds up outside the small synthetic fixture.

- [ ] **Step 5: Commit the regenerated docs**

```bash
cd /Users/user/Programming_Directory/tiphydr
git add man/plot_clade_persistence.Rd NAMESPACE
git commit -m "Regenerate docs for plot_clade_persistence()"
```
