# Plan: `plot_tree_trait()` — Nextstrain-style trait tree visualization

## Context

tiphydr has a complete pipeline for phylodynamic analysis (`explode_tree()` →
`build_tree_df()` → `detect_introductions()`) but no visualization layer. Users
working with `build_tree_df()` output currently have no package-native way to
inspect trait distributions across the tree. This function closes that gap:
given the `tbl_tree` from `build_tree_df()` and a discrete trait column name,
it produces a nextstrain-style ggtree plot with ancestral-state branch coloring,
returning a ggplot2 object the user can further customize.

---

## Pipeline position

```
explode_tree() → build_tree_df() → plot_tree_trait()
```

---

## Function signature

```r
plot_tree_trait <- function(tree_df, trait, palette = "Dark2")
```

| Arg | Type | Required | Description |
|-----|------|----------|-------------|
| `tree_df` | `tbl_tree` | Yes | Output of `build_tree_df()` |
| `trait` | character | Yes | Column name in `tree_df` to color by (discrete only) |
| `palette` | character | No | RColorBrewer palette name; default `"Dark2"` |

Returns a ggplot2 object.

---

## Implementation

### File to create
`R/plot_tree_trait.R`

### Step 1 — Validate inputs

Check for sentinel columns from `build_tree_df()`: `node`, `parent`, `label`,
`branch.length`, `is_tip`. If any are absent:

```r
stop(
  "tree_df does not look like output from build_tree_df(). ",
  "To prepare your tree, run:\n",
  "  build_tree_df(tree, node_probs, metadata)"
)
```

Also check:
- `trait` column exists in `tree_df`
- `trait` column is character or factor (discrete only — numeric triggers `stop()`)

### Step 2 — Compute parent-trait for ancestral-state branch coloring

```r
parent_states <- dplyr::select(tree_df, node, parent_trait = !!rlang::sym(trait))
node_data     <- dplyr::left_join(tree_df, parent_states, by = c("parent" = "node"))
```

Root node gets `NA` for `parent_trait` — handled silently by
`scale_color_manual(na.translate = FALSE)`.

### Step 3 — Build palette

```r
trait_levels <- unique(stats::na.omit(node_data[[trait]]))
n_colors     <- max(3, length(trait_levels))   # brewer.pal minimum is 3
colors       <- RColorBrewer::brewer.pal(n_colors, palette)[seq_len(length(trait_levels))]
names(colors) <- trait_levels
```

### Step 4 — Build plot

```r
phylo <- tidytree::as.phylo(tree_df)

ggtree::ggtree(phylo) %<+% node_data +
  ggplot2::aes(color = parent_trait) +
  ggtree::geom_tippoint(ggplot2::aes(color = .data[[trait]])) +
  ggtree::geom_nodepoint(ggplot2::aes(color = .data[[trait]]), size = 2) +
  ggplot2::scale_color_manual(
    values       = colors,
    name         = trait,
    na.translate = FALSE
  ) +
  ggplot2::theme_classic()
```

- Branches colored by **parent node's trait** (nextstrain ancestral-state style)
- Tips colored by **own trait** via `geom_tippoint()`
- Internal nodes colored by **own trait** via `geom_nodepoint()`

---

## Files to modify

| File | Change |
|------|--------|
| `R/plot_tree_trait.R` | Create — new function |
| `DESCRIPTION` | Add to `Imports:` — `ggtree`, `ggplot2`, `RColorBrewer`, `rlang` |
| `NAMESPACE` | Auto-generated via `devtools::document()` |
| `man/plot_tree_trait.Rd` | Auto-generated via `devtools::document()` |

---

## Tests

File: `tests/testthat/test-plot_tree_trait.R`

Uses existing fixture from `tests/testthat/helper-fixtures.R` — no new test data needed.

| Test | What it checks |
|------|---------------|
| Returns ggplot object | Valid `tbl_tree` + `trait = "inferred_state"` → `inherits(result, "ggplot")` |
| Validation fires — wrong input | Plain tibble missing sentinel columns → `stop()` mentioning `build_tree_df()` |
| Validation fires — missing column | `trait` column name not in `tree_df` → informative `stop()` |
| Validation fires — numeric trait | Numeric column passed as `trait` → `stop()` explaining discrete-only |

---

## Verification

1. `devtools::load_all(".")` — confirm function loads without error
2. Run with a real `build_tree_df()` output and `trait = "inferred_state"` — visually inspect that branches are colored by parent state, tips and nodes by own state
3. `devtools::test()` — all 4 new tests pass
4. `devtools::check()` — no new warnings or errors
