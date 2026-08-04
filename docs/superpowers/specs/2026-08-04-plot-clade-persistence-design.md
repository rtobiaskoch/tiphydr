# Design: `plot_clade_persistence()` — simmap-based persistence-distribution plot

## Context

`plot_explode_tree()` (see `docs/superpowers/plans/2026-06-23-plot_explode_tree.md` and `R/plot_explode_tree.R`) visualizes introduction clades as a timeline: one horizontal segment per clade from `inferred_intro_date` to `last_sample_date`, faceted by deme. It was built to accept either method's clade table — `detect_introductions()`'s node-marginal output, or (per the wnv-ge_dta_stan project's `docs/superpowers/plans/isnt-the-dwell-time-eager-cocoa.md`) the simmap-based `_dta_clade_dwell.tsv`, whose columns were deliberately aligned to `plot_explode_tree()`'s required set so it already accepts that table with zero changes.

That design doc flagged one deferred gap: rendering the simmap method's posterior uncertainty, which `explode_tree()`'s node-marginal output doesn't have. The originally-sketched form was HPD whiskers on the intro/last-sample points. In discussion, this was revised twice:

1. **Whiskers → violin.** A plain HPD line understates what the simmap posterior actually offers — the full shape of the distribution, not just an interval. `persistence_days` is the quantity of direct scientific interest (how long did this introduction actually persist, and how uncertain is that), so the violin is built for `persistence_days` specifically, not for the two endpoint dates.
2. **In-plot → companion plot.** `persistence_days` is a *duration* (units: days); `plot_explode_tree()`'s x-axis is a *calendar date*. The two cannot share one x-scale in one panel, which rules out adding this as a layer inside `plot_explode_tree()`. This is therefore a new, separate function, not a modification to the existing one.

**Why this belongs in `tiphydr`, not the calling pipeline (`wnv-ge_dta_stan`):** `tiphydr` is the shared analysis/plotting surface for both reconstruction methods (node-marginal and simmap-branch); `_dta_clade_dwell.tsv`/`_dta_clade_dwell_raw.tsv` were designed with tiphydr's column vocabulary specifically so downstream plotting could live here rather than as bespoke script code (contrast with `dta_plot.R`'s project-local `_deme_tree_dwell.pdf`, which draws simmap markers directly onto the ggtree phylogeny — a different, out-of-scope visualization not touched by this design).

## Goal

Add `plot_clade_persistence()`, a companion to `plot_explode_tree()`: one horizontal violin per introduction clade, showing the posterior distribution of `persistence_days` from the simmap replicates, with rows ordered/faceted identically to `plot_explode_tree()` so the two plots align when viewed side by side.

## API

```r
plot_clade_persistence(raw_clade_dwell, clade_dwell, min_clade = 1)
```

- `raw_clade_dwell` — per-`(sim_id, intro_node)` tibble, the shape of `_dta_clade_dwell_raw.tsv`. Required columns: `intro_node`, `deme`, `persistence_days`.
- `clade_dwell` — the summarized per-clade tibble, the shape of `_dta_clade_dwell.tsv`. Required columns: `intro_clade_id`, `intro_node`, `deme`, `inferred_intro_date`, `clade_size`. Used only to establish row order and apply the `min_clade` filter — never to source plotted values (those come from `raw_clade_dwell` only, so the violin reflects the actual posterior, not a summary of a summary).
- `min_clade` — same meaning and default (`1`) as `plot_explode_tree()`.

Returns a `ggplot` object.

## Row alignment with `plot_explode_tree()`

`plot_explode_tree()` currently inlines its row-ordering logic:
```r
dplyr::arrange(deme, inferred_intro_date) |>
  dplyr::mutate(local_clade_num = dplyr::row_number(), .by = "deme")
```
Extract this into an internal (unexported) helper, e.g. `.order_clades(clade_tbl, min_clade)`, and call it from both `plot_explode_tree()` and `plot_clade_persistence()`. This is the one shared-code touch point in `plot_explode_tree.R`; nothing else about that function changes.

`plot_clade_persistence()` joins `raw_clade_dwell` to the ordered clade table by `intro_node` (+ `deme` as a safety key) to attach each replicate row's `local_clade_num` and `intro_clade_id`.

## Plot spec

- One violin per clade: `geom_violin(aes(x = persistence_days, y = local_clade_num, group = intro_clade_id), orientation = "y")`.
- Single muted fill/color for every violin — `fill = "grey70", color = "grey40"` — not scaled by any variable (not `posterior_support`, not `deme`, not `persistence` class). Keeps this plot legible as a secondary/detail view next to `plot_explode_tree()`'s already-multi-encoded timeline (color = persistence class, fill = source deme, size = clade size).
- `facet_wrap(~deme, scales = "free_y")` — same faceting convention as `plot_explode_tree()`.
- Fixed row height/lane for every clade regardless of `posterior_support` (no violin-prominence scaling).
- x-axis label: `"Persistence (days)"`. y-axis: no label (matches `plot_explode_tree()`'s `y = "Introduction"` convention — reuse that label for consistency).

## Singleton-support fallback

A clade supported by exactly one simmap replicate has one `persistence_days` value — `geom_violin` cannot compute a kernel density from n = 1. For clades where `raw_clade_dwell` has `< 2` rows for that `intro_node`, draw a single `geom_point()` at that one value on the same row instead of a violin. This mirrors the existing `hpd_or_point()` fallback in `dta_compute.R` (`scripts/dta_compute.R`), which handles the identical n = 1 case for HPD intervals — same rationale, same codebase convention, now mirrored in the plotting layer.

## Testing

New `tests/testthat/test-plot_clade_persistence.R`, following `test-plot_explode_tree.R`'s existing pattern:

1. Returns a `ggplot` for valid input.
2. Errors on missing required columns (from either `raw_clade_dwell` or `clade_dwell`).
3. Produces one facet per destination deme.
4. A singleton-support clade (one row in `raw_clade_dwell` for some `intro_node`) results in a point layer being used for that clade's row rather than a violin — assert via inspecting the built plot's layer data, not just that the function runs.
5. Row alignment: given the same `clade_dwell` input, `local_clade_num` assignment matches what `plot_explode_tree()` would assign for the same clades (guards the shared `.order_clades()` helper against drift between the two functions).

New fixture helper in `tests/testthat/helper-fixtures.R`, alongside the existing `make_intro_*` generators: `make_intro_raw_clade_dwell()` (or similar name, matched to existing naming convention at implementation time), producing a small multi-replicate `raw_clade_dwell` tibble consistent with `make_intro_tree()`/`make_intro_node_probs()`'s existing fixture tree, plus at least one clade with only 1 supporting replicate to exercise the singleton fallback.

## Out of scope (explicitly not building)

- Any modification to `plot_explode_tree()`'s own rendering (no HPD whiskers added there — superseded by this violin approach for `persistence_days`; intro/last-sample date uncertainty is not visualized by this work).
- Generalizing `dta_plot.R`'s tree-overlay dwell plot (`_deme_tree_dwell.pdf`, mid-branch introduction markers on the ggtree phylogeny itself) into `tiphydr`. Different visualization, not addressed here.
- Scaling violin height/prominence by `posterior_support`.
- A combined single-figure/patchwork version merging the timeline and persistence plots into one image. `patchwork` is not added as a dependency; the two plots remain independently returned `ggplot` objects the caller can arrange as they choose.

## Files

- `R/plot_explode_tree.R` — extract `.order_clades()` helper; no other change.
- `R/plot_clade_persistence.R` — new file, exports `plot_clade_persistence()`.
- `NAMESPACE` / `man/plot_clade_persistence.Rd` — regenerated via `devtools::document()`.
- `tests/testthat/test-plot_clade_persistence.R` — new.
- `tests/testthat/helper-fixtures.R` — add one new fixture generator.

## Verification

1. Unit tests above pass via `devtools::test()`.
2. Manual check against a real run's outputs: call `plot_clade_persistence(readr::read_tsv("runs/R028/R028_flyway_1_dta_clade_dwell_raw.tsv"), readr::read_tsv("runs/R028/R028_flyway_1_dta_clade_dwell.tsv"))` and confirm the plot renders, row count per deme facet matches `plot_explode_tree()`'s output for the same `clade_dwell.tsv`, and at least one singleton-support clade (verified earlier in this project: several `intro_node`s in `R028_flyway_1_dta_clade_dwell.tsv` have `posterior_support` corresponding to 1/100 maps) renders as a point, not a violin.
