# Design: `explode_tree()` dual-source (node-marginal + simmap) + `plot_clade_persistence()`

## Context

**Two existing, independent ways to "explode" a dated tree into introduction clades:**

**Node-marginal** (`tiphydr::explode_tree()` / `detect_introductions()`, see `R/explode_tree.R` and `R/detect_introductions.R`): takes `node_probs` — one row per tree node, one marginal state probability per trait state, from a `fitMk` ancestral-state reconstruction. Collapses each node to one `inferred_state` (argmax) + `confidence_state`. Flags an edge as an introduction wherever child and parent `inferred_state` differ; an internal-node transition only counts as **established** if `confidence_state >= confidence` (default 0.5) — below that, only tip-level singleton introductions survive, with a warning. Established clades get every descendant tip whose path back to the introduction node never leaves the destination state ("continued transmission"). Output is one fixed, deterministic clade assignment per tip — no posterior, no mid-branch timing.

**Simmap-branch** (`clade_dwell_from_simmap()` in the calling pipeline's `scripts/dta_compute.R:80-197`, wnv-ge_dta_stan): built from the actual `$maps` branch-segment histories of `n_sim` (100) `phytools::make.simmap()` replicates — each one a full coherent joint history, not a node snapshot. Every segment boundary is an introduction, timed exactly where it falls mid-branch. Aggregated across replicates by `(intro_node, deme)` (comparable because every replicate shares the same fixed topology) into `_dta_clade_dwell.tsv` (one row per clade: `posterior_support` = fraction of replicates in which the introduction occurs at all, median + 95% HPD timing/persistence, modal source) and `_dta_tip_membership.tsv` (one row per `(tipname, intro_node)`: `membership_prob` = fraction of replicates in which *that tip* was swept into *that clade*, `is_modal` = that tip's single highest-`membership_prob` row). Gains mid-branch precision and a real posterior; costs a different output grain — per-clade-with-uncertainty and per-tip-*probability*, not per-tip-single-assignment.

`_dta_clade_dwell.tsv`'s columns were deliberately aligned to `plot_explode_tree()`'s required set (see `docs/superpowers/plans/isnt-the-dwell-time-eager-cocoa.md` in wnv-ge_dta_stan) so the existing timeline plot already accepts it with zero changes. That design doc flagged the reconciliation between the two methods as a deferred gap. This spec closes it in `tiphydr`, the shared analysis surface for both methods, rather than in the calling pipeline.

## Part A — `explode_tree()`: accept either source

### Goal

One function, one return contract (`list(introductions, trees)`), two mutually-exclusive input sources. `plot_explode_tree()` and the calling pipeline's `dta_explode.R` need zero changes regardless of which source produced the result.

### API

```r
explode_tree(
  tree,
  node_probs = NULL,           # existing path
  annotated_tree_path = NULL,  # existing path only (node date lookup)
  clade_dwell = NULL,          # new path
  tip_membership = NULL,       # new path
  confidence = 0.5
)
```

Dispatch by presence, not a `method=` flag: exactly one of `node_probs` or `{clade_dwell, tip_membership}` (both required together) must be supplied. Any other combination (neither, both, or one of the simmap pair missing) errors immediately, naming the valid combinations. `annotated_tree_path` is only used by the node path (`build_tree_df()`'s node-date lookup) — the simmap tables already carry `inferred_intro_date`/`last_sample_date`.

`tip_membership` must be the **full** per-`(tipname, intro_node, deme)` table (every row, not pre-filtered to `is_modal == TRUE`) — the simmap path recomputes each tip's modal pick after establishment filtering (below), so it needs every candidate row to choose from, not just the globally-precomputed modal flag.

### Why `confidence` gates two different things here, on purpose

`confidence` in the node path answers one question: is this clade real (`confidence_state >= confidence`)? The simmap path has **two** distinct uncertainty axes with no node-path analogue collapsing them into one:

- **Clade-level**: is this introduction event real at all? → `posterior_support` (fraction of the 100 replicates in which it's recorded, regardless of which tips it reached). This is the structural analogue of `confidence_state` — same question, same role.
- **Tip-level**: is *this specific tip's* membership in that clade reliable? → `membership_prob` (fraction of replicates in which *this tip* specifically was swept in). The node path has no equivalent — it assigns tips to an established clade deterministically via continued-transmission descent, all-or-nothing, so there's never a separate per-tip confidence number to gate there.

Mathematically, `membership_prob(tip) <= posterior_support(its clade)` always — `clade_dwell_from_simmap()` only ever emits a row for `(intro_node, deme)` in a replicate when that replicate's reached-tip set is non-empty (`scripts/dta_compute.R:153-155`), so the replicates a tip's membership can draw from are a subset of the replicates the clade itself appears in. That bound means a `membership_prob` floor alone is not equivalent to a `posterior_support` floor: a clade can be robustly established (`posterior_support` high) while its tip-level membership is diffuse enough that *no single tip* clears a `membership_prob` floor on its own — with only a tip-level filter, that entire clade silently vanishes from output, even though the node path's equivalent situation (an established transition with uniformly low-confidence descendants) doesn't happen, because the node path has no per-tip uncertainty to begin with. Both filters are needed, applied in sequence, not as alternatives.

**Decision:** one shared `confidence` value gates both steps (default 0.5) rather than two independent parameters — matches the node path's single-parameter API for callers that only ever set one value; can split into a separate `tip_confidence` later if real runs show the two need to diverge.

### Simmap-path algorithm

New internal (unexported) helper, e.g. `.introductions_from_simmap(tree, clade_dwell, tip_membership, confidence)`:

1. **Establish clades**: `established <- dplyr::filter(clade_dwell, posterior_support >= confidence)`. If `clade_dwell` had rows but none clear `confidence`, warn — mirroring `detect_introductions()`'s existing warning for the identical no-established-transitions case (`R/detect_introductions.R`) — and proceed with whatever (possibly zero) rows survive.
2. **Restrict tip candidates to established clades**: inner-join/filter `tip_membership` down to only `(intro_node, deme)` pairs present in `established`, *before* ranking — a tip's rows for non-established clades are never candidates.
3. **Recompute each tip's modal pick among survivors**: per `tipname`, take the row with the highest `membership_prob` among its remaining (established-clade) rows — not the globally-precomputed `is_modal` flag from `_dta_tip_membership.tsv`, which was computed across *all* clades including non-established ones and so can point to a clade that didn't survive step 1.
4. **Apply the tip-level floor**: keep only tips whose recomputed-modal `membership_prob >= confidence`. Tips with no surviving candidate row, or whose best surviving candidate doesn't clear the floor, are dropped (never emitted) — matching the node path's existing behavior for background/never-introduced tips, which likewise never appear in `$introductions`.
5. **Join to `established`** for the remaining columns: `intro_clade_id`, `inferred_intro_date`, `last_sample_date`, `inferred_intro_source`, `inferred_intro_source_probability`, `clade_size`, and map `posterior_support -> intro_confidence_state` (same column name and conceptual role the node path uses). Output columns are identical to `detect_introductions()`'s: `tipname`, `deme`, `intro_clade_id`, `inferred_intro_date`, `last_sample_date`, `inferred_intro_source`, `inferred_intro_source_probability`, `intro_confidence_state`, `clade_size`, `intro_node`.

### Tree extraction — reused unchanged

`explode_tree()`'s existing bottom half (`R/explode_tree.R`, building `$trees` via `ape::extract.clade()` + `ape::keep.tip()`) only reads `intro_clade_id`/`intro_node`/`tipname`/`clade_size` from the `intros` tibble — it never touches `node_probs` directly. No changes needed there; both source paths feed it the same shape.

### Documented limitation (not solved here, flagged for the future)

Step 4's floor is tip-level only; there is currently no output signal for "this clade established (step 1) but zero tips individually cleared the tip-level floor (step 4)" beyond that clade simply contributing no rows to `$introductions`. A future addition could surface this case explicitly (e.g. a `dropped_clades` diagnostic), but it's not built here — YAGNI until a real run shows this happening often enough to matter.

## Part B — `plot_clade_persistence()`

### Goal

A companion to `plot_explode_tree()`: one horizontal box-and-whisker per introduction clade, showing the distribution of `persistence_days` across the simmap replicates. `persistence_days` is a *duration*; `plot_explode_tree()`'s x-axis is a *calendar date* — the two can't share one x-scale in one panel, so this is a new function, not a layer added to the existing one.

### API

```r
plot_clade_persistence(raw_clade_dwell, clade_dwell, min_clade = 1)
```

- `raw_clade_dwell` — per-`(sim_id, intro_node)` tibble, the shape of `_dta_clade_dwell_raw.tsv`. Required columns: `intro_node`, `deme`, `persistence_days`.
- `clade_dwell` — the summarized per-clade tibble, the shape of `_dta_clade_dwell.tsv`. Required columns: `intro_clade_id`, `intro_node`, `deme`, `inferred_intro_date`, `clade_size`. Used only for row ordering and the `min_clade` filter — plotted values come from `raw_clade_dwell` only.
- `min_clade` — same meaning/default as `plot_explode_tree()`.

Returns a `ggplot` object.

### Row alignment with `plot_explode_tree()`

Extract `plot_explode_tree()`'s inlined row-ordering logic (`dplyr::arrange(deme, inferred_intro_date) |> dplyr::mutate(local_clade_num = dplyr::row_number(), .by = "deme")`) into an internal helper, e.g. `.order_clades(clade_tbl, min_clade)`, shared by both functions — the one touch point in `plot_explode_tree.R`; nothing else about that function changes. `plot_clade_persistence()` joins `raw_clade_dwell` to the ordered table by `intro_node` (+ `deme`) to attach `local_clade_num`/`intro_clade_id` to every replicate row.

### Plot spec

- One box per clade: `geom_boxplot(aes(x = persistence_days, y = local_clade_num, group = intro_clade_id), orientation = "y")`.
- Single muted fill/color for every box — `fill = "grey70", color = "grey40"` — not scaled by any variable.
- `facet_wrap(~deme, scales = "free_y")` — same convention as `plot_explode_tree()`.
- Fixed row height for every clade regardless of `posterior_support`.
- x-axis label `"Persistence (days)"`; y-axis unlabeled (matches `plot_explode_tree()`'s `y = "Introduction"`).

No singleton-support special case is needed: `geom_boxplot()` on a single value degrades gracefully to a zero-width box (all quartiles equal that value), so clades with only one supporting replicate render sensibly with no extra code path.

### Testing

New `tests/testthat/test-plot_clade_persistence.R`, following `test-plot_explode_tree.R`'s existing pattern: returns a `ggplot`; errors on missing required columns; one facet per deme; row alignment matches `plot_explode_tree()`'s `local_clade_num` for the same `clade_dwell` input (guards the shared `.order_clades()` helper against drift).

New fixture helper in `tests/testthat/helper-fixtures.R`: `make_intro_raw_clade_dwell()` (name to be matched to convention at implementation time), a small multi-replicate `raw_clade_dwell` tibble consistent with the existing `make_intro_tree()`/`make_intro_node_probs()` fixture tree, including at least one clade with only one supporting replicate to confirm the degenerate-box case renders without error.

## Out of scope (explicitly not building)

- HPD whiskers on `plot_explode_tree()`'s own intro/last-sample points — superseded by the `persistence_days` boxplot approach.
- Generalizing `dta_plot.R`'s tree-overlay dwell plot (`_deme_tree_dwell.pdf`, mid-branch markers on the ggtree phylogeny) into `tiphydr`. Different visualization, not addressed here.
- Scaling box width/prominence by `posterior_support`.
- A combined single-figure/patchwork version merging the timeline and persistence plots. `patchwork` is not added as a dependency; the two plots remain independently returned `ggplot` objects.
- A `dropped_clades` diagnostic for clades that establish (step 1) but produce zero tip assignments (step 4) — flagged as a documented limitation above, not built now.
- Moving `clade_dwell_from_simmap()`'s traversal/aggregation logic from `dta_compute.R` into `tiphydr` — `explode_tree()`'s simmap path consumes already-aggregated tables only; `tiphydr` never touches a raw `make.simmap()` object.

## Files

- `R/explode_tree.R` — add `clade_dwell`/`tip_membership`/dual-source dispatch to `explode_tree()`; internal `.introductions_from_simmap()` helper.
- `R/plot_explode_tree.R` — extract `.order_clades()` helper; no other change.
- `R/plot_clade_persistence.R` — new file, exports `plot_clade_persistence()`.
- `NAMESPACE` / `man/plot_clade_persistence.Rd` — regenerated via `devtools::document()`. `man/explode_tree.Rd` updated for the new parameters.
- `tests/testthat/test-explode_tree.R` — extend with simmap-source cases (establishment filter, tip-level floor, warning on zero established clades, error on invalid input combinations).
- `tests/testthat/test-plot_clade_persistence.R` — new.
- `tests/testthat/helper-fixtures.R` — add fixture generator(s) for `clade_dwell`/`tip_membership`/`raw_clade_dwell` shapes.

## Verification

1. Unit tests above pass via `devtools::test()`.
2. Manual check against a real run: `explode_tree(tree, clade_dwell = readr::read_tsv(".../R028_flyway_1_dta_clade_dwell.tsv"), tip_membership = readr::read_tsv(".../R028_flyway_1_dta_tip_membership.tsv"))` and confirm `$introductions` has the same columns as the node-marginal path's output, `plot_explode_tree(result$introductions)` renders without modification, and no established clade with `posterior_support` near 1 is entirely absent from the output (spot-check against `_dta_clade_dwell.tsv`'s high-`posterior_support` rows).
3. `plot_clade_persistence()` manual check: call against the same run's `_dta_clade_dwell_raw.tsv` + `_dta_clade_dwell.tsv`, confirm row count per deme facet matches `plot_explode_tree()`'s output for the same `clade_dwell.tsv`, and that a known low-`posterior_support` clade (several exist in `R028_flyway_1_dta_clade_dwell.tsv` at 1/100 maps) renders as a degenerate box, not an error.
