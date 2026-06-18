# Plan: `explode_tree.R` — detect viral introductions from a dated tree + DTA node probabilities

## Context

The lab repeatedly parses phylogenetic "introductions" — points where an ancestral
lineage transitions from one deme into another — by hand in ad-hoc Rmd notebooks
(`examples/wnv-foco_phylo_intro_parse_v2.Rmd`). That workflow is base-R, loop-heavy,
hardcoded to a single deme ("CO"), and tied to TreeTime *mugration* output. We want a
reusable, tidyverse-style function in `R/` that works on the cleaner `test25` data
format (a dated `timetree.nwk` + a node×deme probability matrix), detects introductions
into **all** demes, applies the proven continued-transmission definition, and returns
both the split subtrees and a tidy summary table.

An introduction is an edge where the inferred deme of a node differs from its parent's
deme (a state transition toward a destination deme). Reintroductions are allowed: each
independent transition into a deme is its own clade. The continued-transmission filter
keeps a descendant tip only if *every* ancestor on the path back to the introduction
node is also the destination deme — i.e. the lineage never left that deme.

## Approach

A thin orchestrator (`explode_tree()`) reads the two files, parses tip metadata from the
tip labels, and composes four small pure helpers. Tree traversal uses `tidytree`
(`parent`, `ancestor`, `offspring`, `as_tibble`) — the package's stated phylo stack and
the same primitives the example relies on. Node dates come from the timetree itself
(no JSON), derived unit-agnostically. Output is a list: `$introductions` (tidy tibble,
one row per kept tip) and `$trees` (named list of `phylo` subtrees, one per multi-tip
clade, named by `intro_clade_id`).

Rationale for unit-agnostic dates: TreeTime branch lengths may be in days or decimal
years and the `test25` filename (`dayadj`) is ambiguous. Regressing tip calendar dates
(numeric) on root-to-node depth (`ape::node.depth.edgelength`) and predicting internal
nodes makes the slope (units) and intercept (root date) fall out of the data, so the
function is correct regardless of unit. This is a standard molecular-clock linear map
and is verified by checking predicted tip dates ≈ label dates.

## Steps

1. **`R/most_likely_deme.R`** — pure. Input: node-prob tibble (`node` + one column per
   deme). For each row pick the max-probability deme via `which.max` (no `tidyr`).
   Return `tibble(node, inferred_state, confidence_state)`. Ties → first deme (document).

2. **`R/node_dates_from_timetree.R`** — pure. Inputs: `phylo` tree, `tip_dates`
   (tibble `tip_label`, `date` as `Date`). Compute `ape::node.depth.edgelength(tree)`;
   fit `lm(as.numeric(date) ~ depth)` over tips; predict for every node; return
   `tibble(node, inferred_date)` (`Date`). Internal-node ids follow ape numbering
   (`Ntip + seq_len(Nnode)`), tip ids `seq_len(Ntip)`.

3. **`R/build_tree_df.R`** — assemble the annotated `tbl_tree`. Inputs: `phylo`,
   `node_probs`, `tip_meta` (tibble: `tip_label`, `deme`, `date`). Steps:
   - `tidytree::as_tibble(tree)` → base tree table (`node`, `parent`, `label`, ...).
   - Internal nodes: join `most_likely_deme()` output (state + confidence).
   - Tips: `inferred_state = deme` (from label), `confidence_state = 1`.
   - Join `node_dates_from_timetree()`; for tips override with observed label `date`.
   - Add `is_tip = node <= ape::Ntip(tree)`.
   Returns one tidy tibble with `node, parent, label, inferred_state,
   confidence_state, inferred_date, is_tip`.

4. **`R/detect_introductions.R`** — the core logic, pure, all demes, purrr style.
   Input: the `build_tree_df` tibble + `confidence` threshold (default 0.95).
   - **Internal-node intros:** internal nodes with `confidence_state >= confidence`
     whose `inferred_state != parent's inferred_state`. Destination = node state,
     source = parent state, source prob = parent confidence.
   - For each intro node: `tidytree::offspring(df, node, tiponly = TRUE)`; apply
     **continued-transmission filter** — keep tip iff all `ancestor()` rows between the
     tip and the intro node (`inferred_date >= intro_date`) are the destination deme
     (`map`/`keep`, replacing the example's `for` loops).
   - **Single-tip intros:** tips whose `inferred_state != parent's state`
     (`clade_size = 1`, no confidence filter — observed).
   - Combine; assign `intro_clade_id` by ranking distinct intro events on
     `inferred_intro_date` (chronological); compute `clade_size` per clade.
   - Return per-tip tibble: `tipname, deme, intro_clade_id, inferred_intro_date,
     inferred_intro_source, inferred_intro_source_probability, clade_size,
     intro_node` (ape id, kept internally for subtree building).

5. **`R/explode_tree.R`** — orchestrator. Inputs: `tree_path`, `node_probs_path`,
   `confidence = 0.95`, `delim = "|"`, `tip_fields` (default
   `c("strain","date","deme","county","lineage","host")`).
   - `ape::read.tree(tree_path)`; `readr`/`utils::read.delim` the probs tsv.
   - Parse tip labels into `tip_meta` via `stringr::str_split_fixed` (date = field 2,
     deme = field 3); coerce `date` to `Date`.
   - `build_tree_df()` → `detect_introductions()`.
   - For each `intro_clade_id` with `clade_size > 1`: `ape::extract.clade(tree,
     intro_node)` then `ape::keep.tip(subtree, kept_tips)` to enforce the
     continued-transmission membership. Name by `intro_clade_id`.
   - Return `list(introductions = <tibble, intro_node dropped>, trees = <named list>)`.

6. **`DESCRIPTION`** — add `tidytree` to Imports (justified: project's phylo stack;
   provides `as_tibble.phylo`, `parent`, `ancestor`, `offspring`). `ape` already present.

7. **Tests** (`tests/testthat/`) — one file per function, mirroring repo convention.
   Use a tiny hand-built `phylo` + known prob matrix where introductions are obvious,
   so expected clade membership/dates are checkable. Add an end-to-end test on the
   `test25` fixtures (copy into `tests/testthat/testdata/`).

8. **Docs** — roxygen on every function; update `docs/TODO.md` (move `explode_tree`
   to done) and `man/` via `devtools::document()`.

## Critical Files

| File | Action |
|------|--------|
| `R/most_likely_deme.R` | create — prob matrix → most likely deme + prob |
| `R/node_dates_from_timetree.R` | create — node calendar dates via clock regression |
| `R/build_tree_df.R` | create — assemble annotated tbl_tree |
| `R/detect_introductions.R` | create — transition + continued-transmission logic |
| `R/explode_tree.R` | create — orchestrator returning `list(introductions, trees)` |
| `DESCRIPTION` | edit — add `tidytree` to Imports |
| `tests/testthat/test-*.R` (5) | create — unit + e2e tests |
| `docs/TODO.md` | edit — mark `explode_tree` done |

## Reuse

- `R/tree_tip_swap.R:35` — pattern for `fixed = TRUE` matching / purrr over tips.
- `R/extract_metadata.R:73-104` — `stringr::str_split_fixed` delimited-parsing pattern;
  the tip-label parse in `explode_tree()` mirrors it (kept inline — labels are a
  character vector, not an `XStringSet`, so the function isn't directly callable;
  unifying both behind a shared `parse_delimited()` helper is a possible later cleanup).
- `tidytree::parent / ancestor / offspring` — replace the example's manual traversal.
- `ape::node.depth.edgelength / extract.clade / keep.tip / Ntip` — already a dependency.

## Verification

1. `devtools::load_all(".")` then:
   ```r
   res <- explode_tree(
     tree_path       = "temp/test25_dated_dayadj.timetree.nwk",
     node_probs_path = "temp/test25_dta_node_probs.tsv"
   )
   res$introductions   # tidy table, all demes, chronological intro_clade_id
   length(res$trees)   # one phylo per multi-tip clade
   ```
2. **Date sanity (unit check):** inside `node_dates_from_timetree`, predicted tip dates
   must be within a few days of the label dates — confirms the clock regression and that
   branch-length units cancelled correctly.
3. **Logic spot-check:** for one `intro_clade_id`, confirm by hand on the tree that the
   parent deme ≠ destination deme and that every kept tip's ancestors back to the intro
   node are the destination deme.
4. `devtools::document()` then `devtools::test()` — all new + existing tests green.
5. `devtools::check()` — no new NOTEs from the `tidytree` import.
