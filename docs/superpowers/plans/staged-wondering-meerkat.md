# Explicit cross-family combination rows for WNV lineage tables

## Context

`define_lineage()`'s matching algorithm is unchanged by this plan — it stays
strict (a lineage matches only if ALL its rows match) and picks the
candidate with the most matched mutations, exactly as it does today. That
algorithm redesign was explored and validated in an earlier round of this
session (dynamic `required`/optional-mutation matching), but the user opted
to defer it: it's more powerful but changes matching semantics, and for now
they want the current explicit, discrete lineage-row approach extended
instead. It's more verbose, but every match is a literal, auditable table
row rather than a computed label.

The gap being closed: the current 16-row `wnv_mut_gff`/`wnv_mut_orf` tables
already have single-family partial lineages (`SW03_NS4A`, `SW03_NS5`,
`NY10_NS2A`, `NY10_NS4B`) alongside the four base lineages (`NY99`, `WN02`,
`SW03`, `NY10`), but nothing for a sequence carrying markers from **both**
families at once — e.g. all of NY10's mutations plus SW03's `NS5314R` only.
Today such a sequence would match both `NY10` (3 mutations) and `SW03_NS5`
(2 mutations) and silently get labeled just `"NY10"`, losing the SW03
signal. The fix is to add explicit rows for these cross-family
combinations, named `"NY10-SW03_NS5"` style — reusing the existing
alphabetical-`sort()`-then-`"-"`-join convention the tie-break/hybrid logic
already uses elsewhere in `define_lineage()`, applied here as a literal
naming convention for new static rows rather than a computed tie-break.

## Combination space

Four independently-observable markers beyond the shared `env159A/449A`
clade marker: `NS4A85T`, `NS5314R` (SW03 family), `NS2A188K`, `NS4B240M`
(NY10 family). That's 16 presence/absence states; 7 are already covered
(`WN02` = all absent, plus the 6 existing single-family lineages/partials).
This plan adds the remaining 9 cross-family states. Per family, presence is
named via its existing single-family label (full name if both markers
present, `_NS4A`/`_NS5`/`_NS2A`/`_NS4B` suffix if only one); the two
per-family names are then sorted alphabetically and joined with `"-"`
(`"NY10"` < `"SW03..."` alphabetically, so NY10-side always comes first —
matches the user's own example, `NY10-SW03_NS5`).

| # | NS4A | NS5 | NS2A | NS4B | New lineage name | Rows (env + markers) |
|---|:-:|:-:|:-:|:-:|---|---|
| 1 | 1 | 0 | 0 | 1 | `NY10_NS4B-SW03_NS4A` | env, NS4A, NS4B (3) |
| 2 | 1 | 0 | 1 | 0 | `NY10_NS2A-SW03_NS4A` | env, NS4A, NS2A (3) |
| 3 | 1 | 0 | 1 | 1 | `NY10-SW03_NS4A` | env, NS4A, NS2A, NS4B (4) |
| 4 | 0 | 1 | 0 | 1 | `NY10_NS4B-SW03_NS5` | env, NS5, NS4B (3) |
| 5 | 0 | 1 | 1 | 0 | `NY10_NS2A-SW03_NS5` | env, NS5, NS2A (3) |
| 6 | 0 | 1 | 1 | 1 | `NY10-SW03_NS5` | env, NS5, NS2A, NS4B (4) — user's example |
| 7 | 1 | 1 | 0 | 1 | `SW03-NY10_NS4B` | env, NS4A, NS5, NS4B (4) |
| 8 | 1 | 1 | 1 | 0 | `SW03-NY10_NS2A` | env, NS4A, NS5, NS2A (4) |
| 9 | 1 | 1 | 1 | 1 | `NY10-SW03` | env, NS4A, NS5, NS2A, NS4B (5) |

**Update (post-implementation):** naming was refined to order the two family
pieces by mutation count descending (more specific piece first), falling
back to alphabetical on a tie — rows 7 and 8 above are `SW03-NY10_NS4B` and
`SW03-NY10_NS2A` (SW03 has 3 required mutations vs. 2), not the alphabetical
`NY10_...-SW03` order. Rows 1–6 and 9 are unaffected since their two pieces
already tie in mutation count (alphabetical order coincides).

Each new lineage's row count (3–5) exceeds every subset lineage it would
otherwise tie with (`NY10`=3, `SW03`=3, and the four partials=2 each), and
no two of the 9 new lineages share an identical required-residue set, so
under the current strict/most-specific-wins algorithm every one of these
states resolves to a unique, unambiguous winner — no new ties, no changes
needed to `define_lineage()`'s tie-break logic.

Concrete mutation values (already established in the current table):
`env` = pos 159 (gff) / 449 (orf), residue A (the WN02-clade marker; `env`
V = `NY99`, untouched by this plan). `NS4A` = pos 85/2209, residue T.
`NS5` = pos 314/2842, residue R. `NS2A` = pos 188/1331, residue K. `NS4B`
= pos 240/2513, residue M.

## Changes

### 1. `data-raw/wnv_mutations.R` — new file

No generation script currently exists for these tables (confirmed via
earlier exploration — the `.rda`s were hand-built once, with no tracked
source). Create this script defining both full 49-row tables (16 existing
+ 33 new rows across the 9 combinations above) as tibbles, structured so
the combinatorial expansion is easy to read and verify against the table
above — e.g. build each of the 9 new lineages' row-blocks from the
existing single-family building blocks (reusing the same env/NS4A/NS5/
NS2A/NS4B position-residue pairs already implicit in the current 16 rows)
rather than retyping raw numbers nine times. End with
`usethis::use_data(wnv_mut_gff, wnv_mut_orf, overwrite = TRUE)`.

Run the script once to regenerate `data/wnv_mut_gff.rda` and
`data/wnv_mut_orf.rda` in place (49 rows each).

**Packaging mechanics:**
- Add `usethis` to `DESCRIPTION`'s `Suggests:` (not currently listed).
- Add `^data-raw$` to `.Rbuildignore` (not currently listed).

### 2. `R/data.R` — update roxygen docs

Row count changes from 16 to 49 for both `wnv_mut_gff` and `wnv_mut_orf`.
No column changes (schema is unchanged: `lineage`/`gene`/`pos`/`residue`
for gff, `lineage`/`pos`/`residue` for orf). Regenerate `man/*.Rd` via
`devtools::document()`.

### 3. `R/define_lineage.R` — no code change

Confirmed above: the existing strict all-required-match + most-specific-
wins + alphabetical-sort-tie-break logic handles all 9 new lineages
correctly with no modification. Docstring example already cites a
composite-tie case (`"NY10_NS2A-SW03_NS4A"`) consistent with this naming
convention, so no doc update needed here either.

### 4. Tests

No existing test hardcodes the real shipped lineage names (confirmed via
earlier exploration — `tests/testthat/test-define_lineage.R` uses only
synthetic fixtures). Add a small new test file or block exercising the
real package data directly against a few of the 9 new states — at minimum
the user's own example (`NY10-SW03_NS5`) plus the boundary case (`NY10-
SW03`, all 5 markers) and one that should NOT match a new lineage (a
sequence with only the original `NY10` markers, to confirm it still gets
plain `"NY10"` and doesn't accidentally trip a new row). Construct a
minimal synthetic alignment/ref (following the existing
`helper-fixtures.R` pattern) with the right residues at the `wnv_mut_orf`
positions (449, 2209, 2842, 1331, 2513) rather than depending on any real
WNV sequence data.

## Verification

1. `Rscript data-raw/wnv_mutations.R` — regenerate `.rda`s, then
   `load("data/wnv_mut_gff.rda"); nrow(wnv_mut_gff)` / same for orf —
   confirm 49 rows each.
2. `devtools::load_all(".")`, manually construct a synthetic sequence with
   residues matching state 6 (`NS5`, `NS2A`, `NS4B` present, `NS4A`
   absent) and run `define_lineage()` against `wnv_mut_orf` with
   `mut_type = "aa"` — confirm the result is exactly `"NY10-SW03_NS5"`
   with no warning (unique winner, not a tie).
3. `devtools::test()` — full suite passes, including the new test cases.
4. `devtools::document()` — confirm only `man/wnv_mut_gff.Rd` and
   `man/wnv_mut_orf.Rd` change (row-count text).
5. `devtools::check()` — confirm no new `R CMD check` notes/warnings from
   the `data-raw/` addition.

---

# DEFERRED — Dynamic composite lineage labeling (future implementation)

**Status: not scheduled.** The user chose to implement the explicit
cross-family combination rows above instead, for now, because it keeps
every match as a literal, auditable table row rather than a computed
label. This section preserves the fully-designed and validated algorithm
alternative in case the table's combinatorial growth (49 rows today, and
it grows again with every new diagnostic mutation) later makes the
explicit-enumeration approach too unwieldy. Revisit this design at that
point rather than re-deriving it — it has already been through a full
brainstorming session, an Explore-agent codebase audit, and a Plan-agent
mechanical validation (traced the full 2×2×2×2 combinatorial space by hand
and confirmed zero discrepancies against the current data).

## Context

`R/define_lineage.R` assigns lineage labels by requiring ALL of a
lineage's mutations to match. The shipped WNV data (`wnv_mut_gff`,
`wnv_mut_orf`) needs "partial match" cases — e.g. a sequence with SW03's
shared `env159A` marker plus only its `NS4A85T` mutation (not
`NS5314R`) — and today the only way to catch that is to hand-duplicate a
whole new lineage row-set (`SW03_NS4A`) in the lookup table. This is
combinatorial: N independently-optional defining mutations per lineage
need up to 2^N−1 hand-written combo rows, and every new diagnostic
mutation makes it worse (the explicit-enumeration plan above already
needs 33 new rows just for the current 4-marker cross-family case).

Git history (commit `15d048f`) shows the table used to be flat — one row
per (lineage, mutation), no combo rows — before combos were hand-added to
work around `define_lineage()`'s all-or-nothing matching. This design
reverts to that flat structure, plus one new column, and fixes the
matching algorithm instead of continuing to hand-expand the data.

**Design, confirmed with user + validated by a Plan-agent mechanical
review:** each lineage row is tagged `required` (TRUE = "core" mutation
shared with parent clade, must always match) or not (FALSE = "defining"
mutation, independently optional). A sequence matching some-but-not-all
of a lineage's defining mutations gets a dynamically-built composite
label (`SW03_NS4A`), reusing the *existing* tie-break/hybrid-composite
logic (`"-"`-joined names on equal specificity) rather than replacing it.
The Plan-agent traced the full 2×2×2×2 combinatorial space of the
SW03/NY10 families and confirmed the new 8-row tables reproduce identical
outcomes to today's 16-row tables in every case.

Backward compatibility: if a `muts` table has no `required` column (the
existing generic test suite doesn't), every row defaults to `required =
TRUE`, which collapses the new logic to today's exact behavior. No
existing test needs to change for this reason alone.

## Changes

### 1. `R/define_lineage.R` — matching algorithm rewrite

Replace the `purrr::keep(...)` matching block and the `split()` call
(currently drops all columns but `pos`/`residue`) with:

- Before splitting: default-fill `muts$required <- TRUE` and
  `muts$gene <- as.character(muts$pos)` for any missing column. Add a
  one-time (not per-sequence) validation warning if any lineage has zero
  `required` rows (vacuous required-gate footgun — `all(logical(0))` is
  `TRUE` in R).
- `lineage_muts <- split(muts[, c("pos","residue","required","gene")], muts$lineage)`
- New `match_lineage(lin_name, seq_idx)` helper (extracted, replaces the
  inline closure) implementing the gate/score/label logic below. Reuses
  the existing `get_residue()` closure unchanged.
- New `row_matches(seq_idx, rows)` helper factoring out the shared
  `map2_lgl` position/residue comparison, used both for the required gate
  and the optional tally.

Per-lineage logic (validated pseudocode from Plan-agent review):
- `req_rows` (required==TRUE) must ALL match, else lineage excluded
  (`NULL`).
- If `nrow(opt_rows) == 0`: label = lineage name, score = `nrow(req_rows)`.
- Else count `n_opt_matched` and `unique()` the matched rows' `gene` tags:
  - all optional rows matched → label = lineage name (bare, no suffix),
    score = `nrow(req_rows) + n_opt_matched`.
  - zero optional rows matched → excluded (`NULL`) — this is what lets a
    less-specific base lineage like `WN02` win instead of tying with an
    empty-defining-set "SW03".
  - partial → label = `paste0(lineage_name, "_", paste(sort(unique(tags)), collapse = "_"))`,
    same score formula.

In `assign_one()`: build `matched <- purrr::compact(purrr::map(names(lineage_muts), match_lineage, seq_idx))`,
extract scores/labels, `unique()` the candidate label vector before the
existing `max()` + tie/hybrid-join block (guards against two distinct
lineages coincidentally computing the same composite string). The
ambiguous/unknown fallback block (iterates `unique(muts$pos)` across all
rows) is unaffected and stays as-is.

Refined `assign_one()` structure validated by Plan-agent review:

```r
# ── one-time setup, before assign_one is defined ──
if (!"required" %in% names(muts)) muts$required <- TRUE
if (!"gene"     %in% names(muts)) muts$gene     <- as.character(muts$pos)
# (optional new validation) warn/error if any lineage has zero required rows:
#   no_req <- names(which(!vapply(split(muts$required, muts$lineage), any, logical(1))))
#   if (length(no_req)) warning("Lineage(s) with no required rows: ", paste(no_req, collapse=", "))

lineage_muts <- split(muts[, c("pos", "residue", "required", "gene")], muts$lineage)

# shared row-matching helper (replaces the inline map2_lgl closure)
row_matches <- function(seq_idx, rows) {
  purrr::map2_lgl(rows$pos, rows$residue, function(p, r) {
    val <- get_residue(seq_idx, p)
    !is.na(val) && val == toupper(r)
  })
}

# Returns NULL (not a candidate) or list(label=, score=)
match_lineage <- function(lin_name, seq_idx) {
  rows     <- lineage_muts[[lin_name]]
  req_rows <- rows[rows$required, , drop = FALSE]
  opt_rows <- rows[!rows$required, , drop = FALSE]

  if (nrow(req_rows) > 0L && !all(row_matches(seq_idx, req_rows))) {
    return(NULL)                                   # required gate fails
  }
  if (nrow(opt_rows) == 0L) {
    return(list(label = lin_name, score = nrow(req_rows)))
  }

  opt_hit       <- row_matches(seq_idx, opt_rows)
  n_opt_matched <- sum(opt_hit)

  if (n_opt_matched == nrow(opt_rows)) {
    return(list(label = lin_name, score = nrow(req_rows) + n_opt_matched))
  }
  if (n_opt_matched == 0L) {
    return(NULL)                                   # optional-only, nothing present
  }

  tag <- paste(sort(unique(opt_rows$gene[opt_hit])), collapse = "_")
  list(
    label = paste0(lin_name, "_", tag),
    score = nrow(req_rows) + n_opt_matched
  )
}

assign_one <- function(seq_idx) {
  matched <- purrr::compact(
    purrr::map(names(lineage_muts), match_lineage, seq_idx = seq_idx)
  )

  if (length(matched) == 0L) {
    # unchanged ambiguous/unknown fallback (still iterates unique(muts$pos))
    ...
  }

  scores     <- purrr::map_dbl(matched, "score")
  max_score  <- max(scores)
  candidates <- unique(purrr::map_chr(matched, "label")[scores == max_score])  # unique() guards collision edge case

  if (length(candidates) > 1L) {
    hybrid_label <- paste(sort(candidates), collapse = "-")
    warning(
      "Sequence '", seq_names[seq_idx], "' matched lineages with equal specificity: ",
      paste(sort(candidates), collapse = ", "), ". Assigning composite lineage '", hybrid_label, "'."
    )
    return(hybrid_label)
  }
  candidates
}
```

This preserves every existing helper (`get_residue`, the `lineage_muts`
split-by-name pattern, the `.iupac_ambiguous`/fuzzy-codon fallback block
untouched) and only changes the shape of what gets split and how
`matched` is scored/labeled.

Update the function's roxygen docs: describe the new optional `required`
and `gene` columns, backward-compat default, and how composite labels are
now generated dynamically instead of requiring pre-enumerated combo rows.

### 2. `R/resolve_mut_positions.R` — pass through `required`

Currently hardcodes exactly 4 input columns + 3 computed ones in its
output `tibble::tibble(...)` (line ~163), silently dropping anything
else. Add `required = if ("required" %in% names(muts)) muts$required else TRUE`
to that output, inserted right after `residue` (mirrors input-column
order, computed columns last): final column order `lineage, gene, pos,
residue, required, comparison_type, codon_start, aa_start`. `gene`
already passes through unchanged — no change needed there.

`expect_named()` in `tests/testthat/test-resolve_mut_positions.R` defaults
to `ignore.order = FALSE`, so the expected-column vector must be updated
to this exact new order, not just have `required` added anywhere.

### 3. `data-raw/wnv_mutations.R` — canonical (de-duplicated) tables

Proposed canonical replacement content (validated against the current
16-row tables for every biologically meaningful input — Plan-agent traced
the full 2×2×2×2 SW03/NY10 combinatorial space and found zero
discrepancies):

```
wnv_mut_gff (8 rows):
lineage, gene, pos, residue, required
NY99,  env,  159, V, TRUE
WN02,  env,  159, A, TRUE
SW03,  env,  159, A, TRUE
SW03,  NS4A,  85, T, FALSE
SW03,  NS5,  314, R, FALSE
NY10,  env,  159, A, TRUE
NY10,  NS2A, 188, K, FALSE
NY10,  NS4B, 240, M, FALSE

wnv_mut_orf (8 rows, NEW gene column added):
lineage, pos, residue, required, gene
NY99, 449, V, TRUE, env
WN02, 449, A, TRUE, env
SW03, 449, A, TRUE, env
SW03, 2209, T, FALSE, NS4A
SW03, 2842, R, FALSE, NS5
NY10, 449, A, TRUE, env
NY10, 1331, K, FALSE, NS2A
NY10, 2513, M, FALSE, NS4B
```

This is a much smaller table (8 rows vs. 49 in the explicit-enumeration
plan above) because the dynamic algorithm covers all 16 presence/absence
states from just these 8 canonical rows — no combinatorial expansion
needed in the data itself.

Same packaging mechanics apply as the explicit-enumeration plan: add
`usethis` to `Suggests:`, add `^data-raw$` to `.Rbuildignore`.

### 4. `R/data.R` — update roxygen docs

Row count changes from (whatever it is at the time this is picked back
up) to 8; document the new `required` (logical) and, for `wnv_mut_orf`,
the new `gene` (chr) column. Regenerate `man/*.Rd` via
`devtools::document()`.

### 5. Tests

- `tests/testthat/test-resolve_mut_positions.R`: update the
  `expect_named()` assertion to the new 8-column order.
- `tests/testthat/test-define_lineage.R`: add new test cases (existing
  ones are untouched — they don't use a `required` column, so they
  exercise the backward-compat default path). New fixture mirroring the
  SW03-style hierarchy (one required "core" row + two optional "defining"
  rows tagged by gene), covering:
  - zero defining mutations present → falls through to a separate
    base-only lineage (not the family name, not "ambiguous"/"unknown").
  - one defining mutation present → dynamic single-tag composite label.
  - both defining mutations present → collapses to the bare lineage name
    (no suffix).
  - a tie between two dynamically-labeled partial matches still produces
    the existing sorted `"-"`-joined hybrid label with a warning.
  - a lineage with no `gene` column falls back to `pos`-based tags.

## Known edge cases (flagged by Plan-agent review — must be handled, not just noted)

- **(a) Duplicate/colliding dynamic labels.** Unlike today, where matched
  names come from unique `split()` keys, the new scheme computes label
  strings at runtime — two distinct named lineages could coincidentally
  compute the same string. `unique()` the candidate label vector before
  the tie-break `sort()`/`paste()` join (already reflected in the
  pseudocode above). Consider also validating at input time that no
  lineage name collides with `{other_lineage}_{gene}` for any gene tag in
  the table.
- **(b) Duplicate gene tags within one lineage's optional rows** — two
  optional rows sharing a `gene` value would double up in the label
  (`"SW03_NS4A_NS4A"`) without the `unique()` already applied to
  `opt_present_tags` in the pseudocode above.
- **(c) A lineage with zero required rows** — `all(logical(0))` is `TRUE`
  in R, so a lineage with no `required` rows becomes a candidate for
  almost any sequence sharing even one unrelated marker. Add the one-time
  validation warning noted in the setup step above.
- **(d) Sandbox scripts** — `sandbox/testing ds_1.R` hardcodes
  `muts[-4, ]` (a row-index-based drop) and uses a stale `gff=` argument
  name; would need re-checking against whatever the canonical table looks
  like when this is implemented. `sandbox/ds_0.25_check.R` consumes
  `wnv_mut_orf` wholesale with no hardcoded assumptions — should keep
  working unmodified.

## Verification (when implemented)

1. `devtools::load_all(".")` then run the specific reproduction cases the
   Plan-agent traced by hand: a sequence with only `env159A` → `"WN02"`;
   one with `env159A + NS4A85T` only → `"SW03_NS4A"`; one with all three
   SW03 mutations → `"SW03"`; the existing NY10/SW03 tie case from the
   function's own docstring → `"NY10_NS2A-SW03_NS4A"` with a warning.
2. `devtools::test()` — full suite must pass, including the unmodified
   `test-define_lineage.R` cases (confirms backward-compat default path).
3. `devtools::document()` to regenerate `man/*.Rd` — diff to confirm only
   the two data docs and `define_lineage()`'s Rd changed.
4. `Rscript data-raw/wnv_mutations.R` — regenerate `.rda`s, confirm 8 rows
   each, matching the tables above exactly.
5. `devtools::check()` — confirm no new `R CMD check` notes/warnings.
