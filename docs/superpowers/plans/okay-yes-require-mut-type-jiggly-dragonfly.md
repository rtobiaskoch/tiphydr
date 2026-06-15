# Simplify `define_lineage`: require `mut_type`, drop the gff/gene branch

## Context

`define_lineage()` currently carries three code paths: a `nuc` pathway, an `aa`
pathway, and a third `gff`/`gene`-column pathway that resolves gene-relative
amino-acid positions to absolute genomic coordinates *inside* the function
(via `resolve_mut_positions()`) and then does per-row nuc-vs-aa dispatch. This
third path duplicates ~80 lines of assignment logic (`get_residue_row`, a
second `assign_one`) and lets one table mix nuc and aa rows.

Two problems motivated this change:

1. **The original bug.** `mut_type` defaulted to a value, so a caller passing an
   amino-acid mutation table (`data/lineage_mutations_ORF.csv`, residues
   `A K M T V`) without setting `mut_type` silently compared AA letters against
   single nucleotides and returned all `unknown`. A content-based "is this aa or
   nuc?" guard cannot fix this: every residue in that table (`A`, `K`, `M`, `T`,
   `V`) is also a valid IUPAC nucleotide code, so the table is indistinguishable
   from a nuc table by inspection. The robust fix is to **make `mut_type`
   required** so intent is always explicit.

2. **Redundant complexity.** `resolve_mut_positions()` already emits an
   `aa_start` column — the 1-indexed amino-acid position from the CDS start —
   which is exactly the coordinate the `aa` pathway indexes into. So the gff
   path can move *out* of `define_lineage` and become an explicit preprocessing
   step the caller runs first. This deletes the duplicated branch and gives
   clean single-responsibility separation (resolve = coordinate conversion,
   define = assignment).

User decisions (this session): require `mut_type` (no default); mixed nuc+aa
tables are not used and may be dropped.

## Approach

### 1. `R/define_lineage.R` — remove gff branch, require `mut_type`

- Signature becomes `define_lineage(alignment, ref, muts, mut_type, verbose = TRUE)`
  — `mut_type` has **no default**; drop the `gff` parameter.
- Add an explicit guard: error with a clear message if `mut_type` is missing,
  and keep the existing `!mut_type %in% c("nuc","aa")` check (test at
  `test-define_lineage.R:54` expects `"type must be"`).
- Delete the `has_gene_col` detection block (`define_lineage.R:40-68`), keeping
  only the `required_cols <- c("lineage","pos","residue")` validation for `muts`.
- Delete the entire `has_gene_col == TRUE` branch: `get_residue_row` and the
  first `assign_one` (`define_lineage.R:135-214`). Keep only the existing
  nuc/aa `get_residue` closure (`define_lineage.R:91-125`) and the second
  `assign_one` (`define_lineage.R:215-274`).
- Update roxygen: drop the `@param gff`, drop the gene-column alternative in the
  `@param muts` description, and reword `@param mut_type` to state it is required.

### 2. `R/resolve_mut_positions.R` — docs only

No code change. Update the `@description` (`resolve_mut_positions.R:5-6`) so it no
longer says it is "consumed by `define_lineage()` when its `gff` argument is
supplied." It is now a standalone preprocessing step: its `aa_start` output is
fed to `define_lineage(mut_type = "aa")` as the `pos` column. Reuse the existing
`cds_start` argument (`resolve_mut_positions.R:46`, formula at lines 149-157).

### 3. `tests/testthat/test-define_lineage.R`

- Add `mut_type = "nuc"` to every call that currently relies on the removed
  default: lines 3, 13, 20, 30, 41, 49, 80, 272, 284, 301, 339 (the nuc-default
  and ambiguous-warning tests).
- Remove the gff-pathway test block (`test-define_lineage.R:66-183`): the gff
  regression, strict-match, single-mutation, wildtype, "mixed gene + nuc",
  "gene column without gff errors", and "gene column all-nuc fallthrough" tests.
  All depend on `define_lineage`'s removed gff/gene handling.
- Add **one** replacement regression test for the new composition workflow:
  ```r
  resolved <- resolve_mut_positions(make_test_gene_muts(),
                                    test_path("testdata/test_genome.gff"),
                                    cds_start = 1L)
  aa_muts  <- dplyr::mutate(resolved, pos = aa_start)
  result   <- define_lineage(make_test_biostring(), make_test_ref(),
                             aa_muts, mut_type = "aa")
  # expect equal to the nuc-pathway result (make_test_muts, mut_type = "nuc")
  ```
  (For the fixture gff, `gene_A` starts at nt 1 so `cds_start = 1L`; `gene_A AA2`
  → `aa_start = 2`, `gene_B AA1` → `aa_start = 4`, matching `make_test_aa_muts()`.)
- Add a test asserting `define_lineage(...)` **errors when `mut_type` is omitted**.
- Keep the aa-pathway tests (`190-277`) and ambiguous tests (`281-348`) as-is
  apart from the `mut_type = "nuc"` additions noted above.

`make_test_gene_muts()` and `testdata/test_genome.gff` stay — they are still used
by `tests/testthat/test-resolve_mut_positions.R`.

### 4. Sandbox callers

- `sandbox/ds_0.25_check.R:19` — add `mut_type = "aa"` (ORF table is amino-acid).
- `sandbox/testing ds_1.R:18-24` — the `gff = ...` call no longer works; convert
  to the resolve→`mutate(pos = aa_start)`→`define_lineage(mut_type = "aa")`
  composition (the `resolve_mut_positions(...)` call already on line 16 shows the
  inputs).
- `sandbox/interactive_tests.R:15` — add `mut_type = "nuc"` (`m = make_test_muts()`).
  The now-unused `mg` / `gff` lines (12-13) can be left or removed.

### 5. Regenerate docs

Run `devtools::document()` to update `man/define_lineage.Rd` (drops the `gff`
param). NAMESPACE is unaffected — `define_lineage` stays `@export`, no new imports.

## Verification

1. `devtools::load_all(".")`
2. `testthat::test_file("tests/testthat/test-define_lineage.R")` — all pass,
   including the new composition regression and the missing-`mut_type` error test.
3. `testthat::test_file("tests/testthat/test-resolve_mut_positions.R")` — still
   passes (untouched, fixtures intact).
4. Run `sandbox/ds_0.25_check.R` end-to-end — lineage table should match the
   verified `mut_type = "aa"` distribution from this session
   (NY99 / NY10 / WN02 / SW03 / NS2A_only / NS4B_only across 1585 sequences).
5. Optionally `devtools::check()` for the package as a whole.

## Out of scope / dropped

- Mixed nuc+aa mutation tables in a single call (user does not use these).
- A content-based residue/type sniff check — unreliable, since the project's AA
  residues are all valid IUPAC codes (see Context). Requiring `mut_type` is the
  fix instead.
