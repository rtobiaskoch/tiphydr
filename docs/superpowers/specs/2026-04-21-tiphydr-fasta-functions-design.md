# Design: tiphydr FASTA Functions

**Date:** 2026-04-21  
**Scope:** Six new R functions for the tiphydr package covering FASTA I/O, metadata nesting, alignment trimming, and lineage assignment.

---

## 1. `fasta_read`

**File:** `R/fasta_read.R`

**Signature:** `fasta_read(paths, verbose = TRUE)`

**Arguments:**
- `paths`: character vector of one or more file paths to FASTA files
- `verbose`: logical — if `TRUE`, prints a message per file with sequence count

**Behavior:**
- If `paths` is length 1: applies `Biostrings::readDNAStringSet()` and returns a `DNAStringSet` directly
- If `paths` is length > 1: applies `Biostrings::readDNAStringSet()` to each path via `purrr::map()`, returns a named list of `DNAStringSet` objects
- List names (multi-path case) derived from `tools::file_path_sans_ext(basename(paths))`
- Messages via `message()`, suppressible with `suppressMessages()`

**Returns:**
- Single path → `DNAStringSet`
- Multiple paths → named `list` of `DNAStringSet` objects

---

## 2. `fasta_combine`

**File:** `R/fasta_combine.R`

**Signature:** `fasta_combine(fasta_list, verbose = TRUE)`

**Arguments:**
- `fasta_list`: list of `DNAStringSet` objects (e.g. output of `fasta_read()` with multiple paths)
- `verbose`: logical — if `TRUE`, prints combined sequence count

**Behavior:**
- Concatenates all `DNAStringSet` objects via `do.call(c, fasta_list)`
- Warns via `warning()` if any sequence names are duplicated across inputs
- Stops with informative error if `fasta_list` is empty

**Returns:** single `DNAStringSet`

---

## 3. `fasta_nest`

**File:** `R/fasta_nest.R`

**Signature:** `fasta_nest(df, biostring, id_col, delim, verbose = TRUE)`

**Arguments:**
- `df`: metadata dataframe
- `biostring`: `DNAStringSet` of sequences
- `id_col`: string — name of column in `df` whose values are matched against sequence names
- `delim`: string — delimiter character used to split sequence names before matching (e.g. `"|"` splits `"WNV|2021|CO|CSU123"` into parts; the `id_col` value is matched against any part using exact string equality)
- `verbose`: logical — if `TRUE`, prints match summary

**Behavior:**
- For each row, splits every sequence name in `biostring` on `delim` using `stringr::str_split()`, then checks whether the row's `id_col` value exactly matches any resulting part
- Stores the matched length-1 `DNAStringSet` in a new `sequence` list-column (stored as `NULL` if no match)
- Warns via `warning()` on rows with no match
- Warns via `warning()` on rows matching more than one sequence name (first match used)
- Messages via `message()` with count of matched, unmatched rows

**Returns:** tibble with all original columns plus a `sequence` list-column of length-1 `DNAStringSet` objects

**Notes:**
- Analogous to `sf` geometry column — one sequence per row
- Fully compatible with `dplyr::filter()`, `mutate()`, `group_by()`
- `delim` avoids regex metacharacter footguns from raw sequence name matching; matching within split parts is exact (no regex)

---

## 4. `fasta_unnest`

**File:** `R/fasta_unnest.R`

**Signature:** `fasta_unnest(nested_df, verbose = TRUE)`

**Arguments:**
- `nested_df`: tibble with a `sequence` list-column as produced by `fasta_nest()`
- `verbose`: logical — if `TRUE`, prints count of rows dropped due to `NULL` sequences

**Behavior:**
- Drops rows where `sequence` is `NULL`, warns if any are dropped
- Extracts and concatenates the `sequence` list-column into a single `DNAStringSet` via `do.call(c, ...)`
- Returns the metadata (all columns except `sequence`) and the `DNAStringSet` as separate named elements

**Returns:** `list(metadata = tibble, biostring = DNAStringSet)`

---

## 5. `fasta_trim_ref`

**File:** `R/fasta_trim_ref.R`

**Signature:** `fasta_trim_ref(fasta, ref, alignment = NULL, drop = FALSE, verbose = TRUE)`

**Arguments:**
- `fasta`: `DNAStringSet` — sequences to trim
- `ref`: single-sequence `DNAStringSet` — reference used to determine trim coordinates
- `alignment`: optional `DNAStringSet` — pre-computed alignment to add sequences into. If provided, uses MAFFT `--add` profile alignment mode to insert new sequences without disturbing existing columns
- `drop`: logical (default `FALSE`). If `TRUE` and `alignment` is provided, sequences present in `alignment` but absent from `fasta` are removed from output
- `verbose`: logical — if `TRUE`, prints progress messages

**MAFFT dependency check:**
- On load, checks `Sys.which("mafft") != ""` before any alignment call
- If MAFFT is absent, stops with: `"MAFFT not found on PATH. Install from https://mafft.cbrc.jp/alignment/software/ or use the 'msa' Bioconductor package for basic alignment (note: profile alignment via --add requires MAFFT)."`

**Behavior:**
1. Validate `ref` has length exactly 1; stop otherwise
2. Check MAFFT availability
3. If `alignment = NULL`: combine `c(ref, fasta)`, convert to `DNAbin` via `ape::as.DNAbin()`, align with `ips::mafft()`, convert back to character matrix
4. If `alignment` provided: identify sequences in `fasta` not already present in `alignment` by name; use `ips::mafft(..., add = alignment)` to add new sequences into the existing alignment using profile mode
5. Find columns in the alignment where the ref row contains `A/T/C/G` (non-gap positions) — these are the reference coordinates
6. Subset alignment matrix to those columns; after this step, alignment position N equals reference position N (ungapped) — this is the invariant that makes `define_lineage` coordinates unambiguous
7. Remove ref row
8. If `drop = TRUE`: remove sequences not in `fasta` from output
9. Collapse matrix rows to strings, convert to `DNAStringSet`

**Validation:**
- Warns via `warning()` if trimmed sequence length does not equal `nchar(ref)` (ungapped)
- Stops if `ref` has length > 1

**Messages (via `message()`, only if `verbose = TRUE`):**
- Original average sequence length
- Alignment length post-MAFFT
- Trim start and end positions
- Final trimmed length
- If `alignment` provided: count of sequences added; if `drop = TRUE`, count dropped

**Returns:** `DNAStringSet` of trimmed sequences (ref excluded)

---

## 6. `define_lineage`

**File:** `R/define_lineage.R`

**Signature:** `define_lineage(alignment, ref, muts, type = "nuc", verbose = TRUE)`

**Arguments:**
- `alignment`: `DNAStringSet` — pre-aligned sequences, all equal length. Expected to be output of `fasta_trim_ref()` so that alignment positions equal ungapped reference positions
- `ref`: single-sequence `DNAStringSet` — establishes coordinate system; used to locate the first start codon (ATG) for AA translation
- `muts`: long-format dataframe with columns:
  - `lineage`: character — lineage label
  - `pos`: integer — position in **reference coordinates** (ungapped); equivalent to alignment position after `fasta_trim_ref()`
  - `residue`: character — expected nucleotide (e.g. `"A"`) when `type = "nuc"`, or expected amino acid (e.g. `"K"`) when `type = "aa"`
- `type`: `"nuc"` (default) or `"aa"`
- `verbose`: logical — if `TRUE`, prints lineage count table and unknowns

**Coordinate system contract:**
Positions in `muts$pos` are ungapped reference coordinates. Because `fasta_trim_ref()` removes all gap-only columns (positions where ref has a gap), alignment column N corresponds exactly to reference nucleotide N. This function assumes that contract has been satisfied; it does not re-validate it, but documents the assumption clearly in the roxygen header.

**Behavior:**
1. If `type = "aa"`:
   - Find the first ATG in the ref sequence using `Biostrings::matchPattern("ATG", ref[[1]])` and take the start position of the first match as `cds_start`
   - Translate each sequence in `alignment` from `cds_start` using `Biostrings::translate(..., if.fuzzy.codon = "X")`
   - `muts$pos` are now amino acid positions relative to the translated CDS
2. For each unique lineage in `muts`, collect its required `(pos, residue)` pairs
3. For each sequence × lineage: check whether **all** required `(pos, residue)` pairs are present (strict match)
4. Assignment rules:
   - Matches exactly one lineage → assign that lineage
   - Matches zero lineages → assign `"unknown"`
   - Matches multiple lineages → assign the lineage with the most required mutations (most specific); if tied, assign first alphabetically and warn
5. Return dataframe with `strain` and `lineage` columns

**Messages (via `message()`, only if `verbose = TRUE`):**
- Table of lineage counts
- Count of `"unknown"` sequences; print their names if any

**Returns:** dataframe with columns `strain` (character) and `lineage` (character)

---

## Package dependencies introduced

| Package | Use |
|---------|-----|
| `Biostrings` | All sequence operations, `matchPattern()` for ATG detection |
| `purrr` | `map()` in `fasta_read` |
| `stringr` | `str_split()` for delimiter-based name matching in `fasta_nest` |
| `dplyr` | Dataframe manipulation throughout |
| `ape` | `as.DNAbin()` conversion for MAFFT input |
| `ips` | `mafft()` alignment, including `--add` profile mode |
| `tools` | `file_path_sans_ext()` in `fasta_read` |

---

## File layout

```
R/
  fasta_read.R
  fasta_combine.R
  fasta_nest.R
  fasta_unnest.R
  fasta_trim_ref.R
  define_lineage.R
```

Each file contains exactly one exported function with roxygen2 documentation.

---

## RSE notes

- All messages use `message()` (stderr), suppressible with `suppressMessages()`
- All functions are pure — no side effects, no global state modification
- `verbose = TRUE` is a standard parameter across all functions for pipeline use
- `fasta_trim_ref()` output satisfies the coordinate contract required by `define_lineage()`
- MAFFT system dependency is checked at call time with a helpful install message
