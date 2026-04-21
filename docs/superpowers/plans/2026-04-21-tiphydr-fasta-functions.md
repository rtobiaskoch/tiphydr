# tiphydr FASTA Functions — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement six R package functions (`fasta_read`, `fasta_combine`, `fasta_nest`, `fasta_unnest`, `fasta_trim_ref`, `define_lineage`) plus a synthetic WNV test dataset, following TDD throughout.

**Architecture:** Each function lives in its own `R/` file and is pure (no side effects beyond optional `message()` output). Functions compose left-to-right: `fasta_read → fasta_combine → fasta_nest/unnest → fasta_trim_ref → define_lineage`. The critical coordinate invariant is that `fasta_trim_ref()` output guarantees alignment position N = ungapped reference position N, which `define_lineage()` relies on for mutation lookup without additional coordinate mapping. Test fixtures are defined once in `tests/testthat/helper-fixtures.R` and reused across all test files.

**Tech Stack:** R, Biostrings, purrr, stringr, dplyr, ape, ips (MAFFT wrapper), system2 (MAFFT --add profile mode), testthat

---

## File Map

**Create:**
- `data-raw/create_wnv_test_data.R` — script to generate test FASTA and CSV files
- `inst/extdata/wnv_ref.fasta` — 84nt synthetic WNV-like reference (ATG at pos 1)
- `inst/extdata/wnv_seqs.fasta` — 5 test sequences with known mutations
- `inst/extdata/wnv_metadata.csv` — strain metadata (strain_id, date, location, host)
- `inst/extdata/wnv_muts.csv` — mutation definitions for two lineages (NY10, PARTIAL_A)
- `tests/testthat.R` — testthat entry point
- `tests/testthat/helper-fixtures.R` — shared in-memory test fixtures
- `R/fasta_read.R`
- `R/fasta_combine.R`
- `R/fasta_nest.R`
- `R/fasta_unnest.R`
- `R/fasta_trim_ref.R`
- `R/define_lineage.R`
- `tests/testthat/test-fasta_read.R`
- `tests/testthat/test-fasta_combine.R`
- `tests/testthat/test-fasta_nest.R`
- `tests/testthat/test-fasta_unnest.R`
- `tests/testthat/test-fasta_trim_ref.R`
- `tests/testthat/test-define_lineage.R`

**Modify:**
- `DESCRIPTION` — add Imports and Suggests

---

## Task 1: Package infrastructure

**Files:**
- Modify: `DESCRIPTION`
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-fixtures.R`

- [ ] **Step 1: Add dependencies to DESCRIPTION**

Open `DESCRIPTION` and replace the existing `Imports:` block (or add one after `Encoding:`):

```
Imports:
    Biostrings,
    purrr,
    stringr,
    dplyr,
    ape,
    ips,
    tools,
    tibble
Suggests:
    testthat (>= 3.0.0)
```

- [ ] **Step 2: Create tests/testthat.R**

```r
library(testthat)
library(tiphydr)
test_check("tiphydr")
```

- [ ] **Step 3: Create shared test fixtures**

Create `tests/testthat/helper-fixtures.R`. These build in-memory objects so all tests run without file I/O (faster, no path dependencies):

```r
# 84nt synthetic WNV-like reference starting with ATG
# Position 10 = G (wildtype), position 25 = T (wildtype)
# Verified: nchar(.ref_seq) == 84
.ref_seq <- paste0(
  "ATGAAACCCG",  # pos 1-10  — pos 10 = G
  "GGTTTTAAAT",  # pos 11-20
  "GGCTTACGAT",  # pos 21-30 — pos 25 = T
  "TCCAAAGCTT",  # pos 31-40
  "GCGAAATTTG",  # pos 41-50
  "GCCCAAATTC",  # pos 51-60
  "CTTGCATGCA",  # pos 61-70
  "AAGCTTGGAA"   # pos 71-80
)  # 80nt — add 4 to round to 84
.ref_seq <- paste0(.ref_seq, "ATCG")  # 81-84

# Build a test sequence by substituting at specific positions
.mutate_seq <- function(base = .ref_seq, pos10 = "G", pos25 = "T") {
  s <- base
  substr(s, 10, 10) <- pos10
  substr(s, 25, 25) <- pos25
  s
}

#' In-memory DNAStringSet with 5 test sequences
make_test_biostring <- function() {
  Biostrings::DNAStringSet(c(
    `WNV|2021|CO|NY10_001`    = .mutate_seq(pos10 = "A", pos25 = "C"),  # NY10: both mutations
    `WNV|2020|CO|NS2A_002`    = .mutate_seq(pos10 = "A", pos25 = "T"),  # pos10 only
    `WNV|2022|WY|NS4B_003`    = .mutate_seq(pos10 = "G", pos25 = "C"),  # pos25 only
    `WNV|2019|CO|OTHER_004`   = .mutate_seq(pos10 = "G", pos25 = "T"),  # no mutations
    `WNV|2021|NM|EXTRA_005`   = .mutate_seq(pos10 = "A", pos25 = "C")   # not in metadata
  ))
}

#' In-memory reference DNAStringSet (length 1)
make_test_ref <- function() {
  r <- Biostrings::DNAStringSet(.ref_seq)
  names(r) <- "WNV_REF_KU877344"
  r
}

#' Metadata dataframe — 4 rows, no EXTRA_005
make_test_metadata <- function() {
  tibble::tibble(
    strain_id = c("NY10_001", "NS2A_002", "NS4B_003", "OTHER_004"),
    date      = as.Date(c("2021-07-15", "2020-08-01", "2022-06-20", "2019-05-10")),
    location  = c("Colorado", "Colorado", "Wyoming", "Colorado"),
    host      = c("Culex tarsalis", "Culex pipiens", "Culex tarsalis", "Bird")
  )
}

#' Mutation definitions: NY10 (2 mutations, strict), PARTIAL_A (1 mutation)
make_test_muts <- function() {
  tibble::tibble(
    lineage = c("NY10",  "NY10",  "PARTIAL_A"),
    pos     = c(10L,     25L,     10L),
    residue = c("A",     "C",     "A")
  )
}
```

- [ ] **Step 4: Commit**

```bash
git add DESCRIPTION tests/
git commit -m "feat: add test infrastructure for tiphydr FASTA functions"
```

---

## Task 2: WNV test dataset

**Files:**
- Create: `data-raw/create_wnv_test_data.R`
- Create: `inst/extdata/wnv_ref.fasta`, `wnv_seqs.fasta`, `wnv_metadata.csv`, `wnv_muts.csv`

- [ ] **Step 1: Create data-raw generation script**

Create `data-raw/create_wnv_test_data.R`:

```r
# Generates synthetic WNV test data for tiphydr examples and tests.
# Run once: source("data-raw/create_wnv_test_data.R")
# Outputs go to inst/extdata/ and are committed to the package.
#
# Sequence design:
#   - 84nt WNV-like reference, ATG at position 1
#   - Mutation site 1: position 10 (G=wildtype, A=NY10/PARTIAL_A)
#   - Mutation site 2: position 25 (T=wildtype, C=NY10)
#   - Lineage NY10 requires BOTH pos10=A AND pos25=C (strict match)
#   - Lineage PARTIAL_A requires only pos10=A

library(Biostrings)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

# ── Reference ─────────────────────────────────────────────────────────────────
ref_seq <- paste0(
  "ATGAAACCCG",  # pos 1-10  — pos 10 = G (wildtype)
  "GGTTTTAAAT",  # pos 11-20
  "GGCTTACGAT",  # pos 21-30 — pos 25 = T (wildtype)
  "TCCAAAGCTT",  # pos 31-40
  "GCGAAATTTG",  # pos 41-50
  "GCCCAAATTC",  # pos 51-60
  "CTTGCATGCA",  # pos 61-70
  "AAGCTTGGAA",  # pos 71-80
  "ATCG"         # pos 81-84
)
stopifnot(nchar(ref_seq) == 84)

ref <- DNAStringSet(c(WNV_REF_KU877344 = ref_seq))
writeXStringSet(ref, "inst/extdata/wnv_ref.fasta")
message("Wrote wnv_ref.fasta (", nchar(ref_seq), " nt)")

# ── Test sequences ─────────────────────────────────────────────────────────────
make_seq <- function(pos10, pos25) {
  s <- ref_seq
  substr(s, 10, 10) <- pos10
  substr(s, 25, 25) <- pos25
  s
}

seqs <- DNAStringSet(c(
  `WNV|2021|CO|NY10_001`     = make_seq("A", "C"),  # NY10: both mutations
  `WNV|2020|CO|NS2A_002`     = make_seq("A", "T"),  # pos10=A only → PARTIAL_A, not NY10
  `WNV|2022|WY|NS4B_003`     = make_seq("G", "C"),  # pos25=C only → unknown
  `WNV|2019|CO|OTHER_004`    = make_seq("G", "T"),  # no mutations → unknown
  `WNV|2021|CO|PARTIAL_005`  = make_seq("A", "T")   # pos10=A only → PARTIAL_A
))
writeXStringSet(seqs, "inst/extdata/wnv_seqs.fasta")
message("Wrote wnv_seqs.fasta (", length(seqs), " sequences)")

# ── Metadata ──────────────────────────────────────────────────────────────────
metadata <- data.frame(
  strain_id = c("NY10_001",       "NS2A_002",       "NS4B_003",       "OTHER_004",    "PARTIAL_005"),
  date      = c("2021-07-15",     "2020-08-01",     "2022-06-20",     "2019-05-10",   "2021-09-15"),
  location  = c("Colorado",       "Colorado",       "Wyoming",        "Colorado",     "Colorado"),
  host      = c("Culex tarsalis", "Culex pipiens",  "Culex tarsalis", "Bird",         "Culex tarsalis"),
  stringsAsFactors = FALSE
)
write.csv(metadata, "inst/extdata/wnv_metadata.csv", row.names = FALSE)
message("Wrote wnv_metadata.csv (", nrow(metadata), " rows)")

# ── Mutation definitions ───────────────────────────────────────────────────────
# Positions are reference coordinates (ungapped).
# After fasta_trim_ref(), alignment position N = reference position N.
muts <- data.frame(
  lineage = c("NY10",  "NY10",  "PARTIAL_A"),
  pos     = c(10L,     25L,     10L),
  residue = c("A",     "C",     "A"),
  stringsAsFactors = FALSE
)
write.csv(muts, "inst/extdata/wnv_muts.csv", row.names = FALSE)
message("Wrote wnv_muts.csv (", nrow(muts), " rows)")

message("\nDone. Expected define_lineage() assignments with type='nuc':")
message("  NY10_001   -> NY10 (both mutations; NY10 is more specific than PARTIAL_A)")
message("  NS2A_002   -> PARTIAL_A (pos10=A only)")
message("  NS4B_003   -> unknown (pos25=C but not pos10=A)")
message("  OTHER_004  -> unknown (no mutations)")
message("  PARTIAL_005 -> PARTIAL_A (pos10=A only)")
```

- [ ] **Step 2: Run the script**

```r
source("data-raw/create_wnv_test_data.R")
```

Expected output:
```
Wrote wnv_ref.fasta (84 nt)
Wrote wnv_seqs.fasta (5 sequences)
Wrote wnv_metadata.csv (5 rows)
Wrote wnv_muts.csv (3 rows)

Done. Expected define_lineage() assignments with type='nuc':
  NY10_001   -> NY10 (both mutations; NY10 is more specific than PARTIAL_A)
  NS2A_002   -> PARTIAL_A (pos10=A only)
  NS4B_003   -> unknown (pos25=C but not pos10=A)
  OTHER_004  -> unknown (no mutations)
  PARTIAL_005 -> PARTIAL_A (pos10=A only)
```

Verify:
```bash
ls inst/extdata/
# wnv_metadata.csv  wnv_muts.csv  wnv_ref.fasta  wnv_seqs.fasta
```

- [ ] **Step 3: Commit**

```bash
git add inst/extdata/ data-raw/
git commit -m "feat: add synthetic WNV test dataset (84nt, NY10/PARTIAL_A lineages)"
```

---

## Task 3: fasta_read

**Files:**
- Create: `R/fasta_read.R`
- Create: `tests/testthat/test-fasta_read.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-fasta_read.R`:

```r
test_that("single path returns DNAStringSet directly (not a list)", {
  path   <- system.file("extdata", "wnv_ref.fasta", package = "tiphydr")
  result <- suppressMessages(fasta_read(path))
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), 1L)
})

test_that("multiple paths return named list of DNAStringSets", {
  ref_path <- system.file("extdata", "wnv_ref.fasta",  package = "tiphydr")
  seq_path <- system.file("extdata", "wnv_seqs.fasta", package = "tiphydr")
  result   <- suppressMessages(fasta_read(c(ref_path, seq_path)))
  expect_type(result, "list")
  expect_length(result, 2L)
  expect_named(result, c("wnv_ref", "wnv_seqs"))
  expect_s4_class(result[["wnv_ref"]],  "DNAStringSet")
  expect_s4_class(result[["wnv_seqs"]], "DNAStringSet")
  expect_equal(length(result[["wnv_seqs"]]), 5L)
})

test_that("non-existent file throws informative error", {
  expect_error(fasta_read("ghost.fasta"), "Files not found")
})

test_that("verbose=FALSE suppresses all messages", {
  path <- system.file("extdata", "wnv_ref.fasta", package = "tiphydr")
  expect_silent(fasta_read(path, verbose = FALSE))
})
```

- [ ] **Step 2: Run tests — expect all 4 to FAIL**

```bash
Rscript -e "devtools::test(filter = 'fasta_read')"
```
Expected: `Error in fasta_read(...) : could not find function "fasta_read"`

- [ ] **Step 3: Implement fasta_read**

Create `R/fasta_read.R`:

```r
#' Read one or more FASTA files into DNAStringSet objects
#'
#' @param paths character vector of file paths to FASTA files
#' @param verbose logical; if TRUE prints sequence count per file (default TRUE)
#'
#' @return Single path: a DNAStringSet. Multiple paths: a named list of
#'   DNAStringSets, list names derived from filenames without extension.
#' @export
fasta_read <- function(paths, verbose = TRUE) {

  if (!is.character(paths) || length(paths) == 0) {
    stop("paths must be a non-empty character vector")
  }

  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }

  read_one <- function(path) {
    seqs <- Biostrings::readDNAStringSet(path)
    if (verbose) message("Read ", length(seqs), " sequences from ", basename(path))
    seqs
  }

  if (length(paths) == 1L) return(read_one(paths))

  result       <- purrr::map(paths, read_one)
  names(result) <- tools::file_path_sans_ext(basename(paths))
  result
}
```

- [ ] **Step 4: Run tests — expect all 4 to PASS**

```bash
Rscript -e "devtools::test(filter = 'fasta_read')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ]`

- [ ] **Step 5: Commit**

```bash
git add R/fasta_read.R tests/testthat/test-fasta_read.R
git commit -m "feat: add fasta_read()"
```

---

## Task 4: fasta_combine

**Files:**
- Create: `R/fasta_combine.R`
- Create: `tests/testthat/test-fasta_combine.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-fasta_combine.R`:

```r
test_that("combines a list of DNAStringSets into one", {
  ref  <- make_test_ref()
  seqs <- make_test_biostring()
  result <- suppressMessages(fasta_combine(list(ref, seqs)))
  expect_s4_class(result, "DNAStringSet")
  expect_equal(length(result), length(ref) + length(seqs))
})

test_that("warns on duplicate sequence names", {
  seqs <- make_test_biostring()
  expect_warning(
    suppressMessages(fasta_combine(list(seqs, seqs))),
    "Duplicate sequence names"
  )
})

test_that("stops on empty list", {
  expect_error(fasta_combine(list()), "non-empty list")
})

test_that("verbose=FALSE suppresses messages", {
  expect_silent(fasta_combine(list(make_test_biostring()), verbose = FALSE))
})
```

- [ ] **Step 2: Run tests — expect all 4 to FAIL**

```bash
Rscript -e "devtools::test(filter = 'fasta_combine')"
```

- [ ] **Step 3: Implement fasta_combine**

Create `R/fasta_combine.R`:

```r
#' Combine a list of DNAStringSet objects into one
#'
#' @param fasta_list list of DNAStringSet objects (e.g. output of fasta_read()
#'   with multiple paths)
#' @param verbose logical; if TRUE prints combined sequence count (default TRUE)
#'
#' @return DNAStringSet
#' @export
fasta_combine <- function(fasta_list, verbose = TRUE) {

  if (!is.list(fasta_list) || length(fasta_list) == 0) {
    stop("fasta_list must be a non-empty list of DNAStringSet objects")
  }

  combined <- do.call(c, fasta_list)

  dupe_names <- names(combined)[duplicated(names(combined))]
  if (length(dupe_names) > 0) {
    warning("Duplicate sequence names found: ", paste(unique(dupe_names), collapse = ", "))
  }

  if (verbose) {
    message("Combined ", length(fasta_list), " DNAStringSets into ",
            length(combined), " total sequences")
  }

  combined
}
```

- [ ] **Step 4: Run tests — expect all 4 to PASS**

```bash
Rscript -e "devtools::test(filter = 'fasta_combine')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 4 ]`

- [ ] **Step 5: Commit**

```bash
git add R/fasta_combine.R tests/testthat/test-fasta_combine.R
git commit -m "feat: add fasta_combine()"
```

---

## Task 5: fasta_nest

**Files:**
- Create: `R/fasta_nest.R`
- Create: `tests/testthat/test-fasta_nest.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-fasta_nest.R`:

```r
test_that("adds sequence list-column to metadata", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|")
  ))
  expect_true("sequence" %in% names(result))
  expect_equal(nrow(result), nrow(meta))
  expect_s4_class(result$sequence[[1]], "DNAStringSet")
  expect_equal(length(result$sequence[[1]]), 1L)
})

test_that("matched sequence name matches id_col value via delimiter split", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|")
  ))
  # NY10_001 should match WNV|2021|CO|NY10_001
  expect_equal(names(result$sequence[[1]]), "WNV|2021|CO|NY10_001")
})

test_that("unmatched row stores NULL and raises warning", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  # EXTRA_005 is in seqs but not meta; OTHER_004 is in meta — it will match
  # Remove a row whose id is NOT in seqs to force a NULL
  meta_extra <- tibble::add_row(meta, strain_id = "GHOST_999",
                                 date = as.Date("2021-01-01"),
                                 location = "Colorado", host = "Bird")
  expect_warning(
    suppressMessages(fasta_nest(meta_extra, seqs, id_col = "strain_id", delim = "|")),
    regexp = "no matching sequence"
  )
  result <- suppressMessages(suppressWarnings(
    fasta_nest(meta_extra, seqs, id_col = "strain_id", delim = "|")
  ))
  expect_null(result$sequence[[nrow(result)]])
})

test_that("wrong id_col name throws informative error", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  expect_error(
    fasta_nest(meta, seqs, id_col = "bad_col", delim = "|"),
    "not found in df columns"
  )
})

test_that("verbose=FALSE suppresses messages", {
  seqs <- make_test_biostring()
  meta <- make_test_metadata()
  expect_no_message(suppressWarnings(
    fasta_nest(meta, seqs, id_col = "strain_id", delim = "|", verbose = FALSE)
  ))
})
```

- [ ] **Step 2: Run tests — expect all 5 to FAIL**

```bash
Rscript -e "devtools::test(filter = 'fasta_nest')"
```

- [ ] **Step 3: Implement fasta_nest**

Create `R/fasta_nest.R`:

```r
#' Nest a DNAStringSet as a list-column in a metadata dataframe
#'
#' Analogous to an sf geometry column: one sequence per metadata row.
#' The resulting `sequence` list-column is compatible with dplyr::filter(),
#' mutate(), and group_by(). Use fasta_unnest() to reverse.
#'
#' Matching is done by splitting each sequence name on `delim` and checking
#' whether the row's `id_col` value exactly equals any resulting part.
#' This avoids regex metacharacter issues common in FASTA sequence names.
#'
#' @param df dataframe with metadata
#' @param biostring DNAStringSet of sequences
#' @param id_col string; name of column in df to match against sequence name parts
#' @param delim string; delimiter used to split sequence names before matching
#'   (e.g. "|" splits "WNV|2021|CO|CSU123" into c("WNV","2021","CO","CSU123"))
#' @param verbose logical; if TRUE prints match summary (default TRUE)
#'
#' @return tibble with all original columns plus a `sequence` list-column of
#'   length-1 DNAStringSet objects (NULL where no match found)
#' @export
fasta_nest <- function(df, biostring, id_col, delim, verbose = TRUE) {

  if (!id_col %in% names(df)) {
    stop("id_col '", id_col, "' not found in df columns: ",
         paste(names(df), collapse = ", "))
  }

  seq_names   <- names(biostring)
  # Pre-split all sequence names once — reused for every row
  split_names <- stringr::str_split(seq_names, stringr::fixed(delim))

  match_sequence <- function(id_value) {
    matched_idx <- which(
      purrr::map_lgl(split_names, ~ id_value %in% .x)
    )

    if (length(matched_idx) == 0L) return(NULL)

    if (length(matched_idx) > 1L) {
      warning("id '", id_value, "' matched ", length(matched_idx),
              " sequence names; using first: ", seq_names[matched_idx[1]])
    }

    biostring[matched_idx[1]]
  }

  id_values <- df[[id_col]]
  sequences <- purrr::map(id_values, match_sequence)

  n_matched   <- sum(!purrr::map_lgl(sequences, is.null))
  n_unmatched <- nrow(df) - n_matched

  if (n_unmatched > 0L) {
    unmatched_ids <- id_values[purrr::map_lgl(sequences, is.null)]
    warning(n_unmatched, " rows had no matching sequence: ",
            paste(unmatched_ids, collapse = ", "))
  }

  if (verbose) {
    message("Matched ", n_matched, " of ", nrow(df), " rows to sequences")
  }

  dplyr::mutate(df, sequence = sequences)
}
```

- [ ] **Step 4: Run tests — expect all 5 to PASS**

```bash
Rscript -e "devtools::test(filter = 'fasta_nest')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ]`

- [ ] **Step 5: Commit**

```bash
git add R/fasta_nest.R tests/testthat/test-fasta_nest.R
git commit -m "feat: add fasta_nest()"
```

---

## Task 6: fasta_unnest

**Files:**
- Create: `R/fasta_unnest.R`
- Create: `tests/testthat/test-fasta_unnest.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-fasta_unnest.R`:

```r
test_that("returns list with metadata and biostring elements", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  result <- suppressMessages(fasta_unnest(nested))
  expect_type(result, "list")
  expect_named(result, c("metadata", "biostring"))
  expect_s3_class(result$metadata, "data.frame")
  expect_s4_class(result$biostring, "DNAStringSet")
})

test_that("metadata does not contain sequence column", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  result <- suppressMessages(fasta_unnest(nested))
  expect_false("sequence" %in% names(result$metadata))
})

test_that("biostring length equals number of non-NULL rows", {
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(make_test_metadata(), make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  # All 4 metadata rows match — biostring should have 4 sequences
  result <- suppressMessages(fasta_unnest(nested))
  expect_equal(length(result$biostring), 4L)
})

test_that("NULL rows are dropped and warned", {
  meta_extra <- tibble::add_row(make_test_metadata(),
                                 strain_id = "GHOST_999",
                                 date = as.Date("2021-01-01"),
                                 location = "Colorado", host = "Bird")
  nested <- suppressMessages(suppressWarnings(
    fasta_nest(meta_extra, make_test_biostring(),
               id_col = "strain_id", delim = "|")
  ))
  expect_message(fasta_unnest(nested), "Dropping 1 rows")
  result <- suppressMessages(fasta_unnest(nested))
  expect_equal(nrow(result$metadata), 4L)
})

test_that("stops if sequence column is missing", {
  meta <- make_test_metadata()
  expect_error(fasta_unnest(meta), "must have a 'sequence' list-column")
})
```

- [ ] **Step 2: Run tests — expect all 5 to FAIL**

```bash
Rscript -e "devtools::test(filter = 'fasta_unnest')"
```

- [ ] **Step 3: Implement fasta_unnest**

Create `R/fasta_unnest.R`:

```r
#' Extract metadata and sequences from a nested sequence dataframe
#'
#' Inverse of fasta_nest(). Rows with NULL sequences (no FASTA match) are
#' dropped with a warning. Returns metadata and sequences as separate elements.
#'
#' @param nested_df tibble with a `sequence` list-column as produced by
#'   fasta_nest()
#' @param verbose logical; if TRUE prints count of rows dropped (default TRUE)
#'
#' @return named list: metadata (tibble, sequence column removed) and
#'   biostring (DNAStringSet)
#' @export
fasta_unnest <- function(nested_df, verbose = TRUE) {

  if (!"sequence" %in% names(nested_df)) {
    stop("nested_df must have a 'sequence' list-column (produced by fasta_nest())")
  }

  null_mask <- purrr::map_lgl(nested_df$sequence, is.null)
  n_null    <- sum(null_mask)

  if (n_null > 0L) {
    if (verbose) message("Dropping ", n_null, " rows with NULL sequences")
    nested_df <- nested_df[!null_mask, ]
  }

  biostring <- do.call(c, nested_df$sequence)
  metadata  <- dplyr::select(nested_df, -sequence)

  list(metadata = metadata, biostring = biostring)
}
```

- [ ] **Step 4: Run tests — expect all 5 to PASS**

```bash
Rscript -e "devtools::test(filter = 'fasta_unnest')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ]`

- [ ] **Step 5: Commit**

```bash
git add R/fasta_unnest.R tests/testthat/test-fasta_unnest.R
git commit -m "feat: add fasta_unnest()"
```

---

## Task 7: fasta_trim_ref

**Files:**
- Create: `R/fasta_trim_ref.R`
- Create: `tests/testthat/test-fasta_trim_ref.R`

**Note:** All tests in this task require MAFFT on PATH. They use `testthat::skip_if()` to skip gracefully when MAFFT is absent. The profile alignment path (`--add`) uses `system2("mafft", ...)` with temp files because MAFFT's `--add` flag requires controlling the command directly.

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-fasta_trim_ref.R`:

```r
skip_if_no_mafft <- function() {
  testthat::skip_if(nchar(Sys.which("mafft")) == 0, "MAFFT not on PATH")
}

test_that("returns DNAStringSet with ref excluded", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()
  result <- suppressMessages(fasta_trim_ref(seqs, ref))
  expect_s4_class(result, "DNAStringSet")
  expect_false(names(ref) %in% names(result))
})

test_that("all trimmed sequences equal reference length", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()
  result <- suppressMessages(fasta_trim_ref(seqs, ref))
  # Reference is 84nt; all trimmed seqs must also be 84nt (gaps allowed)
  expect_true(all(Biostrings::width(result) == Biostrings::width(ref)))
})

test_that("sequences with leading insertions trim to reference coordinates", {
  skip_if_no_mafft()
  ref <- make_test_ref()
  # Add 5 extra nt before the reference region in one sequence
  padded_seq <- Biostrings::DNAStringSet(
    c(test_padded = paste0("AAAAA", as.character(ref[[1]])))
  )
  result <- suppressMessages(fasta_trim_ref(padded_seq, ref))
  # After trimming, the padded sequence should equal the reference content
  # (the 5 extra nt are insertions relative to ref and get removed)
  expect_equal(Biostrings::width(result[[1]]), Biostrings::width(ref))
})

test_that("ref of length > 1 throws error", {
  two_refs <- c(make_test_ref(), make_test_ref())
  expect_error(fasta_trim_ref(make_test_biostring(), two_refs), "length 1")
})

test_that("absent MAFFT gives informative stop message", {
  # Mock Sys.which to return empty string — only run if MAFFT actually present
  # so we can test the error path without MAFFT
  skip_if_no_mafft()
  with_mocked_bindings(
    `Sys.which` = function(x) if (x == "mafft") "" else base::Sys.which(x),
    expect_error(fasta_trim_ref(make_test_biostring(), make_test_ref()),
                 "MAFFT not found"),
    .package = "base"
  )
})

test_that("drop=TRUE removes extra sequences from alignment", {
  skip_if_no_mafft()
  ref    <- make_test_ref()
  seqs   <- make_test_biostring()

  # First align all 5 sequences
  full_aln <- suppressMessages(fasta_trim_ref(seqs, ref))

  # Now provide a subset of 2 sequences as fasta, use full alignment
  # With drop=TRUE, only the 2 sequences should appear in output
  subset_seqs <- seqs[1:2]
  result <- suppressMessages(
    fasta_trim_ref(subset_seqs, ref, alignment = full_aln, drop = TRUE)
  )
  expect_equal(length(result), 2L)
  expect_true(all(names(result) %in% names(subset_seqs)))
})
```

- [ ] **Step 2: Run tests — expect 4 to SKIP (no MAFFT) or FAIL (MAFFT present)**

```bash
Rscript -e "devtools::test(filter = 'fasta_trim_ref')"
```

- [ ] **Step 3: Implement fasta_trim_ref**

Create `R/fasta_trim_ref.R`:

```r
#' Align sequences to a reference and trim to reference coordinates
#'
#' Uses MAFFT for alignment. When fasta is the first alignment, all sequences
#' plus the reference are aligned together with ips::mafft(). When an existing
#' alignment is provided, new sequences are inserted via MAFFT --add profile
#' mode, which preserves existing alignment columns exactly.
#'
#' Coordinate invariant: after trimming, alignment position N equals ungapped
#' reference position N. This is the contract required by define_lineage().
#'
#' @param fasta DNAStringSet; sequences to align and trim
#' @param ref length-1 DNAStringSet; reference sequence defining trim coordinates
#' @param alignment optional DNAStringSet; existing alignment to add sequences
#'   into via MAFFT --add profile mode (NULL = build from scratch)
#' @param drop logical; if TRUE and alignment provided, remove sequences present
#'   in alignment but absent from fasta (default FALSE)
#' @param verbose logical; if TRUE prints progress messages (default TRUE)
#'
#' @return DNAStringSet of trimmed sequences (reference excluded)
#' @export
fasta_trim_ref <- function(fasta, ref, alignment = NULL, drop = FALSE,
                           verbose = TRUE) {

  if (length(ref) != 1L) stop("ref must be a single-sequence DNAStringSet (length 1)")

  if (nchar(Sys.which("mafft")) == 0L) {
    stop(
      "MAFFT not found on PATH.\n",
      "Install from: https://mafft.cbrc.jp/alignment/software/\n",
      "Alternative for basic alignment (no --add): BiocManager::install('msa')\n",
      "Profile alignment via --add requires MAFFT."
    )
  }

  ref_name <- names(ref)

  if (verbose) {
    message("Input: ", length(fasta), " sequences, avg length ",
            round(mean(Biostrings::width(fasta)), 1), " nt")
    message("Reference: ", ref_name, " (", Biostrings::width(ref), " nt)")
  }

  # ── Build alignment ──────────────────────────────────────────────────────────
  if (is.null(alignment)) {
    # Full alignment from scratch
    combined     <- c(ref, fasta)
    combined_bin <- ape::as.DNAbin(combined)
    aligned_bin  <- ips::mafft(combined_bin)
    aligned_mat  <- toupper(as.character(aligned_bin))

  } else {
    # Profile alignment: add new sequences via MAFFT --add
    existing_names <- names(alignment)
    new_seqs       <- fasta[!names(fasta) %in% existing_names]
    n_new          <- length(new_seqs)

    if (verbose) message("Adding ", n_new, " new sequences to existing alignment")

    if (n_new > 0L) {
      # Write temp files for MAFFT --add
      tmp_new   <- tempfile(fileext = ".fasta")
      tmp_exist <- tempfile(fileext = ".fasta")
      tmp_out   <- tempfile(fileext = ".fasta")
      on.exit(unlink(c(tmp_new, tmp_exist, tmp_out)), add = TRUE)

      # Include ref in the existing alignment so we can trim from it
      aln_with_ref <- if (ref_name %in% existing_names) {
        alignment
      } else {
        c(ref, alignment)
      }

      Biostrings::writeXStringSet(new_seqs,    tmp_new)
      Biostrings::writeXStringSet(aln_with_ref, tmp_exist)

      exit_code <- system2(
        "mafft",
        args   = c("--quiet", "--add", tmp_new, tmp_exist),
        stdout = tmp_out
      )
      if (exit_code != 0L) stop("MAFFT --add alignment failed (exit code ", exit_code, ")")

      added_aln   <- Biostrings::readDNAStringSet(tmp_out)
      aligned_mat <- toupper(as.character(ape::as.DNAbin(added_aln)))

    } else {
      # No new sequences — use existing alignment directly (ensure ref is present)
      aln_with_ref <- if (ref_name %in% existing_names) {
        alignment
      } else {
        c(ref, alignment)
      }
      aligned_mat <- toupper(as.character(ape::as.DNAbin(aln_with_ref)))
    }
  }

  if (verbose) message("Alignment length: ", ncol(aligned_mat), " nt")

  # ── Trim to reference coordinates ───────────────────────────────────────────
  ref_row         <- aligned_mat[ref_name, ]
  valid_positions <- which(ref_row %in% c("A", "T", "C", "G"))

  if (verbose) {
    message("Trimming columns ", valid_positions[1], " to ",
            tail(valid_positions, 1), " (", length(valid_positions), " reference nt)")
  }

  aligned_trimmed <- aligned_mat[, valid_positions, drop = FALSE]

  # Remove reference row
  aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) != ref_name, , drop = FALSE]

  # Optionally drop sequences not in fasta
  if (drop && !is.null(alignment)) {
    keep_names      <- names(fasta)
    n_dropped       <- sum(!rownames(aligned_trimmed) %in% keep_names)
    if (verbose) message("Dropping ", n_dropped, " sequences not in fasta")
    aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) %in% keep_names, , drop = FALSE]
  }

  # Collapse each row to a sequence string
  seqs         <- apply(aligned_trimmed, 1, paste0, collapse = "")
  trimmed_set  <- Biostrings::DNAStringSet(seqs)
  names(trimmed_set) <- rownames(aligned_trimmed)

  # Validate coordinate invariant
  expected_len <- Biostrings::width(ref)
  actual_lens  <- unique(Biostrings::width(trimmed_set))
  if (!all(actual_lens == expected_len)) {
    warning("Trimmed length (", paste(actual_lens, collapse = ","),
            " nt) does not match reference length (", expected_len, " nt)")
  }

  if (verbose) {
    message("Output: ", length(trimmed_set), " sequences, ",
            expected_len, " nt (= reference length)")
  }

  trimmed_set
}
```

- [ ] **Step 4: Run tests — expect PASS (or SKIP if MAFFT absent)**

```bash
Rscript -e "devtools::test(filter = 'fasta_trim_ref')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 5 | PASS 0 ]` (no MAFFT) or `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]` (MAFFT present)

- [ ] **Step 5: Commit**

```bash
git add R/fasta_trim_ref.R tests/testthat/test-fasta_trim_ref.R
git commit -m "feat: add fasta_trim_ref() with MAFFT --add profile alignment"
```

---

## Task 8: define_lineage

**Files:**
- Create: `R/define_lineage.R`
- Create: `tests/testthat/test-define_lineage.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-define_lineage.R`:

```r
test_that("returns dataframe with strain and lineage columns", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_s3_class(result, "data.frame")
  expect_named(result, c("strain", "lineage"))
  expect_equal(nrow(result), length(make_test_biostring()))
})

test_that("strict match: NY10 requires both mutations", {
  # NY10_001 has pos10=A and pos25=C → matches NY10 (more specific than PARTIAL_A)
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("partial match assigns less-specific lineage, not unknown", {
  # NS2A_002 has pos10=A but pos25=T → matches PARTIAL_A (1 mutation) but not NY10 (2 mutations)
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "PARTIAL_A")
})

test_that("no match assigns unknown", {
  # OTHER_004 has no mutations
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|OTHER_004"], "unknown")
  # NS4B_003 has pos25=C but not pos10=A → neither lineage fully matched
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "unknown")
})

test_that("multiple-match: assigns most-specific (most mutations)", {
  # NY10_001 matches both NY10 (2 muts) and PARTIAL_A (1 mut) → NY10 wins
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("muts missing required columns throws error", {
  bad_muts <- tibble::tibble(lineage = "NY10", position = 10L, base = "A")
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), bad_muts),
    "lineage.*pos.*residue"
  )
})

test_that("invalid type throws error", {
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts(),
                   type = "protein"),
    "type must be"
  )
})
```

- [ ] **Step 2: Run tests — expect all 7 to FAIL**

```bash
Rscript -e "devtools::test(filter = 'define_lineage')"
```

- [ ] **Step 3: Implement define_lineage**

Create `R/define_lineage.R`:

```r
#' Assign lineage labels to sequences based on diagnostic mutations
#'
#' Requires the alignment to have been produced by fasta_trim_ref() so that
#' alignment position N equals ungapped reference position N. Positions in
#' `muts$pos` must use this same reference coordinate system.
#'
#' Assignment is strict: all mutations for a lineage must be present. If a
#' sequence matches multiple lineages, the one with the most required mutations
#' (most specific) is assigned. Ties are broken alphabetically with a warning.
#'
#' @param alignment DNAStringSet; pre-aligned sequences, all equal length.
#'   Must be output of fasta_trim_ref() to satisfy coordinate contract.
#' @param ref length-1 DNAStringSet; reference sequence. Used to locate the
#'   first ATG start codon when type = "aa".
#' @param muts dataframe with columns: lineage (chr), pos (int), residue (chr).
#'   pos is in reference coordinates (ungapped).
#' @param type "nuc" (default) or "aa"
#' @param verbose logical; if TRUE prints lineage count table (default TRUE)
#'
#' @return dataframe with columns strain (chr) and lineage (chr)
#' @export
define_lineage <- function(alignment, ref, muts, type = "nuc", verbose = TRUE) {

  required_cols <- c("lineage", "pos", "residue")
  if (!all(required_cols %in% names(muts))) {
    stop("muts must have columns: ", paste(required_cols, collapse = ", "))
  }

  if (!type %in% c("nuc", "aa")) stop("type must be 'nuc' or 'aa'")

  seq_names <- names(alignment)

  # ── Build residue lookup function ────────────────────────────────────────────
  if (type == "aa") {
    # Find first ATG in reference to establish CDS start
    atg_matches <- Biostrings::matchPattern("ATG", ref[[1]])
    if (length(atg_matches) == 0L) stop("No ATG start codon found in ref")
    cds_start <- Biostrings::start(atg_matches)[1]

    # Translate all sequences from CDS start
    cds_seqs   <- Biostrings::subseq(alignment, start = cds_start)
    translated <- Biostrings::translate(cds_seqs, if.fuzzy.codon = "X")

    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(translated[seq_idx])) return(NA_character_)
      as.character(Biostrings::subseq(translated[seq_idx], pos, pos))
    }

  } else {
    # nuc: alignment position = reference position (coordinate contract)
    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(alignment[seq_idx])) return(NA_character_)
      toupper(as.character(Biostrings::subseq(alignment[seq_idx], pos, pos)))
    }
  }

  # ── Group mutations by lineage ───────────────────────────────────────────────
  lineage_muts <- split(
    muts[, c("pos", "residue"), drop = FALSE],
    muts$lineage
  )

  # ── Assign lineage to each sequence ─────────────────────────────────────────
  assign_one <- function(seq_idx) {
    matched <- purrr::keep(names(lineage_muts), function(lin) {
      lin_muts <- lineage_muts[[lin]]
      all(purrr::map2_lgl(
        lin_muts$pos, lin_muts$residue,
        ~ !is.na(get_residue(seq_idx, .x)) &&
            get_residue(seq_idx, .x) == toupper(.y)
      ))
    })

    if (length(matched) == 0L) return("unknown")
    if (length(matched) == 1L) return(matched)

    # Multiple: pick the lineage with the most required mutations (most specific)
    n_muts      <- purrr::map_int(matched, ~ nrow(lineage_muts[[.x]]))
    max_n       <- max(n_muts)
    candidates  <- matched[n_muts == max_n]

    if (length(candidates) > 1L) {
      warning(
        "Sequence '", seq_names[seq_idx],
        "' matched lineages with equal specificity: ",
        paste(sort(candidates), collapse = ", "),
        ". Assigning first alphabetically."
      )
      candidates <- sort(candidates)[1]
    }

    candidates
  }

  result_lineages <- purrr::map_chr(seq_along(alignment), assign_one)

  result <- data.frame(
    strain  = seq_names,
    lineage = result_lineages,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("\nLineage assignments:")
    print(table(result$lineage))

    n_unknown <- sum(result$lineage == "unknown")
    if (n_unknown > 0L) {
      message("\n", n_unknown, " unknown sequences:")
      message(paste(result$strain[result$lineage == "unknown"], collapse = "\n"))
    }
  }

  result
}
```

- [ ] **Step 4: Run tests — expect all 7 to PASS**

```bash
Rscript -e "devtools::test(filter = 'define_lineage')"
```
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ]`

- [ ] **Step 5: Run full test suite**

```bash
Rscript -e "devtools::test()"
```
Expected: all tests PASS (fasta_trim_ref tests SKIP if MAFFT absent, otherwise PASS)

- [ ] **Step 6: Commit**

```bash
git add R/define_lineage.R tests/testthat/test-define_lineage.R
git commit -m "feat: add define_lineage() with strict lineage assignment"
```

---

## Self-Review

**Spec coverage check:**

| Spec section | Task | Status |
|---|---|---|
| `fasta_read` — single/multi path, `verbose` | Task 3 | ✓ |
| `fasta_combine` — concat, dupe warn, empty stop | Task 4 | ✓ |
| `fasta_nest` — delimiter split, NULL on no match, warnings | Task 5 | ✓ |
| `fasta_unnest` — drops NULL, returns list | Task 6 | ✓ |
| `fasta_trim_ref` — full align, `--add` profile, `drop`, MAFFT check | Task 7 | ✓ |
| `define_lineage` — strict match, `"unknown"`, most-specific tie-break, `type="aa"` | Task 8 | ✓ |
| WNV test dataset in `inst/extdata/` | Task 2 | ✓ |
| Shared test fixtures | Task 1 | ✓ |
| `verbose=TRUE` across all functions | All tasks | ✓ |
| `message()` not `cat()` | All tasks | ✓ |

**Placeholder scan:** No TBDs. All steps have complete code.

**Type consistency:** `make_test_biostring()` / `make_test_ref()` / `make_test_metadata()` / `make_test_muts()` defined in Task 1 and used by name consistently in Tasks 3–8. Function signatures match spec exactly.
