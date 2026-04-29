# extract_metadata() Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `extract_metadata(biostring, names, delim)` — parses pipe-delimited FASTA sequence headers into a tidy metadata tibble.

**Architecture:** Single pure function in `R/extract_metadata.R`. Validates inputs, splits each sequence name with `stringr::str_split_fixed()` into a character matrix, pads short headers with `NA` and warns, then returns a named tibble. TDD throughout: write failing tests first, then implement.

**Tech Stack:** R, Biostrings (`XStringSet`), stringr (`str_split_fixed`, `fixed`, `str_count`), tibble (`as_tibble`), methods (`is`), testthat 3

---

## File Map

| File | Action |
|------|--------|
| `R/extract_metadata.R` | **Create** — full function implementation with roxygen docs |
| `tests/testthat/test-extract_metadata.R` | **Create** — all unit tests |
| `man/extract_metadata.Rd` | **Auto-generated** by `devtools::document()` — do not edit manually |
| `NAMESPACE` | **Auto-updated** by `devtools::document()` |

---

## Task 1: Input validation — tests + implementation

**Files:**
- Create: `tests/testthat/test-extract_metadata.R`
- Create: `R/extract_metadata.R`

- [ ] **Step 1: Create the test file with failing validation tests**

Create `tests/testthat/test-extract_metadata.R` with this content:

```r
# test-extract_metadata.R — unit tests for extract_metadata()
# Shared fixture: make_test_biostring() from helper-fixtures.R
# Returns DNAStringSet with 6 sequences named "WNV|YEAR|STATE|STRAIN_ID"

test_that("errors when biostring is not an XStringSet", {
  expect_error(
    extract_metadata("not_a_biostring", names = c("virus", "year")),
    "biostring must be an XStringSet"
  )
})

test_that("errors when names is not a character vector", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = 1:3),
    "names must be a non-empty character vector"
  )
})

test_that("errors when names is empty", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = character(0)),
    "names must be a non-empty character vector"
  )
})

test_that("errors when names contains duplicates", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = c("virus", "virus", "year")),
    "names contains duplicate column names"
  )
})

test_that("errors when delim is not a single non-empty string", {
  bs <- make_test_biostring()
  expect_error(
    extract_metadata(bs, names = c("virus", "year"), delim = ""),
    "delim must be a single non-empty string"
  )
  expect_error(
    extract_metadata(bs, names = c("virus", "year"), delim = c("|", "/")),
    "delim must be a single non-empty string"
  )
})
```

- [ ] **Step 2: Run tests to confirm they all fail**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::test(filter = 'extract_metadata')"
```

Expected: 5 failures — `could not find function "extract_metadata"`

- [ ] **Step 3: Create R/extract_metadata.R with validation only**

Create `R/extract_metadata.R`:

```r
#' Extract metadata from FASTA sequence headers into a tibble
#'
#' Splits delimited metadata embedded in sequence header names (e.g. GISAID-style
#' `hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05`) into a tidy tibble with
#' one row per sequence and one column per field. Returns metadata only — pair
#' with \code{fasta_nest()} if you also need the sequences.
#'
#' @param biostring An \code{XStringSet} (e.g. \code{DNAStringSet}) whose
#'   \code{names()} contain delimited metadata fields.
#' @param names Character vector of column names for the output tibble. Length
#'   controls how many fields are extracted; extra header fields are silently
#'   dropped. The first element is conventionally the sequence ID or strain name.
#' @param delim Single string used as the field delimiter. Default \code{"|"} for
#'   GISAID-style headers. Treated as a fixed string, not a regular expression.
#'
#' @return A \code{tibble} with \code{length(names)} columns and one row per
#'   sequence. Columns are typed as \code{character}; missing fields are
#'   \code{NA_character_}.
#'
#' @export
#'
#' @examples
#' library(Biostrings)
#' seqs <- DNAStringSet(c(
#'   "hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05" = "ATCG",
#'   "hCoV-19/USA/CA-CDPH/2021|EPI_ISL_999999|2021-06-15" = "GCTA"
#' ))
#' extract_metadata(seqs, names = c("strain", "accession", "date"), delim = "|")
extract_metadata <- function(biostring, names, delim = "|") {

  # ── Input validation ─────────────────────────────────────────────────────────

  # biostring must be a Biostrings XStringSet (DNAStringSet, AAStringSet, etc.)
  if (!methods::is(biostring, "XStringSet")) {
    stop("biostring must be an XStringSet (e.g. DNAStringSet)")
  }

  # names must be a non-empty character vector — controls output column count
  if (!is.character(names) || length(names) == 0L) {
    stop("names must be a non-empty character vector")
  }

  # Duplicate column names would produce an ambiguous tibble
  dupe_names <- names[duplicated(names)]
  if (length(dupe_names) > 0L) {
    stop("names contains duplicate column names: ", paste(dupe_names, collapse = ", "))
  }

  # delim must be a single non-empty string — passed to stringr::fixed()
  if (!is.character(delim) || length(delim) != 1L || !nzchar(delim)) {
    stop("delim must be a single non-empty string")
  }

  # ── Core logic ───────────────────────────────────────────────────────────────

  # Placeholder — implemented in Task 2
  NULL
}
```

- [ ] **Step 4: Load package and run tests to confirm validation tests pass**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::load_all('.'); devtools::test(filter = 'extract_metadata')"
```

Expected: 5 PASS (validation), 0 FAIL. The function returns NULL for valid input — that is intentional at this stage.

- [ ] **Step 5: Commit**

```bash
git add R/extract_metadata.R tests/testthat/test-extract_metadata.R
git commit -m "feat: scaffold extract_metadata with input validation and tests"
```

---

## Task 2: Core extraction — happy path

**Files:**
- Modify: `tests/testthat/test-extract_metadata.R` — append new tests
- Modify: `R/extract_metadata.R` — add core logic after the validation block

- [ ] **Step 1: Append happy-path tests to test-extract_metadata.R**

Open `tests/testthat/test-extract_metadata.R` and append:

```r
test_that("returns a tibble with correct dimensions for exact field count", {
  # make_test_biostring() names: "WNV|YEAR|STATE|STRAIN_ID" (4 pipe-fields)
  bs     <- make_test_biostring()
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 6L)   # 6 sequences in fixture
  expect_equal(ncol(result), 4L)
  expect_named(result, c("virus", "year", "state", "strain_id"))
})

test_that("extracts correct values from pipe-delimited headers", {
  bs     <- make_test_biostring()
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  # First sequence: "WNV|2021|CO|NY10_001"
  expect_equal(result$virus[1],     "WNV")
  expect_equal(result$year[1],      "2021")
  expect_equal(result$state[1],     "CO")
  expect_equal(result$strain_id[1], "NY10_001")
})

test_that("extra header fields beyond length(names) are silently dropped", {
  bs     <- make_test_biostring()
  # Request only 2 columns — 4th and 3rd fields are dropped
  result <- suppressWarnings(
    extract_metadata(bs, names = c("virus", "year"), delim = "|")
  )

  expect_equal(ncol(result), 2L)
  expect_named(result, c("virus", "year"))
  expect_equal(result$virus[1], "WNV")
  expect_equal(result$year[1],  "2021")
})

test_that("works with a single sequence", {
  bs     <- make_test_biostring()[1]   # subset to 1 sequence
  result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")

  expect_equal(nrow(result), 1L)
  expect_equal(result$strain_id, "NY10_001")
})

test_that("works with a custom delimiter", {
  # Build a DNAStringSet with slash-delimited names
  bs <- Biostrings::DNAStringSet(c(
    "WNV/2021/CO" = "ATCG",
    "WNV/2020/WY" = "GCTA"
  ))
  result <- extract_metadata(bs, names = c("virus", "year", "state"), delim = "/")

  expect_equal(result$state, c("CO", "WY"))
})
```

- [ ] **Step 2: Run to confirm new tests fail**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::load_all('.'); devtools::test(filter = 'extract_metadata')"
```

Expected: 5 PASS (validation), 5 FAIL (happy path — function returns NULL)

- [ ] **Step 3: Replace the placeholder NULL in R/extract_metadata.R with core logic**

Find the `# ── Core logic ───` section (after the validation block) and replace everything from that comment to the end of the function with:

```r
  # ── Core logic ───────────────────────────────────────────────────────────────

  # Pull raw header strings from the XStringSet names slot
  headers <- base::names(biostring)

  # Split each header into exactly length(names) fields.
  # str_split_fixed returns a [n_sequences × length(names)] character matrix.
  # - Headers with MORE fields: extra fields are silently dropped.
  # - Headers with FEWER fields: missing positions are padded with "".
  mat <- stringr::str_split_fixed(
    string  = headers,
    pattern = stringr::fixed(delim),
    n       = length(names)
  )

  # Detect short headers: count delimiters; fewer than (length(names) - 1) means
  # the header cannot fill all requested fields.
  n_delims  <- stringr::str_count(headers, stringr::fixed(delim))
  short_idx <- which(n_delims < (length(names) - 1L))

  if (length(short_idx) > 0L) {
    warning(
      "The following sequence headers have fewer fields than `names` ",
      "— missing fields set to NA:\n",
      paste(headers[short_idx], collapse = "\n"),
      call. = FALSE
    )
  }

  # Convert "" padding produced by str_split_fixed for short headers to NA
  mat[mat == ""] <- NA_character_

  # Coerce to tibble and apply caller-supplied column names
  result           <- tibble::as_tibble(mat, .name_repair = "minimal")
  colnames(result) <- names

  result
```

- [ ] **Step 4: Run tests to confirm all pass**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::load_all('.'); devtools::test(filter = 'extract_metadata')"
```

Expected: 10 PASS, 0 FAIL

- [ ] **Step 5: Commit**

```bash
git add R/extract_metadata.R tests/testthat/test-extract_metadata.R
git commit -m "feat: implement extract_metadata core extraction logic"
```

---

## Task 3: Short-header edge cases

**Files:**
- Modify: `tests/testthat/test-extract_metadata.R` — append short-header tests
- No changes to `R/extract_metadata.R` — logic is already in place

- [ ] **Step 1: Append short-header tests**

Append to `tests/testthat/test-extract_metadata.R`:

```r
test_that("pads with NA and warns when a header has fewer fields than names", {
  # Build a set where one sequence is missing the 4th pipe-field
  bs <- Biostrings::DNAStringSet(c(
    "WNV|2021|CO|NY10_001" = "ATCG",   # 4 fields — OK
    "WNV|2020|WY"          = "GCTA"    # 3 fields — short by 1
  ))

  expect_warning(
    result <- extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|"),
    "fewer fields than `names`"
  )

  # Short row: strain_id should be NA
  expect_equal(result$strain_id[1], "NY10_001")
  expect_true(is.na(result$strain_id[2]))
})

test_that("warning message lists the affected sequence names", {
  bs <- Biostrings::DNAStringSet(c(
    "WNV|2021|CO|NY10_001" = "ATCG",
    "WNV|2020|WY"          = "GCTA"
  ))

  expect_warning(
    extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|"),
    "WNV\\|2020\\|WY"   # the short header should appear in the warning text
  )
})

test_that("no warning when all headers have exactly the right number of fields", {
  bs <- make_test_biostring()   # all 6 sequences have 4 pipe-fields
  expect_no_warning(
    extract_metadata(bs, names = c("virus", "year", "state", "strain_id"), delim = "|")
  )
})

test_that("no warning when headers have more fields than names", {
  bs <- make_test_biostring()   # 4 pipe-fields; request only 2
  expect_no_warning(
    extract_metadata(bs, names = c("virus", "year"), delim = "|")
  )
})
```

- [ ] **Step 2: Run to confirm new tests pass (logic already handles this)**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::load_all('.'); devtools::test(filter = 'extract_metadata')"
```

Expected: 14 PASS, 0 FAIL

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-extract_metadata.R
git commit -m "test: add short-header edge-case tests for extract_metadata"
```

---

## Task 4: Documentation and export

**Files:**
- No source edits — roxygen block already written in Task 1
- Run `devtools::document()` to regenerate `man/` and `NAMESPACE`

- [ ] **Step 1: Regenerate docs**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::document()"
```

Expected output includes:
```
ℹ Updating tiphydr documentation
Writing extract_metadata.Rd
Writing NAMESPACE
```

- [ ] **Step 2: Verify man page was created**

```bash
ls /Users/user/Programming_Directory/tiphydr/man/extract_metadata.Rd
```

Expected: file exists

- [ ] **Step 3: Run full test suite to confirm nothing regressed**

```bash
cd /Users/user/Programming_Directory/tiphydr
Rscript -e "devtools::test()"
```

Expected: all tests pass, 0 failures

- [ ] **Step 4: Commit**

```bash
git add R/extract_metadata.R man/extract_metadata.Rd NAMESPACE
git commit -m "docs: add roxygen documentation and export for extract_metadata"
```

---

## Self-review checklist

- [x] **Spec coverage:** All 8 test scenarios from the spec are covered across Tasks 1–3
- [x] **No placeholders:** All steps include complete code
- [x] **Type consistency:** `extract_metadata` signature is identical across all tasks; `make_test_biostring()` fixture used consistently
- [x] **Validation errors:** All 4 spec-required error messages match exactly
- [x] **Warning behaviour:** Warns with short header names, no warning for exact/extra field counts
- [x] **`base::names()` shadowing:** Called as `base::names(biostring)` before the `names` argument is used in logic, avoiding the shadowing issue
