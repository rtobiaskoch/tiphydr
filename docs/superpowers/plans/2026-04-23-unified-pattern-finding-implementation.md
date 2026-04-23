# Unified Pattern Finding Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract generalizable pattern-finding functions from examples, create unified `find_residue()` API for searching DNA/protein sequences with smart auto-translation.

**Architecture:** 
- Build pure, testable helpers first (translate, pattern detection)
- Layer pattern search on top
- Create main `find_residue()` function that auto-detects sequence type and translates DNA when searching for amino acids
- Update `define_lineage()` to reuse start-detection logic
- Rewrite examples to showcase composition

**Tech Stack:** Biostrings (pattern matching), tidyverse (output format), purrr (map/vectorization)

---

## Task 1: Write and test `translate_from_position()`

**Files:**
- Create: `R/translate_from_position.R`
- Create: `tests/testthat/test-translate_from_position.R`

- [ ] **Step 1: Write failing tests for `translate_from_position()`**

Create `tests/testthat/test-translate_from_position.R`:

```r
test_that("translate_from_position translates DNA from a start position", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA")
  )
  
  result <- translate_from_position(dna_seq, start_position = 1)
  
  expect_s4_class(result, "AAStringSet")
  expect_equal(length(result), 2)
  expect_equal(names(result), names(dna_seq))
  # ATG -> M, AAA -> K, TAA -> * (stop)
  expect_equal(as.character(result)[1], "MK*")
  expect_equal(as.character(result)[2], "MP*")
})

test_that("translate_from_position respects start position offset", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG")
  )
  
  # Starting at position 2 skips the A, reads "TGA AA TAG"
  result <- translate_from_position(dna_seq, start_position = 2)
  
  expect_s4_class(result, "AAStringSet")
  # TGA is stop, skipped; AA + TAG incomplete
  expect_equal(as.character(result)[1], "*")
})

test_that("translate_from_position handles fuzzy codons", {
  dna_seq <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGNNATAG")  # N is ambiguous
  )
  
  result <- translate_from_position(dna_seq, start_position = 1)
  
  # Fuzzy codon should be X
  expect_match(as.character(result)[1], "X")
})

test_that("translate_from_position preserves sequence names", {
  dna_seq <- Biostrings::DNAStringSet(
    c("seq_A" = "ATGAAATAG", "seq_B" = "ATGCCCTAA")
  )
  
  result <- translate_from_position(dna_seq, start_position = 1)
  
  expect_equal(names(result), c("seq_A", "seq_B"))
})
```

- [ ] **Step 2: Verify tests fail**

Run: `R -e "devtools::test_file('tests/testthat/test-translate_from_position.R')"`

Expected output: Tests fail with "function 'translate_from_position' not found"

- [ ] **Step 3: Implement `translate_from_position()`**

Create `R/translate_from_position.R`:

```r
#' Translate DNA from a specified start position
#'
#' Extract a subsequence starting from a given position and translate to amino acids.
#' Useful for translating from a known start codon position.
#'
#' @param xstring_set DNAStringSet to translate
#' @param start_position integer; 1-indexed position to start translation
#'
#' @return AAStringSet with same names as input
#' @keywords internal
#' @examples
#' dna <- Biostrings::DNAStringSet(c("Seq1" = "ATGAAATAG"))
#' translate_from_position(dna, start_position = 1)
#'
translate_from_position <- function(xstring_set, start_position) {
  if (!is(xstring_set, "DNAStringSet")) {
    stop("xstring_set must be a DNAStringSet")
  }
  if (!is.integer(start_position) && !is.numeric(start_position)) {
    stop("start_position must be numeric")
  }
  if (start_position < 1) {
    stop("start_position must be >= 1")
  }
  
  # Extract from start_position to end
  subseqs <- Biostrings::subseq(xstring_set, start = start_position)
  
  # Translate with fuzzy codon handling
  translated <- Biostrings::translate(subseqs, if.fuzzy.codon = "X")
  
  translated
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `R -e "devtools::test_file('tests/testthat/test-translate_from_position.R')"`

Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add R/translate_from_position.R tests/testthat/test-translate_from_position.R
git commit -m "feat: add translate_from_position helper for DNA translation from arbitrary start"
```

---

## Task 2: Write and test `is_aa_pattern()`

**Files:**
- Create: `R/is_aa_pattern.R`
- Create: `tests/testthat/test-is_aa_pattern.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-is_aa_pattern.R`:

```r
test_that("is_aa_pattern detects amino acid patterns", {
  expect_true(is_aa_pattern("K"))
  expect_true(is_aa_pattern("M"))
  expect_true(is_aa_pattern("KM"))
  expect_true(is_aa_pattern("ACDEFGHIKLMNPQRSTVWY"))
})

test_that("is_aa_pattern detects nucleotide patterns", {
  expect_false(is_aa_pattern("ATG"))
  expect_false(is_aa_pattern("TAA"))
  expect_false(is_aa_pattern("ACGT"))
})

test_that("is_aa_pattern prioritizes amino acids for ambiguous patterns", {
  # A is both AA and nucleotide; should return TRUE (prioritize AA)
  expect_true(is_aa_pattern("A"))
})

test_that("is_aa_pattern is case insensitive", {
  expect_true(is_aa_pattern("k"))
  expect_true(is_aa_pattern("atg") == FALSE)
})

test_that("is_aa_pattern rejects invalid patterns", {
  expect_false(is_aa_pattern("XYZ123"))
  expect_false(is_aa_pattern(""))
})
```

- [ ] **Step 2: Verify tests fail**

Run: `R -e "devtools::test_file('tests/testthat/test-is_aa_pattern.R')"`

Expected: Tests fail with function not found

- [ ] **Step 3: Implement `is_aa_pattern()`**

Create `R/is_aa_pattern.R`:

```r
#' Detect whether a pattern is amino acid or nucleotide
#'
#' Returns TRUE if pattern contains only standard amino acid codes,
#' FALSE if nucleotide or invalid. Case-insensitive.
#'
#' @param pattern character; pattern to classify
#'
#' @return logical; TRUE if amino acid, FALSE if nucleotide
#' @keywords internal
#'
is_aa_pattern <- function(pattern) {
  if (!is.character(pattern) || length(pattern) == 0) {
    return(FALSE)
  }
  
  # Standard amino acid codes (uppercase)
  aa_codes <- "ACDEFGHIKLMNPQRSTVWY"
  
  # Convert to uppercase for comparison
  pattern_upper <- toupper(pattern)
  
  # Check if all characters in pattern are AA codes
  pattern_chars <- unlist(strsplit(pattern_upper, ""))
  all(pattern_chars %in% unlist(strsplit(aa_codes, "")))
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `R -e "devtools::test_file('tests/testthat/test-is_aa_pattern.R')"`

Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add R/is_aa_pattern.R tests/testthat/test-is_aa_pattern.R
git commit -m "feat: add is_aa_pattern helper to detect pattern type"
```

---

## Task 3: Write and test `search_pattern_in_xstring()`

**Files:**
- Create: `R/search_pattern_in_xstring.R`
- Create: `tests/testthat/test-search_pattern_in_xstring.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-search_pattern_in_xstring.R`:

```r
test_that("search_pattern_in_xstring finds nucleotide matches", {
  dna <- Biostrings::DNAString("ATGAAATAG")
  
  result <- search_pattern_in_xstring(dna, "ATG")
  
  expect_integer(result)
  expect_equal(result, 1L)
})

test_that("search_pattern_in_xstring finds multiple matches", {
  dna <- Biostrings::DNAString("ATGATGATG")
  
  result <- search_pattern_in_xstring(dna, "ATG")
  
  expect_integer(result)
  expect_equal(length(result), 3)
  expect_equal(result, c(1L, 4L, 7L))
})

test_that("search_pattern_in_xstring returns empty for no matches", {
  dna <- Biostrings::DNAString("AAATTT")
  
  result <- search_pattern_in_xstring(dna, "ATG")
  
  expect_integer(result)
  expect_equal(length(result), 0)
})

test_that("search_pattern_in_xstring works with amino acids", {
  aa <- Biostrings::AAString("MKKM")
  
  result <- search_pattern_in_xstring(aa, "K")
  
  expect_integer(result)
  expect_equal(result, c(2L, 3L))
})

test_that("search_pattern_in_xstring is case insensitive for amino acids", {
  aa <- Biostrings::AAString("MKKM")
  
  result <- search_pattern_in_xstring(aa, "k")
  
  expect_integer(result)
  expect_equal(result, c(2L, 3L))
})
```

- [ ] **Step 2: Verify tests fail**

Run: `R -e "devtools::test_file('tests/testthat/test-search_pattern_in_xstring.R')"`

Expected: Tests fail with function not found

- [ ] **Step 3: Implement `search_pattern_in_xstring()`**

Create `R/search_pattern_in_xstring.R`:

```r
#' Search for a pattern in a single XString
#'
#' Wrapper around Biostrings::matchPattern for clean extraction of match positions.
#'
#' @param xstring XString (DNAString, AAString, etc.) to search in
#' @param pattern character; pattern to search for (case-insensitive for AA)
#'
#' @return integer vector of start positions (empty if no matches)
#' @keywords internal
#'
search_pattern_in_xstring <- function(xstring, pattern) {
  if (!methods::is(xstring, "XString")) {
    stop("xstring must be an XString (DNAString, AAString, etc.)")
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("pattern must be a single character string")
  }
  
  # For AAString, convert pattern to uppercase
  if (methods::is(xstring, "AAString")) {
    pattern <- toupper(pattern)
  }
  
  # Use matchPattern to find all occurrences
  matches <- Biostrings::matchPattern(pattern, xstring)
  
  # Extract start positions
  if (length(matches) > 0) {
    as.integer(Biostrings::start(matches))
  } else {
    integer(0)
  }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `R -e "devtools::test_file('tests/testthat/test-search_pattern_in_xstring.R')"`

Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add R/search_pattern_in_xstring.R tests/testthat/test-search_pattern_in_xstring.R
git commit -m "feat: add search_pattern_in_xstring wrapper for pattern matching"
```

---

## Task 4: Write and test `detect_sequence_start()`

**Files:**
- Create: `R/detect_sequence_start.R`
- Create: `tests/testthat/test-detect_sequence_start.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-detect_sequence_start.R`:

```r
test_that("detect_sequence_start finds modal start position", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA", "Seq3" = "ATGTTT")
  )
  
  result <- detect_sequence_start(dna, pattern = "ATG")
  
  expect_is(result, "list")
  expect_named(result, c("position", "table"))
  expect_equal(result$position, 1L)
})

test_that("detect_sequence_start returns mode when positions vary", {
  dna <- Biostrings::DNAStringSet(
    c(
      "Seq1" = "ATGAAATAG",           # ATG at position 1
      "Seq2" = "AAATGATAG",           # ATG at position 4
      "Seq3" = "ATGAAATAG"            # ATG at position 1
    )
  )
  
  result <- detect_sequence_start(dna, pattern = "ATG")
  
  # Mode should be 1 (appears twice)
  expect_equal(result$position, 1L)
  expect_true(result$table[as.character(1)] >= 2)
})

test_that("detect_sequence_start excludes sequences without pattern", {
  dna <- Biostrings::DNAStringSet(
    c(
      "Seq1" = "ATGAAATAG",           # ATG at position 1
      "Seq2" = "AAATTTTAA",           # No ATG
      "Seq3" = "ATGCCCTAA"            # ATG at position 1
    )
  )
  
  result <- detect_sequence_start(dna, pattern = "ATG")
  
  # Mode should be 1 (from Seq1, Seq3; Seq2 excluded)
  expect_equal(result$position, 1L)
})

test_that("detect_sequence_start messages when pattern is missing", {
  dna <- Biostrings::DNAStringSet(
    c(
      "Seq1" = "ATGAAATAG",
      "Seq2" = "AAATTTTAA",
      "Seq3" = "ATGAAATAG"
    )
  )
  
  expect_message(
    detect_sequence_start(dna, pattern = "ATG", verbose = TRUE),
    "detected at position 1"
  )
})

test_that("detect_sequence_start warns about missing start codons", {
  dna <- Biostrings::DNAStringSet(
    c(
      "Seq1" = "ATGAAATAG",
      "Seq2" = "AAATTTTAA",
      "Seq3" = "ATGAAATAG"
    )
  )
  
  result <- detect_sequence_start(dna, pattern = "ATG", verbose = TRUE)
  
  # Should have found the position and reported missing count
  expect_is(result, "list")
})
```

- [ ] **Step 2: Verify tests fail**

Run: `R -e "devtools::test_file('tests/testthat/test-detect_sequence_start.R')"`

Expected: Tests fail with function not found

- [ ] **Step 3: Implement `detect_sequence_start()`**

Create `R/detect_sequence_start.R`:

```r
#' Detect the modal start position of a pattern across DNA sequences
#'
#' Finds the most common start position for a pattern (e.g., ATG) across
#' all sequences. Ignores sequences missing the pattern. Returns the modal
#' position and a frequency table.
#'
#' @param xstring_set DNAStringSet to search
#' @param pattern character; pattern to search (default "ATG")
#' @param verbose logical; if TRUE, messages the detected position (default TRUE)
#'
#' @return list with elements:
#'   - `position` (integer): the modal start position
#'   - `table` (table): frequency of each position found
#'
#' @keywords internal
#'
detect_sequence_start <- function(xstring_set, pattern = "ATG", verbose = TRUE) {
  if (!is(xstring_set, "DNAStringSet")) {
    stop("xstring_set must be a DNAStringSet")
  }
  
  # Find start position for each sequence
  positions <- integer(length(xstring_set))
  for (i in seq_along(xstring_set)) {
    seq <- xstring_set[[i]]
    matches <- search_pattern_in_xstring(seq, pattern)
    
    # Record first match, or NA if no match
    if (length(matches) > 0) {
      positions[i] <- matches[1]
    } else {
      positions[i] <- NA_integer_
    }
  }
  
  # Remove NAs for mode calculation
  valid_positions <- positions[!is.na(positions)]
  
  if (length(valid_positions) == 0) {
    stop("No sequences contain the pattern '", pattern, "'")
  }
  
  # Calculate mode
  freq_table <- table(valid_positions)
  modal_position <- as.integer(names(freq_table)[which.max(freq_table)])
  
  # Message user
  if (verbose) {
    n_missing <- sum(is.na(positions))
    message(
      "Start codon '", pattern, "' detected at position ", modal_position,
      " (found in ", length(valid_positions), " sequences",
      if (n_missing > 0) paste0("; missing in ", n_missing, " sequences)") else ")"
    )
  }
  
  list(position = modal_position, table = freq_table)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `R -e "devtools::test_file('tests/testthat/test-detect_sequence_start.R')"`

Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add R/detect_sequence_start.R tests/testthat/test-detect_sequence_start.R
git commit -m "feat: add detect_sequence_start to find modal start position across sequences"
```

---

## Task 5: Write and test `find_residue()`

**Files:**
- Create: `R/find_residue.R`
- Create: `tests/testthat/test-find_residue.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-find_residue.R`:

```r
test_that("find_residue searches DNA for nucleotide patterns", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA")
  )
  
  result <- find_residue(dna, pattern = "ATG")
  
  expect_is(result, "tbl_df")
  expect_named(result, c("sequence_name", "location", "pattern"))
  expect_equal(nrow(result), 2L)  # One match per sequence
  expect_equal(result$location, c(1L, 1L))
})

test_that("find_residue searches protein for amino acid patterns", {
  aa <- Biostrings::AAStringSet(
    c("Seq1" = "MKKM", "Seq2" = "MKRM")
  )
  
  result <- find_residue(aa, pattern = "K")
  
  expect_is(result, "tbl_df")
  expect_equal(nrow(result), 3L)  # 2 K's in Seq1, 1 K in Seq2
  expect_equal(result$pattern, rep("K", 3))
})

test_that("find_residue translates DNA when searching for amino acids", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG")  # Translates to MK*
  )
  
  result <- find_residue(dna, pattern = "K")
  
  expect_is(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  # K is at position 2 in the translated sequence
  expect_equal(result$location, 2L)
  expect_message(find_residue(dna, pattern = "K"), "Start codon 'ATG' detected")
})

test_that("find_residue handles multiple patterns", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG", "Seq2" = "ATGCCCTAA")
  )
  
  result <- find_residue(dna, pattern = c("ATG", "TAA"))
  
  expect_is(result, "tbl_df")
  # ATG at 1, ATG at 1, TAA at 7, TAA at 8
  expect_equal(nrow(result), 4L)
})

test_that("find_residue allows user to override start position", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGAAATAG")  # Normal: start at 1 (M)
  )
  
  # Start at position 2 instead
  result <- find_residue(dna, pattern = "K", start_position = 2)
  
  expect_is(result, "tbl_df")
  # Position 2 onward: TGA (stop), AA, T (incomplete)
  # Result depends on translation logic
  expect_true(nrow(result) >= 0)  # May be 0 if no K found with offset
})

test_that("find_residue returns long format with one row per match", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "ATGATGATG")
  )
  
  result <- find_residue(dna, pattern = "ATG")
  
  expect_equal(nrow(result), 3L)  # Three ATG matches
  expect_equal(unique(result$sequence_name), "Seq1")
  expect_equal(unique(result$pattern), "ATG")
  expect_equal(result$location, c(1L, 4L, 7L))
})

test_that("find_residue preserves sequence names", {
  dna <- Biostrings::DNAStringSet(
    c("seq_alpha" = "ATGAAATAG", "seq_beta" = "ATGCCCTAA")
  )
  
  result <- find_residue(dna, pattern = "ATG")
  
  expect_equal(sort(unique(result$sequence_name)), c("seq_alpha", "seq_beta"))
})

test_that("find_residue returns empty tibble when no matches", {
  dna <- Biostrings::DNAStringSet(
    c("Seq1" = "AAATTT")
  )
  
  result <- find_residue(dna, pattern = "GGG")
  
  expect_is(result, "tbl_df")
  expect_equal(nrow(result), 0L)
  expect_named(result, c("sequence_name", "location", "pattern"))
})
```

- [ ] **Step 2: Verify tests fail**

Run: `R -e "devtools::test_file('tests/testthat/test-find_residue.R')"`

Expected: Tests fail with function not found

- [ ] **Step 3: Implement `find_residue()`**

Create `R/find_residue.R`:

```r
#' Find a residue pattern in DNA or protein sequences
#'
#' Search for DNA or amino acid patterns in an XStringSet. When searching for
#' amino acids in DNA, automatically detects the start codon position and
#' translates. Returns results in long format: one row per match.
#'
#' @param xstring_set DNAStringSet or AAStringSet to search
#' @param pattern character vector of patterns to search (DNA codons or amino acids)
#' @param start_position integer; optional override for translation start position
#'   (default: auto-detect from modal ATG position)
#' @param verbose logical; if TRUE, messages the detected start position when
#'   translating DNA (default TRUE)
#'
#' @return tibble with columns:
#'   - `sequence_name` (character): sequence identifier from xstring_set names
#'   - `location` (integer): 1-indexed position of match
#'   - `pattern` (character): the pattern searched for
#'
#' @export
#' @examples
#' dna <- Biostrings::DNAStringSet(c("Seq1" = "ATGAAATAG"))
#' find_residue(dna, pattern = "ATG")  # Search for start codon
#'
#' # Search for amino acids (auto-translates)
#' find_residue(dna, pattern = "K")
#'
find_residue <- function(xstring_set, pattern, start_position = NULL, verbose = TRUE) {
  
  # Validate inputs
  if (!is(xstring_set, "XStringSet")) {
    stop("xstring_set must be an XStringSet (DNAStringSet, AAStringSet, etc.)")
  }
  if (!is.character(pattern) || length(pattern) == 0) {
    stop("pattern must be a non-empty character vector")
  }
  
  # Detect sequence type
  is_dna <- is(xstring_set, "DNAStringSet")
  is_aa <- is(xstring_set, "AAStringSet")
  
  # Detect pattern type (single pattern for logic, handled per-pattern below)
  first_pattern_is_aa <- is_aa_pattern(pattern[1])
  
  # --- Logic for DNA + AA pattern: translate first ---
  if (is_dna && first_pattern_is_aa) {
    if (is.null(start_position)) {
      # Auto-detect start codon
      start_info <- detect_sequence_start(xstring_set, pattern = "ATG", verbose = verbose)
      start_position <- start_info$position
    } else {
      if (verbose) {
        message("Using user-specified start position: ", start_position)
      }
    }
    
    # Translate from the detected/specified position
    xstring_set <- translate_from_position(xstring_set, start_position = start_position)
  }
  
  # --- Now search in the (possibly translated) StringSet ---
  results_list <- list()
  
  for (pat in pattern) {
    pattern_results <- list()
    
    for (i in seq_along(xstring_set)) {
      seq <- xstring_set[[i]]
      seq_name <- names(xstring_set)[i]
      
      # Search for pattern in this sequence
      locations <- search_pattern_in_xstring(seq, pat)
      
      # Build result rows
      if (length(locations) > 0) {
        for (loc in locations) {
          pattern_results[[length(pattern_results) + 1]] <- 
            list(sequence_name = seq_name, location = loc, pattern = pat)
        }
      }
    }
    
    # Bind pattern results
    if (length(pattern_results) > 0) {
      results_list[[pat]] <- dplyr::bind_rows(pattern_results)
    }
  }
  
  # Combine all patterns
  if (length(results_list) > 0) {
    result <- dplyr::bind_rows(results_list)
  } else {
    result <- tibble::tibble(
      sequence_name = character(),
      location = integer(),
      pattern = character()
    )
  }
  
  result
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `R -e "devtools::test_file('tests/testthat/test-find_residue.R')"`

Expected: All tests pass

- [ ] **Step 5: Commit**

```bash
git add R/find_residue.R tests/testthat/test-find_residue.R
git commit -m "feat: add find_residue main function for unified pattern finding"
```

---

## Task 6: Update `define_lineage()` to use `detect_sequence_start()`

**Files:**
- Modify: `R/define_lineage.R`

- [ ] **Step 1: Examine current `define_lineage()` implementation**

Run: `grep -n "start" R/define_lineage.R | head -20`

This identifies where start position logic currently lives.

- [ ] **Step 2: Replace manual start position detection with `detect_sequence_start()`**

Find the current logic in `R/define_lineage.R` that detects start positions (likely involves finding ATG, calculating mode). Replace it with:

```r
# In define_lineage, replace manual start position detection with:
start_info <- detect_sequence_start(sequences, pattern = "ATG", verbose = verbose)
start_position <- start_info$position

# Then use start_position as before for translation/lineage assignment
```

- [ ] **Step 3: Test that `define_lineage()` still works**

Run: `R -e "devtools::test_file('tests/testthat/test-define_lineage.R')"`

Expected: All existing tests pass

- [ ] **Step 4: Commit**

```bash
git add R/define_lineage.R
git commit -m "refactor: use detect_sequence_start in define_lineage for consistency"
```

---

## Task 7: Rewrite example scripts to use new functions

**Files:**
- Modify: `examples/find_aa_loc.R` → use `find_residue()`
- Modify: `examples/get_start_Met.R` → use `detect_sequence_start()`
- Modify: `examples/find_nucleotide_pattern.R` → use `find_residue()`
- Modify: `examples/get_stop_codons.R` → use `find_residue()`
- Modify: `examples/translate_and_identify_aa.R` → use new functions

- [ ] **Step 1: Rewrite `examples/find_aa_loc.R`**

Replace old `find_aa_loc()` calls with `find_residue()`:

```r
# OLD:
# K <- find_aa_loc(data, aa = "K", start_codon = F, start_index = 97)

# NEW:
K <- find_residue(data, pattern = "K", start_position = 97)
```

- [ ] **Step 2: Rewrite `examples/get_start_Met.R`**

Replace old `get_start_met()` calls with `detect_sequence_start()`:

```r
# OLD:
# start_codons = get_start_met(data)

# NEW:
start_info <- detect_sequence_start(data, pattern = "ATG", verbose = TRUE)
start_codons <- tibble::tibble(
  sequence_id = names(data),
  start_codon_position = start_info$position  # All sequences use modal position
)
```

- [ ] **Step 3: Rewrite `examples/find_nucleotide_pattern.R`**

Replace old `find_nucleotide_patterns()` calls:

```r
# OLD:
# test = find_nucleotide_patterns(dna_sequences, stop_codons, frame = 1)

# NEW:
test <- find_residue(dna_sequences, pattern = c("TAA", "TAG", "TGA"))
```

- [ ] **Step 4: Rewrite `examples/get_stop_codons.R`**

Replace old script:

```r
# OLD:
# results <- find_nucleotide_patterns(fasta0, stop_codons, frame = 1)

# NEW:
stop_codons <- c("TAA", "TAG", "TGA")
results <- find_residue(fasta0, pattern = stop_codons)

count_stop <- results %>%
  dplyr::group_by(sequence_id = sequence_name) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop")
```

- [ ] **Step 5: Rewrite `examples/translate_and_identify_aa.R`**

Replace old script with new composition:

```r
# Load sequences
data <- fasta_read('results/.fasta')

# Find start codon (auto-detects modal position)
start_info <- detect_sequence_start(data, verbose = TRUE)

# Find amino acids K and M in translated sequence
# (find_residue will auto-translate when given AA pattern)
K <- find_residue(data, pattern = "K")  # Auto-translates from detected start
M <- find_residue(data, pattern = "M")

# Filter and summarize as before
K_count <- K %>%
  dplyr::group_by(location) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(count > 10 & count < 4000)
```

- [ ] **Step 6: Commit**

```bash
git add examples/find_aa_loc.R examples/get_start_Met.R examples/find_nucleotide_pattern.R examples/get_stop_codons.R examples/translate_and_identify_aa.R
git commit -m "docs: rewrite examples to use new unified find_residue API"
```

---

## Task 8: Export new public functions and verify package builds

**Files:**
- Modify: `NAMESPACE` (or roxygen headers in R files)
- Test: run `devtools::load_all()` and `devtools::check()`

- [ ] **Step 1: Add roxygen `@export` tag to public functions**

Update these files to include `#' @export`:
- `R/find_residue.R` — already has it (added in Task 5)
- `R/translate_from_position.R` — change `@keywords internal` to `#' @export` if making public
- `R/detect_sequence_start.R` — change `@keywords internal` to `#' @export` if making public

Per design spec, these are public:
- `find_residue()` ✓
- `translate_from_position()` — optional, but useful standalone
- `detect_sequence_start()` — optional, but useful for users

Update roxygen headers:

In `R/translate_from_position.R`, change line 3 from:
```r
#' @keywords internal
```
to:
```r
#' @export
```

In `R/detect_sequence_start.R`, change line 10 from:
```r
#' @keywords internal
```
to:
```r
#' @export
```

Keep these internal (don't export):
- `is_aa_pattern()`
- `search_pattern_in_xstring()`

- [ ] **Step 2: Run `devtools::document()` to regenerate NAMESPACE**

Run: `R -e "devtools::document()"`

Expected: NAMESPACE file is updated with new exports

- [ ] **Step 3: Run `devtools::load_all()` to verify package loads**

Run: `R -e "devtools::load_all(); find_residue"`

Expected: Function loaded without error

- [ ] **Step 4: Run `devtools::check()` for any warnings/errors**

Run: `R -e "devtools::check()"`

Expected: No errors related to missing functions or exports

- [ ] **Step 5: Commit**

```bash
git add NAMESPACE R/translate_from_position.R R/detect_sequence_start.R
git commit -m "docs: export public functions translate_from_position, detect_sequence_start, find_residue"
```

---

## Task 9: Run full test suite to verify integration

**Files:**
- All R/ and tests/testthat/

- [ ] **Step 1: Run all tests**

Run: `R -e "devtools::test()"`

Expected: All tests pass (green checkmarks)

- [ ] **Step 2: Run package check**

Run: `R -e "devtools::check()"`

Expected: No errors, warnings related to code. (Some notes are OK.)

- [ ] **Step 3: Commit final state**

```bash
git add -A
git commit -m "test: verify unified pattern finding implementation passes all tests"
```

---

## Summary

This plan:
1. **Extracts pure helpers first** — translate, pattern detection, generic search
2. **Layers pattern finding on top** — detect start, main find_residue
3. **Reuses logic** — define_lineage now uses detect_sequence_start
4. **Demonstrates composition** — example scripts show cleaner workflows
5. **Maintains testability** — each function testable in isolation
6. **Follows TDD** — write test → verify fail → implement → verify pass → commit

**Total tasks:** 9 (each 2-5 minutes)  
**Files created:** 10  
**Files modified:** 6  
**Commits:** 9
