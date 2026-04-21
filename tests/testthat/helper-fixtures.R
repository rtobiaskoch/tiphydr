# helper-fixtures.R — shared in-memory test objects for all testthat tests
# Loaded automatically by testthat before any test file runs.

# ---------------------------------------------------------------------------
# Private reference sequence (dot-prefix keeps it out of test namespaces)
# ---------------------------------------------------------------------------

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
)  # 80nt — add 4 to reach 84
.ref_seq <- paste0(.ref_seq, "ATCG")  # positions 81-84

# Sanity-check at load time so any accidental edit is caught immediately
stopifnot(nchar(.ref_seq) == 84L)

# ---------------------------------------------------------------------------
# Private helper: substitute bases at positions 10 and 25
# ---------------------------------------------------------------------------

#' Build a variant of .ref_seq by substituting at positions 10 and/or 25.
#'
#' @param base   Character scalar — template sequence (default: .ref_seq).
#' @param pos10  Single character — base to place at position 10.
#' @param pos25  Single character — base to place at position 25.
#' @return Character scalar of the same length as `base`.
.mutate_seq <- function(base = .ref_seq, pos10 = "G", pos25 = "T") {
  s <- base
  substr(s, 10, 10) <- pos10
  substr(s, 25, 25) <- pos25
  s
}

# ---------------------------------------------------------------------------
# Public factory functions — called by individual test files
# ---------------------------------------------------------------------------

#' In-memory DNAStringSet with 5 test sequences.
#'
#' Sequence layout (positions 10 and 25):
#'   NY10_001  : A / C  — both mutations relative to reference
#'   NS2A_002  : A / T  — pos 10 mutation only
#'   NS4B_003  : G / C  — pos 25 mutation only
#'   OTHER_004 : G / T  — no mutations (matches reference)
#'   EXTRA_005 : A / C  — both mutations; intentionally absent from metadata
#'
#' @return Biostrings::DNAStringSet with names in pipe-delimited format
#'   \code{VIRUS|YEAR|STATE|STRAIN_ID}
make_test_biostring <- function() {
  Biostrings::DNAStringSet(c(
    `WNV|2021|CO|NY10_001`    = .mutate_seq(pos10 = "A", pos25 = "C"),  # NY10: both mutations
    `WNV|2020|CO|NS2A_002`    = .mutate_seq(pos10 = "A", pos25 = "T"),  # pos10 only
    `WNV|2022|WY|NS4B_003`    = .mutate_seq(pos10 = "G", pos25 = "C"),  # pos25 only
    `WNV|2019|CO|OTHER_004`   = .mutate_seq(pos10 = "G", pos25 = "T"),  # no mutations
    `WNV|2021|NM|EXTRA_005`   = .mutate_seq(pos10 = "A", pos25 = "C")   # not in metadata
  ))
}

#' In-memory reference DNAStringSet (length 1).
#'
#' @return Biostrings::DNAStringSet, single sequence named "WNV_REF_KU877344"
make_test_ref <- function() {
  Biostrings::DNAStringSet(c(WNV_REF_KU877344 = .ref_seq))
}

#' Metadata tibble — 4 rows, deliberately excludes EXTRA_005.
#'
#' Columns: strain_id (chr), date (Date), location (chr), host (chr).
#' strain_id values match the 4th pipe-field of the DNAStringSet names
#' returned by \code{make_test_biostring()}.
#'
#' @return tibble::tibble with 4 rows
make_test_metadata <- function() {
  tibble::tibble(
    strain_id = c("NY10_001", "NS2A_002", "NS4B_003", "OTHER_004"),
    date      = as.Date(c("2021-07-15", "2020-08-01", "2022-06-20", "2019-05-10")),
    location  = c("Colorado", "Colorado", "Wyoming", "Colorado"),
    host      = c("Culex tarsalis", "Culex pipiens", "Culex tarsalis", "Bird")
  )
}

#' Mutation definitions tibble — 3 rows covering two lineages.
#'
#' NY10     requires mutations at pos 10 (A) AND pos 25 (C) — strict/both.
#' PARTIAL_A requires only pos 10 (A) — single-position lineage.
#'
#' Columns: lineage (chr), pos (int), residue (chr).
#'
#' @return tibble::tibble with 3 rows
make_test_muts <- function() {
  tibble::tibble(
    lineage = c("NY10",  "NY10",  "PARTIAL_A"),
    pos     = c(10L,     25L,     10L),
    residue = c("A",     "C",     "A")
  )
}
