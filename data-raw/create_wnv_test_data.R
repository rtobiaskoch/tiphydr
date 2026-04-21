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
