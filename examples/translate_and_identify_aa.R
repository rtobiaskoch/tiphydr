# Example: Translate DNA and identify amino acid locations
#
# This example shows a complete workflow using detect_sequence_start()
# and find_residue() to identify amino acid positions in DNA sequences.

library(Biostrings)
library(dplyr)

# Create example DNA sequences (starting with ATG)
dna_sequences <- DNAStringSet(c(
  "Seq1" = "ATGAAAAAATTTTAG",       # M K K F stop
  "Seq2" = "ATGAAAAAATAG",          # M K K stop
  "Seq3" = "ATGAAAAAAAAATTAG"       # M K K K F stop
))

# ========================================================================
# Step 1: Detect the consensus start codon position
# ========================================================================
# Find the modal ATG position across all sequences
start_info <- detect_sequence_start(
  dna_sequences,
  pattern = "ATG",
  verbose = TRUE
)

modal_start <- start_info$position

# ========================================================================
# Step 2: Find amino acid locations
# ========================================================================
# Method 1: Auto-detect start position (uses ATG by default)
K_locations_auto <- find_residue(
  dna_sequences,
  pattern = "K",
  verbose = TRUE
)

print("Lysine (K) locations (auto-start):")
print(K_locations_auto)

# Method 2: Specify a custom start position
K_locations_custom <- find_residue(
  dna_sequences,
  pattern = "K",
  start_position = 1,
  verbose = FALSE
)

print("Lysine (K) locations (position 1):")
print(K_locations_custom)

# ========================================================================
# Step 3: Find multiple amino acids
# ========================================================================
multi_aa <- find_residue(
  dna_sequences,
  pattern = c("K", "M", "F"),
  verbose = FALSE
)

print("Multiple amino acids:")
print(multi_aa)

# ========================================================================
# Step 4: Analysis and filtering
# ========================================================================
# Count amino acids per sequence
aa_counts <- multi_aa %>%
  group_by(sequence_name, pattern) %>%
  summarise(count = n(), .groups = "drop")

print("Amino acid counts per sequence:")
print(aa_counts)

# Find sequences with K at position 2
K_at_pos_2 <- K_locations_auto %>%
  filter(location == 2)

print("Sequences with K at position 2:")
print(K_at_pos_2)

# ========================================================================
# Step 5: Get translated sequences for validation
# ========================================================================
aa_sequences <- translate_from_position(
  dna_sequences,
  start_position = modal_start
)

print("Translated sequences:")
print(aa_sequences)
