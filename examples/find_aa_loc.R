# Example: Find amino acid locations in DNA sequences using find_residue()
#
# This example shows how to use find_residue() to search for amino acids
# in DNA sequences by auto-translating from a start position.

library(Biostrings)
library(dplyr)

# Create example DNA sequences (starting with ATG)
dna_sequences <- DNAStringSet(c(
  "Seq1" = "ATGAAAAAACTAG",        # ATG (M), AAA (K), AAA (K), CTAstop
  "Seq2" = "ATGAAATTCTAA"          # ATG (M), AAA (K), TTC (F), TAA stop
))

# Find all lysine (K) residues by auto-translating from the detected start codon
# find_residue() automatically:
# 1. Detects that "K" is an amino acid pattern
# 2. Finds the modal ATG start position across all sequences
# 3. Translates DNA to amino acids
# 4. Searches the translated sequences for "K"

result <- find_residue(dna_sequences, pattern = "K", verbose = TRUE)

# Result is a tibble in long format:
# - sequence_name: ID of the sequence
# - location: position in the translated protein
# - pattern: the amino acid searched for (K)
print(result)

# Find multiple amino acids at once
result_multi <- find_residue(
  dna_sequences,
  pattern = c("K", "M"),
  verbose = FALSE
)
print(result_multi)

# Override the start position if you know a different one
result_custom <- find_residue(
  dna_sequences,
  pattern = "K",
  start_position = 1,
  verbose = TRUE
)
print(result_custom)

# Filter results to a specific position (e.g., position 2)
K_position_2 <- result %>%
  filter(pattern == "K", location == 2)
print(K_position_2)
