# Example: Find nucleotide patterns using find_residue()
#
# This example shows how to use find_residue() to search for DNA patterns
# (like stop codons or specific motifs) in DNA sequences.

library(Biostrings)
library(dplyr)

# Create example DNA sequences
dna_sequences <- DNAStringSet(c(
  "Seq1" = "AATGCCGTAGTTGAATAA",  # Contains TGA at position 11
  "Seq2" = "CATGGGATGCAA",        # Contains TAG at position 10
  "Seq3" = "TTTATTAA"             # Contains TAA at position 6
))

# Find all stop codons (TAA, TAG, TGA) using find_residue()
stop_codons <- c("TAA", "TAG", "TGA")

# find_residue() recognizes these as nucleotide patterns (not AA)
# and searches the DNA sequences directly
result <- find_residue(dna_sequences, pattern = stop_codons, verbose = FALSE)

# Result is in long format: one row per match per sequence per pattern
print(result)

# Count occurrences by sequence
stop_codon_summary <- result %>%
  group_by(sequence_name) %>%
  summarise(count = n(), .groups = "drop")

print(stop_codon_summary)

# Filter for specific stop codons
taa_codons <- result %>%
  filter(pattern == "TAA")

print(taa_codons)

# Find start codons
atg_result <- find_residue(dna_sequences, pattern = "ATG", verbose = FALSE)
print(atg_result)
