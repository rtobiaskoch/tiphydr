# Example: Find stop codons using find_residue()
#
# This example shows how to use find_residue() to locate all stop codons
# (TAA, TAG, TGA) in a set of DNA sequences.

library(Biostrings)
library(dplyr)

# Create example DNA sequences
dna_sequences <- DNAStringSet(c(
  "Seq1" = "AATGCCGTAGTTGAATAA",  # TAA at position 17, TGA at position 11
  "Seq2" = "CATGGGATGCAA",        # TAG at position 10
  "Seq3" = "TTTATTAA"             # TAA at position 6
))

# Define stop codon patterns
stop_codons <- c("TAA", "TAG", "TGA")

# Find all stop codons in the sequences
result <- find_residue(dna_sequences, pattern = stop_codons, verbose = FALSE)

# Result is a tibble with:
# - sequence_name: ID of the sequence
# - location: position in the DNA sequence
# - pattern: which stop codon (TAA, TAG, or TGA)
print(result)

# Summary: count stop codons per sequence
stop_codon_count <- result %>%
  group_by(sequence_name) %>%
  summarise(
    total_stops = n(),
    .groups = "drop"
  )

print(stop_codon_count)

# Advanced: count by stop codon type
stop_by_type <- result %>%
  group_by(sequence_name, pattern) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = pattern,
    values_from = count,
    values_fill = 0
  )

print(stop_by_type)
