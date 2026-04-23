# Example: Detect the modal start codon position using detect_sequence_start()
#
# This example shows how to use detect_sequence_start() to find the most
# common (modal) ATG start codon position across a set of DNA sequences.

library(Biostrings)
library(dplyr)

# Create example DNA sequences with ATG start codons
dna_sequences <- DNAStringSet(c(
  "Seq1" = "ATGAAAAAACTAG",      # ATG at position 1
  "Seq2" = "ATGAAATTCTAA",       # ATG at position 1
  "Seq3" = "GGGATGCCCTAG",       # ATG at position 4
  "Seq4" = "ATGCCCCAA"           # ATG at position 1
))

# Find the modal (most common) start codon position
# Returns a list with:
#   - position: the most frequently occurring position
#   - table: frequency table of all positions found
start_info <- detect_sequence_start(
  dna_sequences,
  pattern = "ATG",
  verbose = TRUE
)

# Access the modal position
modal_position <- start_info$position
cat("Modal ATG position:", modal_position, "\n")

# View the frequency table
cat("\nFrequency table of start codon positions:\n")
print(start_info$table)

# Common workflow: use this to establish a consensus translation start
# then translate all sequences from that position
aa_sequences <- translate_from_position(
  dna_sequences,
  start_position = modal_position
)
print(aa_sequences)
