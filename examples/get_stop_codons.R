#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# --------------------------------  H E A D E R -------------------------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load required libraries

library(Biostrings)
library(dplyr)
library(stringr)


#pulls get_aa_mut function

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#  R E A D 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fasta0 = readDNAStringSet("results/alignment_downsampled_trimmed.fasta")

dna_sequences <- DNAStringSet(c(
  "Seq1" = "AATGCCGTAGTTGAATAA",  # ATG in frame 2 (pos 2), TGA in frame 1 (pos 11)
  "Seq2" = "CATGGGATGCAA",         # Two ATGs (frame 1 pos 2, frame 3 pos 8)
  "Seq3" = "TTTATTAA"          # ATG in frame 3 (pos 4)
))

# Search for all stop codons
stop_codons <- c("TAA", "TAG", "TGA")

test = find_nucleotide_patterns(dna_sequences, stop_codons, frame = 1)
results <- find_nucleotide_patterns(fasta0, stop_codons, frame = 1)


count_stop = results %>%
  group_by(sequence_id) %>%
  summarise(count = sum(count))
