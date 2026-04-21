#!/usr/bin/env Rscript

# utils/fasta_trim_ref.R
library(Biostrings)
library(ape)
library(ips)

fasta_trim_ref <- function(input_fasta_path, ref_fasta_path, output_fasta_path = "trimmed.fasta") {
  # Read input and reference sequences
  cat("Reading FASTA files...\n")
  ref_fasta <- readDNAStringSet(ref_fasta_path)
  input_fasta <- readDNAStringSet(input_fasta_path)

  cat("Reference fasta length: ", nchar(ref_fasta), "\n")
  
  # Combine sequences
  combined_fasta <- c(ref_fasta, input_fasta)

  cat("Avg length in Combined FASTA: ", mean(nchar(combined_fasta)), "\n")

  # Convert to DNAbin for ips::mafft
  combined_bin <- ape::as.DNAbin(combined_fasta)

  # Align using mafft (in-memory)
  aligned_bin <- ips::mafft(combined_bin)

  # Convert back to DNAStringSet
  aligned_dna <- as.character(aligned_bin)

  cat("Alignment length: ", ncol(aligned_dna), "\n")

  # Extract reference ID and alignment positions
  ref_id <- names(ref_fasta)


  # get the reference row and convert to upper
  ref_row <- toupper(aligned_dna[ref_id, ])
  valid_positions <- which(ref_row %in% c("A", "T", "C", "G"))


  start_pos <- valid_positions[1]
  end_pos <- tail(valid_positions, 1)

  cat("\nTrimming from position ", start_pos, " to ", end_pos, " in the alignment\n")


  # Subset the matrix to only keep valid columns (non-gap positions in ref)
  aligned_trimmed <- aligned_dna[, valid_positions, drop = FALSE]


  # Remove the 'ref' row
  aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) != "ref", , drop = FALSE]

  

  cat("Trimmed alignment length: ", ncol(aligned_trimmed), "\n")


  if (ncol(aligned_trimmed) != nchar(ref_fasta)) {
    warning("Trimmed alignment is not the same length as the reference\n")
  }


  # Collapse matrix rows into sequences
  seqs <- apply(aligned_trimmed, 1, paste0, collapse = "")


  # Convert to DNAStringSet
  trimmed_set <- DNAStringSet(seqs)
  names(trimmed_set) <- rownames(aligned_trimmed)

  cat("Saving trimmed FASTA to: ", output_fasta_path, "\n")
  writeXStringSet(trimmed_set, filepath = output_fasta_path)


}

# Example usage:
# fasta_trim_ref("input.fasta", "reference.fasta", "output.fasta")