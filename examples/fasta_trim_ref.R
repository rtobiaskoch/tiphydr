#!/usr/bin/env Rscript

# utils/fasta_trim_ref.R
library(Biostrings)
library(ape)
library(ips)

fasta_trim_ref <- function(input_fasta, ref_fasta) {
  
  # Require input fasta and ref fasta to be DNAstringsets
  if (class(input_fasta) != "DNAStringSet") {
    stop("Input FASTA must be a DNAStringSet")
  }
  if (class(ref_fasta) != "DNAStringSet") {
    stop("Reference FASTA must be a DNAStringSet")}

  # Read input and reference sequences
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

  cat("\nTrimming from position ", start_pos, " to ", end_pos, " in the alignment\n",
      "Total trim length should be ", end_pos-start_pos+1, "\n")


  # Subset the matrix to only keep valid columns (non-gap positions in ref)
  aligned_trimmed <- aligned_dna[, valid_positions, drop = FALSE]


  # Remove the 'ref' row
  aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) != names(ref_fasta), ,drop = FALSE]


  if (ncol(aligned_trimmed) != nchar(ref_fasta)) {
    warning("Trimmed alignment is not the same length as the reference\n")
  }

  # Collapse matrix rows into sequences
  seqs <- apply(aligned_trimmed, 1, paste0, collapse = "")

  # Convert to DNAStringSet
  trimmed_set <- DNAStringSet(seqs)
  names(trimmed_set) <- rownames(aligned_trimmed)

  return(trimmed_set)

}


# Example usage:
# fasta_trim_ref("input.fasta", "reference.fasta", "output.fasta")
# fasta_trim_ref(DNAStringSet(c("ATCG", "ATGC")), DNAStringSet(c("ATCG", "TGCA")))
