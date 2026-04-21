
# Optimized function to check for specific amino acid mutations in DNA sequences
get_aa_mut <- function(fasta, position, amino_acid, start_index = 1) {
  # Load required library
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package is required. Please install it with BiocManager::install('Biostrings')")
  }
  
  # Pre-allocate vectors for results
  sequence_names <- names(fasta)
  n_sequences <- length(fasta)
  actual_aa_vector <- character(n_sequences)
  mutation_found_vector <- logical(n_sequences)
  
  # Translate all sequences at once
  trim_fasta = subseq(fasta, start = start_index)
  aa_seqs <- Biostrings::translate(trim_fasta, if.fuzzy.codon = "X")
  
  # Vectorized extraction of amino acids at the specified position
  # Create a position filter - only check sequences long enough
  valid_positions <- width(aa_seqs) >= position
  
  # Initialize all as FALSE/NA
  mutation_found_vector[] <- FALSE
  actual_aa_vector[] <- NA_character_
  
  if (any(valid_positions)) {
    # Extract the amino acid at the specified position for all valid sequences
    # Using Biostrings::subseq is much faster than iterating
    actual_aas <- as.character(Biostrings::subseq(aa_seqs[valid_positions], 
                                                  start = position, 
                                                  end = position))
    
    # Store the actual amino acids
    actual_aa_vector[valid_positions] <- actual_aas
    
    # Check which ones match the target amino acid
    mutation_found_vector[valid_positions] <- actual_aas == amino_acid
    
    # Mark sequences that are too short
    actual_aa_vector[!valid_positions] <- "Position out of bounds"
  }
  
  # Create the results data frame
  results_df <- data.frame(
    sequence_name = sequence_names,
    position = position,
    target_amino_acid = amino_acid,
    actual_amino_acid = actual_aa_vector,
    mutation_found = mutation_found_vector,
    stringsAsFactors = FALSE) %>%
    mutate(mutation_found = if_else(actual_amino_acid == "X", 
                                    "unknown", 
                                    as.character(mutation_found))
           )
  
  return(results_df)
}

# Example usage:
# library(Biostrings)
# fasta0 <- readDNAStringSet("your_sequences.fasta")
# results <- check_amino_acid_mutation(fasta0, 1331, "K")
