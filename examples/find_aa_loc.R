# Function to count occurrences and locations of an amino acid in a DNAStringSet
find_aa_loc <- function(dna_string_set, aa, start_codon = TRUE, start_index = NULL) {
  # Initialize vectors to store results
  sequence_ids <- c()
  amino_acids <- c()
  locations <- c()
  
  # Loop through each sequence in the DNAStringSet
  for (i in seq_along(dna_string_set)) {
    seq <- dna_string_set[[i]]
    seq_name <- names(dna_string_set)[i]
    
    # Determine the starting position for translation
    if (start_codon) {
      # Find the position of the first start codon (ATG)
      start_codon_seq <- "ATG"
      start_pos <- matchPattern(start_codon_seq, seq)
      
      # If a start codon is found, use its position
      if (length(start_pos) > 0) {
        start_index_seq <- start(start_pos)[1]  # Get the start position of the first match
      } else {
        # If no start codon is found, skip this sequence
        next
      }
    } else if (!is.null(start_index)) {
      # Use the provided start_index
      start_index_seq <- start_index
    } else {
      # If start_codon is FALSE and start_index is NULL, start from the beginning
      start_index_seq <- 1
    }
    
    # Extract the sequence from the starting position
    cds_seq <- subseq(seq, start = start_index_seq)
    
    # Translate the sequence to amino acids
    aa_seq <- translate(cds_seq, if.fuzzy.codon = "X")
    
    # Convert the AA sequence to a character vector
    aa_seq_char <- as.character(aa_seq)
    
    # Find indices where the amino acid occurs
    loc <- which(unlist(strsplit(aa_seq_char, "")) == aa)
    
    # If the amino acid is found, add each location as a separate row
    if (length(loc) > 0) {
      sequence_ids <- c(sequence_ids, rep(seq_name, length(loc)))
      amino_acids <- c(amino_acids, rep(aa, length(loc)))
      locations <- c(locations, loc)
    }
  }
  
  # Create a data frame from the results
  results_df <- data.frame(
    sequence_id = sequence_ids,
    amino_acid = amino_acids,
    location = locations,
    stringsAsFactors = FALSE
  )
  
  # Return the data frame
  return(results_df)
}