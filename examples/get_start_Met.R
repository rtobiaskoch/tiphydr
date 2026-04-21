get_start_met <- function(dna_string_set) {
  # Initialize vectors to store results
  sequence_ids <- c()
  start_positions <- c()
  frames <- c()
  
  # Loop through each sequence in the DNAStringSet
  for (i in seq_along(dna_string_set)) {
    seq <- dna_string_set[[i]]
    seq_name <- names(dna_string_set)[i]
    
    # Find the position of the first start codon (ATG)
    start_codon <- "ATG"
    start_pos <- matchPattern(start_codon, seq)
    
    # If a start codon is found, record its position and frame
    if (length(start_pos) > 0) {
      start_index <- start(start_pos)[1]
      frame <- (start_index - 1) %% 3  # Calculate the frame (0, 1, or 2)
      sequence_ids <- c(sequence_ids, seq_name)
      start_positions <- c(start_positions, start_index)
      frames <- c(frames, frame)
    } else {
      # If no start codon is found, record NA
      sequence_ids <- c(sequence_ids, seq_name)
      start_positions <- c(start_positions, NA)
      frames <- c(frames, NA)
    }
  }
  
  # Create a data frame from the results
  results_df <- data.frame(
    sequence_id = sequence_ids,
    start_codon_position = start_positions,
    frame = frames,
    stringsAsFactors = FALSE
  )
  
  # Return the data frame
  return(results_df)
}

# Example usage
#results_df <- get_start_met(dna_sequences)
#print(results_df)