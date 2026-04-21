find_nucleotide_patterns <- function(dna_string_set, 
                                     patterns = "ATG", 
                                     frame = 1) {
  # Validate frame argument
  if (!frame %in% 1:3) {
    stop("Frame must be 1, 2, or 3")
  }
  
  # Initialize results storage
  results_list <- list()
  
  # Loop through each pattern
  for (pattern in patterns) {
    pattern_results <- data.frame(
      sequence_id = names(dna_string_set),
      pattern = pattern,
      count = 0,
      positions = I(vector("list", length(dna_string_set))),
      frame = frame,
      stringsAsFactors = FALSE
    )
    
    # Find matches for the current pattern
    for (i in seq_along(dna_string_set)) {
      seq <- dna_string_set[[i]]
      matches <- matchPattern(pattern, seq)
      
      if (length(matches) > 0) {
        # Correct frame calculation: (position - frame) %% 3 == 0
        match_starts <- start(matches)
        in_frame <- (match_starts - frame) %% 3 == 0
        frame_matches <- matches[in_frame]
        
        # Store results
        pattern_results$count[i] <- length(frame_matches)
        pattern_results$positions[[i]] <- if (length(frame_matches) > 0) {
          start(frame_matches)
        } else {
          NA
        }
      } else {
        pattern_results$positions[[i]] <- NA
      }
    }
    
    results_list[[pattern]] <- pattern_results
  }
  
  # Combine all results
  combined_results <- do.call(rbind, results_list)
  rownames(combined_results) <- NULL
  
  return(combined_results)
}