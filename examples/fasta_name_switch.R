fasta_name_switch <- function(fasta, match_col, replace_col) {
  # Get current names from fasta object
  fasta_names <- names(fasta)
  
  # Find matches between fasta names and accessions
  matches <- match(fasta_names, match_col)
  
  # Replace matched names with corresponding strains
  # Keep original names where there is no match (where matches is NA)
  names(fasta)[!is.na(matches)] <-  replace_col[matches[!is.na(matches)]]
  
  return(fasta)
}