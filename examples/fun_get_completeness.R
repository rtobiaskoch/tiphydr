require(Biostrings)
require(tibble)
require(dplyr)

genome_cmplt = function(fasta, threshold = 0.9) {
  
  nucs = c("A", "T", "C", "G")
  
  cmplt = as_tibble(Biostrings::alphabetFrequency(fasta)) %>%
    rowwise() %>%
    mutate(ambig = sum(c_across(-nucs), na.rm = TRUE)) %>%
    dplyr::select(all_of(nucs), ambig) %>%
    mutate(length = sum(c_across(everything()), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(pct_ambig = round(ambig/length,4)) %>%
    mutate(pct_cmplt = round((length - ambig)/max(length),4)) %>%
    mutate(strain = names(fasta)) %>%
    dplyr::filter(pct_cmplt > threshold)
  
  cat("Starting number of sequences: ", length(fasta), "\n\n")
  cat(length(fasta) - nrow(cmplt), "Sequences removed for being less than", threshold, "coverage")
  
  return(cmplt)
}