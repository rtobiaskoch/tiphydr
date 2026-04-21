# R script to translate sequences and identify amino acid location

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}

rm(list = ls())
source(".Rprofile") #loads get_start_met and fina_aa_loc
# Load the Biostrings package
library(Biostrings)
library(tidyverse)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#USER INPUT
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
data <- readDNAStringSet('results/.fasta')
ref = "nc_009942"
K_loc_filter = 1331
M_loc_filter

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ID STARTING INDEX FOR TRANSLATION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

start_codons = get_start_met(data)

start_index_mode = start_codons %>%
  group_by(start_codon_position) %>%
  count %>%
  ungroup %>%
  filter(n == max(n))

start_index_ref = start_codons %>%
  filter(str_detect(tolower(sequence_id), ref))

if(start_index_ref$start_codon_position != start_index_mode$start_codon_position) {
  cat("\n Warning: start codon position mode doesn't match the reference")
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#CALL FUNCTION AND ID 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Call the function
K <- find_aa_loc(data, aa = "K", 
                 start_codon = F, 
                 start_index = 97)

K_filter = K %>%
  filter(location == K_loc_filter)

K_count = K %>%
  group_by(amino_acid, location) %>%
  count() %>%
  filter(n > 10 & n <4000)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#CALL FUNCTION AND ID 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Call the function
M <- find_aa_loc(data, aa = "M", 
                 start_codon = F, 
                 start_index = 97)

  
