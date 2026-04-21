#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# --------------------------------  H E A D E R -------------------------------------------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Load required libraries

library(Biostrings)
library(dplyr)
library(stringr)


rm(list = ls())
source(".Rprofile")
#pulls get_aa_mut function

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#  R E A D 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fasta0 = readDNAStringSet("results/alignment_downsampled_trimmed.fasta")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#--------------I D   S T A R T I N G   I N D E X --------------
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

start_codons = get_start_met(fasta0)

start_index_mode = start_codons %>%
  group_by(start_codon_position) %>%
  count %>%
  ungroup %>%
  filter(n == max(n))

start_index_ref = start_codons %>%
  dplyr::filter(str_detect(tolower(sequence_id), ref))

if(start_index_ref$start_codon_position != start_index_mode$start_codon_position) {
  cat("\n Warning: start codon position mode doesn't match the reference")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# P U L L 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#get sample list true false for a mutation
ns2a_188 = get_aa_mut(fasta0, 1331, "K", start_index = 97)
ns4b_240 = get_aa_mut(fasta0, 2513, "M", start_index = 97)

ny10 = left_join(ns2a_188, ns4b_240, by = "sequence_name") %>%
  mutate(ny10 = case_when(mutation_found.x & mutation_found.y ~ "NY10",
                          mutation_found.x ~ "NS2A_188K",
                          mutation_found.y ~ "NS4B_240M",
                          T ~ "other"))
table(ny10$ny10)      
         

ny10_out = ny10 %>%
  rename(strain = sequence_name) %>%
 dplyr::select(strain, ny10)
         

write.csv(ny10_out, "data_output/ny10.csv")
write.csv(ny10_out, "../3_persistence_MODELING/data/ny10.csv")
