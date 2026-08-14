library(Biostrings)
library(tiphydr)
library(dplyr)

mdata_cols <- c("accession", "date", "deme", "county", "lineage", "host")
ref <- Biostrings::readDNAStringSet("data/wnv_ref.fasta")
muts <- tiphydr::wnv_mut_orf

seq <- Biostrings::readDNAStringSet("sandbox/beast_ds_0.25_US.fasta")


lin <- tiphydr::define_lineage(seq, ref = ref, muts = muts, mut_type = "aa")

mdata <- tiphydr::extract_metadata(
  seq,
  mdata_cols
) %>%
  mutate(header = names(seq)) %>%
  mutate(lineage = lin$lineage) %>%
  tidyr::unite("strain", all_of(mdata_cols), sep = "|", remove = FALSE)

seq <- tiphydr::fasta_name_switch(
  seq,
  match_col = mdata$accession,
  replace_col = mdata$strain
)
names(seq)
