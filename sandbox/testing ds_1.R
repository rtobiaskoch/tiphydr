library(tidyverse)
library(Biostrings)
devtools::load_all(".")

muts <- wnv_mut_orf
muts_gff <- wnv_mut_gff
ref <- fasta_read("data/wnv_ref.fasta")
trimmed <- fasta_read('data/wnv_sequences.fasta')
trimmed <- Biostrings::chartr(".-", "NN", trimmed)


wn02 <- find_residue(trimmed, pattern = "A") %>%
  dplyr::filter(location == 449)
#
# wn02 %>% dplyr::group_by(location) %>% dplyr::count()

# Resolve gene-relative AA mutations to alignment AA coordinates, then feed
# aa_start in as `pos` for the mut_type = "aa" pathway (gff is no longer a
# define_lineage argument).
muts_gff_resolved <- resolve_mut_positions(
  muts_gff,
  gff = "data/genomic.gff",
  cds_start = 97L
) %>%
  dplyr::mutate(pos = aa_start)

lineages <- define_lineage(
  trimmed[1:10],
  ref,
  muts_gff_resolved,
  mut_type = "aa"
)

lineages <- define_lineage(trimmed, ref, muts, mut_type = "aa")

lin_sum <- lineages %>%
  dplyr::count(lineage)

write.csv(lin_sum, "sandbox/define_lineage_ds_1_output.csv", row.names = FALSE)

unknown <- c(
  "VCTR0010000464|2021-08-25|West|AZ/Maricopa|Unknown|mosquito-culex",
  "PQ005956|2022-06-01|West|CA/LosAngeles|Unknown|mosquito-culex",
  "OZ201316|2018-07-18|Northeast|CT/Unknown|Unknown|mosquito-culex",
  "OZ201866|2018-07-12|Northeast|CT/Unknown|Unknown|mosquito-culex"
)

unknown <- purrr::keep(trimmed, names(trimmed) %in% unknown)
unknown <- Biostrings::translate(unknown, if.fuzzy.codon = "X")
as.character(subseq(unknown, start = 449, end = 449))
as.character(subseq(unknown, start = 2513, end = 2513))
