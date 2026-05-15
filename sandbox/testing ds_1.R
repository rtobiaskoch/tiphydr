library(tidyverse)
devtools::load_all(".")

muts = read.csv("data/lineage_mutations_ORF.csv")
muts_gff = read.csv("data/lineage_mutations_gff.csv")
ref <- fasta_read("data/wnv_ref.fasta")
trimmed <- fasta_read('data/wnv_sequences.fasta')
trimmed = Biostrings::chartr(".-", "NN", trimmed)


wn02 = find_residue(trimmed, pattern = "A") %>%
  dplyr::filter(location == 449)
#
# wn02 %>% dplyr::group_by(location) %>% dplyr::count()

resolve_mut_positions(muts_gff, gff = "data/genomic.gff", cds_start = 97L)

lineages = define_lineage(
  trimmed[1:10],
  ref,
  muts_gff,
  input_type = "aa",
  gff = "data/genomic.gff"
)

lineages = define_lineage(trimmed, ref, muts, input_type = "aa")

lineages %>%
  dplyr::group_by(lineage) %>%
  dplyr::count()


muts


unknown = c(
  "VCTR0010000464|2021-08-25|West|AZ/Maricopa|Unknown|mosquito-culex",
  "PQ005956|2022-06-01|West|CA/LosAngeles|Unknown|mosquito-culex",
  "OZ201316|2018-07-18|Northeast|CT/Unknown|Unknown|mosquito-culex",
  "OZ201866|2018-07-12|Northeast|CT/Unknown|Unknown|mosquito-culex"
)

unknown = purrr::keep(trimmed, names(trimmed) %in% unknown)
unknown = translate(unknown, if.fuzzy.codon = "X")
as.character(subseq(unknown, start = 449, end = 449))
as.character(subseq(unknown, start = 2513, end = 2513))


# Extract position 85 for every sequence in the AAStringSet

t = subseq(trimmed, start = 4087 + 1, end = 4089 + 1)
#t = translate(t, if.fuzzy.codon = "X")

Biostrings::consensusMatrix(t)
Biostrings::codons(t)

# This returns a character vector like: c("A", "A", "T", "A"...)
