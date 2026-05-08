devtools::load_all(".")
library(dplyr)

#test
x = Biostrings::readDNAStringSet("tests/testthat/testdata/wnv_seqs.fasta")
ref = Biostrings::readDNAStringSet("tests/testthat/testdata/wnv_ref.fasta")
x = fasta_trim_ref(x, ref)

#make dummy tree
m <- ape::dist.dna(ape::as.DNAbin(x), model = "K80")
tree <- ape::nj(m)
plot(tree)

#extract metadata to make a new test columm
mdata = extract_metadata(
  x,
  c("id", "year", "state", "lineage"),
  delim = "|"
) %>%
  mutate(newname = 1:nrow(.))


newtree = tree_tip_swap(tree, mdata$lineage, mdata$newname)

plot(newtree)
