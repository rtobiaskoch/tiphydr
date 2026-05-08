#test
x = Biostrings::readDNAStringSet("tests/testthat/testdata/wnv_seqs.fasta")
ref = Biostrings::readDNAStringSet("tests/testthat/testdata/wnv_ref.fasta")
x = fasta_trim_ref(x, ref)

unname(as.character(x))
# 2. phyDat expects a named character vector — as.character() provides that
m <- ape::dist.dna(ape::as.DNAbin(x), model = "K80")
tree <- ape::nj(m)
plot(tree)
