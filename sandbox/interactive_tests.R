devtools::load_all(".")

f = make_test_biostring()
mdata = extract_metadata(biostring = f, names= c("id","date", "state", "lineage"), delim = "|")
r = make_test_ref()
m = make_test_muts()
ma = make_test_aa_muts()
mg = make_test_gene_muts()
gff = readLines("tests/testthat/testdata/test_genome.gff")

result <- define_lineage(
  f, r, m,
  verbose = TRUE
)
result


result2 <- define_lineage(
  f, r, ma, type = "aa",
  verbose = TRUE
)

result2
