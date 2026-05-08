gff_path <- test_path("testdata/test_genome.gff")

test_that("returns tibble with correct columns", {
  result <- resolve_mut_positions(make_test_gene_muts(), gff_path)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("lineage", "gene", "pos", "residue", "comparison_type", "codon_start", "aa_start"))
  expect_equal(nrow(result), 4L)
})

test_that("AA rows get correct codon_start from GFF", {
  result <- resolve_mut_positions(make_test_gene_muts(), gff_path)
  # NY10/gene_A AA2 → 1  + (2-1)*3 = 4
  # NY10/gene_B AA1 → 10 + (1-1)*3 = 10
  # NS4A/gene_A AA2 → 1  + (2-1)*3 = 4
  # NS4B/gene_B AA1 → 10 + (1-1)*3 = 10
  expect_equal(result$comparison_type, c("aa", "aa", "aa", "aa"))
  expect_equal(result$codon_start,     c(4L, 10L, 4L, 10L))
})

test_that("nuc rows pass through with NA codon_start and unchanged pos", {
  nuc_muts <- tibble::tibble(lineage = "L1", gene = "nuc", pos = 42L, residue = "A")
  result   <- resolve_mut_positions(nuc_muts, gff_path)
  expect_equal(result$comparison_type, "nuc")
  expect_true(is.na(result$codon_start))
  expect_equal(result$pos, 42L)
})

test_that("mixed gene + nuc table resolves each row independently", {
  mixed <- tibble::tibble(
    lineage = c("NY10",    "NY10"),
    gene    = c("gene_A",  "nuc"),
    pos     = c(2L,        10L),
    residue = c("E",       "A")
  )
  result <- resolve_mut_positions(mixed, gff_path)
  expect_equal(result$comparison_type, c("aa",  "nuc"))
  # gene_A start=1, AA2 → codon_start = 4
  expect_equal(result$codon_start,     c(4L,   NA_integer_))
})

test_that("unrecognised gene name errors with gene name in message", {
  bad <- tibble::tibble(lineage = "L", gene = "nonexistent_protein", pos = 1L, residue = "A")
  expect_error(resolve_mut_positions(bad, gff_path), "nonexistent_protein")
})

test_that("ambiguous gene name errors listing the candidates", {
  # Write a temp GFF with two non-slash products that both contain "NS"
  tmp <- tempfile(fileext = ".gff")
  writeLines(c(
    "##gff-version 3",
    paste(c("s", ".", "mature_protein_region_of_CDS", "1",  "300", ".", "+", ".", "product=nonstructural protein A"), collapse = "\t"),
    paste(c("s", ".", "mature_protein_region_of_CDS", "301","600", ".", "+", ".", "product=nonstructural protein B"), collapse = "\t")
  ), tmp)
  ambig <- tibble::tibble(lineage = "L", gene = "NS", pos = 1L, residue = "A")
  expect_error(resolve_mut_positions(ambig, tmp), "matches multiple products")
  unlink(tmp)
})

test_that("slash products are deprioritised — NS4B resolves to canonical protein", {
  # Write a temp GFF mirroring the real WNV case: one slash product, one canonical
  tmp <- tempfile(fileext = ".gff")
  writeLines(c(
    "##gff-version 3",
    paste(c("s", ".", "mature_protein_region_of_CDS", "6916", "7311", ".", "+", ".", "product=N-NS4B/WARF4"),            collapse = "\t"),
    paste(c("s", ".", "mature_protein_region_of_CDS", "6916", "7680", ".", "+", ".", "product=nonstructural protein NS4B"), collapse = "\t")
  ), tmp)
  muts <- tibble::tibble(lineage = "L", gene = "NS4B", pos = 1L, residue = "A")
  result <- resolve_mut_positions(muts, tmp)
  # Should resolve to nonstructural protein NS4B (nt_start 6916), codon_start = 6916
  expect_equal(result$comparison_type, "aa")
  expect_equal(result$codon_start, 6916L)
  unlink(tmp)
})

test_that("aa_start is NA when cds_start is not supplied", {
  result <- resolve_mut_positions(make_test_gene_muts(), gff_path)
  expect_true(all(is.na(result$aa_start)))
})

test_that("aa_start computed correctly when cds_start is supplied", {
  # cds_start = 1 (test genome CDS starts at nt 1)
  # gene_A start=1, AA2: codon_start=4  → ceiling((4-1)/3+1)  = ceiling(2)  = 2
  # gene_B start=10, AA1: codon_start=10 → ceiling((10-1)/3+1) = ceiling(4)  = 4
  result <- resolve_mut_positions(make_test_gene_muts(), gff_path, cds_start = 1L)
  expect_equal(result$aa_start, c(2L, 4L, 2L, 4L))
})

test_that("nuc rows always have NA aa_start even when cds_start is supplied", {
  nuc_muts <- tibble::tibble(lineage = "L1", gene = "nuc", pos = 42L, residue = "A")
  result   <- resolve_mut_positions(nuc_muts, gff_path, cds_start = 97L)
  expect_true(is.na(result$aa_start))
})

test_that("cds_start validation rejects non-scalar or non-numeric input", {
  expect_error(resolve_mut_positions(make_test_gene_muts(), gff_path, cds_start = c(1L, 2L)), "single integer")
  expect_error(resolve_mut_positions(make_test_gene_muts(), gff_path, cds_start = NA_integer_), "single integer")
})

test_that("missing required column errors informatively", {
  bad <- tibble::tibble(lineage = "L", gene = "env", residue = "V")  # no pos
  expect_error(resolve_mut_positions(bad, gff_path), "pos")
})

test_that("nonexistent GFF path errors informatively", {
  expect_error(
    resolve_mut_positions(make_test_gene_muts(), "does_not_exist.gff"),
    "not found"
  )
})
