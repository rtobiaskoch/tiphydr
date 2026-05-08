test_that("returns dataframe with strain and lineage columns", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_s3_class(result, "data.frame")
  expect_named(result, c("strain", "lineage"))
  expect_equal(nrow(result), length(make_test_biostring()))
})

test_that("strict match: NY10 requires both mutations", {
  # NY10_001 has pos4=G and pos10=C → matches NY10 (more specific than NS4A or NS4B)
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("single-mutation sequences assign their single-gene lineage", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  # NS2A_002: pos4=G only → NS4A
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "NS4A")
  # NS4B_003: pos10=C only → NS4B
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "NS4B")
})

test_that("wildtype sequence assigns unknown", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|OTHER_004"], "unknown")
})

test_that("multiple-match: assigns most-specific (most mutations)", {
  # NY10_001 matches NY10 (2 muts), NS4A (1 mut), and NS4B (1 mut) → NY10 wins
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("muts missing required columns throws error", {
  bad_muts <- tibble::tibble(lineage = "NY10", position = 10L, base = "A")
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), bad_muts),
    "lineage.*pos.*residue"
  )
})

test_that("invalid type throws error", {
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts(),
                   type = "protein"),
    "type must be"
  )
})

# ── gff pathway tests ──────────────────────────────────────────────────────────

test_that("gff pathway assigns same lineages as nuc pathway (regression)", {
  gff_path <- test_path("testdata/test_genome.gff")

  result_gff <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_gene_muts(), gff = gff_path)
  )
  result_nuc <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )

  # Align by strain name so order does not matter
  gff_lineages <- result_gff$lineage[order(result_gff$strain)]
  nuc_lineages <- result_nuc$lineage[order(result_nuc$strain)]
  expect_equal(gff_lineages, nuc_lineages)
})

test_that("gff pathway: strict match requires all AA mutations", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_gene_muts(), gff = test_path("testdata/test_genome.gff"))
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("gff pathway: single-mutation sequences assign their single-gene lineage", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_gene_muts(), gff = test_path("testdata/test_genome.gff"))
  )
  # NS2A_002: gene_A AA2=E only → NS4A
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "NS4A")
  # NS4B_003: gene_B AA1=Q only → NS4B
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "NS4B")
})

test_that("gff pathway: wildtype assigns unknown", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_gene_muts(), gff = test_path("testdata/test_genome.gff"))
  )
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|OTHER_004"], "unknown")
})

test_that("mixed gene + nuc muts table assigns lineages correctly", {
  # One AA row (gene_A AA2 = Glu) + one nuc row (pos10 = C) — both must match for NY10
  mixed_muts <- tibble::tibble(
    lineage = c("NY10",    "NY10"),
    gene    = c("gene_A",  "nuc"),
    pos     = c(2L,        10L),   # gene_A AA2=Glu; nuc pos10=C required for NY10
    residue = c("E",       "C")
  )
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   mixed_muts, gff = test_path("testdata/test_genome.gff"))
  )
  # NY10_001 has gene_A AA2=Glu and nuc10=C → NY10
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
  # NS2A_002 has gene_A AA2=Glu but nuc10=A → unknown
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "unknown")
})

test_that("gene column without gff and non-nuc genes errors informatively", {
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_gene_muts()),
    "gff"
  )
})

test_that("gene column all-nuc without gff falls through to nuc path", {
  # Inline lineage name "PARTIAL_A" is arbitrary here — tests the pathway, not the lineage names
  nuc_gene_muts <- tibble::tibble(
    lineage = c("NY10",  "NY10",  "PARTIAL_A"),
    gene    = c("nuc",   "nuc",   "nuc"),
    pos     = c(4L,      10L,     4L),
    residue = c("G",     "C",     "G")
  )
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), nuc_gene_muts)
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "PARTIAL_A")
})

# ── type = "aa" pathway tests ─────────────────────────────────────────────────
# AA positions counted from first ATG in reference (pos 1):
#   AA 2 = pos 4-6  (diagnostic: A→G → Glu/E)
#   AA 4 = pos 10-12 (diagnostic: A→C → Gln/Q)

test_that("aa pathway: strict match requires all AA mutations", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("aa pathway: single AA mutation assigns single-gene lineage", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  # NS2A_002: AA2=E only → NS4A; AA4=K (wildtype)
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "NS4A")
  # NS4B_003: AA4=Q only → NS4B; AA2=K (wildtype)
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "NS4B")
})

test_that("aa pathway: wildtype assigns unknown", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|OTHER_004"], "unknown")
})

test_that("aa pathway: most-specific lineage wins on multiple AA matches", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  # NY10_001 matches NY10 (2 AA muts), NS4A (1), and NS4B (1) → NY10 wins
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("aa pathway: fuzzy start codon does not block downstream AA lineage assignment", {
  result <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  # NST_NS2A_004: pos3=N → codon 1 = ATN → AA1 = X (fuzzy), but codon 2 = GAA → AA2 = E
  # The sequencing error at the start codon does not prevent NS4A assignment
  expect_equal(result$lineage[result$strain == "WNV|2019|CO|NST_NS2A_004"], "NS4A")
})

test_that("aa pathway assigns same lineages as nuc pathway (regression)", {
  result_aa <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(),
                   make_test_aa_muts(), type = "aa")
  )
  result_nuc <- suppressMessages(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts())
  )
  aa_lineages  <- result_aa$lineage[order(result_aa$strain)]
  nuc_lineages <- result_nuc$lineage[order(result_nuc$strain)]
  expect_equal(aa_lineages, nuc_lineages)
})
