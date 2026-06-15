test_that("returns dataframe with strain and lineage columns", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  expect_s3_class(result, "data.frame")
  expect_named(result, c("strain", "lineage"))
  expect_equal(nrow(result), length(make_test_biostring()))
})

test_that("strict match: NY10 requires both mutations", {
  # NY10_001 has pos4=G and pos10=C → matches NY10 (more specific than NS4A or NS4B)
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("single-mutation sequences assign their single-gene lineage", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  # NS2A_002: pos4=G only → NS4A
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "NS4A")
  # NS4B_003: pos10=C only → NS4B
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "NS4B")
})

test_that("wildtype sequence assigns unknown", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  expect_equal(
    result$lineage[result$strain == "WNV|2019|CO|OTHER_004"],
    "unknown"
  )
})

test_that("multiple-match: assigns most-specific (most mutations)", {
  # NY10_001 matches NY10 (2 muts), NS4A (1 mut), and NS4B (1 mut) → NY10 wins
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("muts missing required columns throws error", {
  bad_muts <- tibble::tibble(lineage = "NY10", position = 10L, base = "A")
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), bad_muts, mut_type = "nuc"),
    "lineage.*pos.*residue"
  )
})

test_that("invalid type throws error", {
  expect_error(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "protein"
    ),
    "type must be"
  )
})

test_that("AAStringSet input assigns same lineages as the DNA-translate path", {
  # User already has translated sequences; feed them directly with mut_type = "aa".
  # ref is not needed because no nucleotide CDS start has to be located.
  aa_aln <- Biostrings::translate(
    Biostrings::chartr(".-", "NN", make_test_biostring()),
    if.fuzzy.codon = "X"
  )

  result_aa_input <- suppressMessages(
    define_lineage(aa_aln, muts = make_test_aa_muts(), mut_type = "aa")
  )
  result_dna <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )

  aa_lineages <- result_aa_input$lineage[order(result_aa_input$strain)]
  dna_lineages <- result_dna$lineage[order(result_dna$strain)]
  expect_equal(aa_lineages, dna_lineages)
})

test_that("AAStringSet input with mut_type = 'nuc' errors", {
  aa_aln <- Biostrings::AAStringSet(c(s1 = "MKK", s2 = "MEK"))
  expect_error(
    define_lineage(aa_aln, muts = make_test_muts(), mut_type = "nuc"),
    "nucleotide alignment"
  )
})

test_that("nucleotide alignment without ref errors", {
  expect_error(
    define_lineage(make_test_biostring(), muts = make_test_muts(), mut_type = "nuc"),
    "ref is required"
  )
})

test_that("mut_type is required (errors when omitted)", {
  expect_error(
    define_lineage(make_test_biostring(), make_test_ref(), make_test_muts()),
    "mut_type is required"
  )
})

# ── resolve_mut_positions() → define_lineage() composition ──────────────────────
# Gene-relative AA mutations are resolved to alignment AA coordinates first, then
# fed in as mut_type = "aa". This replaces the former in-function gff pathway.

test_that("resolved gene muts assign same lineages as the nuc pathway (regression)", {
  # gene_A starts at nt 1 in the fixture GFF, so the CDS begins at position 1.
  resolved <- resolve_mut_positions(
    make_test_gene_muts(),
    test_path("testdata/test_genome.gff"),
    cds_start = 1L
  )
  # aa_start is the AA position from CDS start; feed it as `pos` for the aa path
  aa_muts <- dplyr::mutate(resolved, pos = aa_start)

  result_resolved <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      aa_muts,
      mut_type = "aa"
    )
  )
  result_nuc <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )

  # Align by strain name so order does not matter
  resolved_lineages <- result_resolved$lineage[order(result_resolved$strain)]
  nuc_lineages <- result_nuc$lineage[order(result_nuc$strain)]
  expect_equal(resolved_lineages, nuc_lineages)
})

# ── mut_type = "aa" pathway tests ─────────────────────────────────────────────────
# AA positions counted from first ATG in reference (pos 1):
#   AA 2 = pos 4-6  (diagnostic: A→G → Glu/E)
#   AA 4 = pos 10-12 (diagnostic: A→C → Gln/Q)

test_that("aa pathway: strict match requires all AA mutations", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("aa pathway: single AA mutation assigns single-gene lineage", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  # NS2A_002: AA2=E only → NS4A; AA4=K (wildtype)
  expect_equal(result$lineage[result$strain == "WNV|2020|CO|NS2A_002"], "NS4A")
  # NS4B_003: AA4=Q only → NS4B; AA2=K (wildtype)
  expect_equal(result$lineage[result$strain == "WNV|2022|WY|NS4B_003"], "NS4B")
})

test_that("aa pathway: wildtype assigns unknown", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  expect_equal(
    result$lineage[result$strain == "WNV|2019|CO|OTHER_004"],
    "unknown"
  )
})

test_that("aa pathway: most-specific lineage wins on multiple AA matches", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  # NY10_001 matches NY10 (2 AA muts), NS4A (1), and NS4B (1) → NY10 wins
  expect_equal(result$lineage[result$strain == "WNV|2021|CO|NY10_001"], "NY10")
})

test_that("aa pathway: fuzzy start codon does not block downstream AA lineage assignment", {
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  # NST_NS2A_004: pos3=N → codon 1 = ATN → AA1 = X (fuzzy), but codon 2 = GAA → AA2 = E
  # The sequencing error at the start codon does not prevent NS4A assignment
  expect_equal(
    result$lineage[result$strain == "WNV|2019|CO|NST_NS2A_004"],
    "NS4A"
  )
})

test_that("aa pathway assigns same lineages as nuc pathway (regression)", {
  result_aa <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa"
    )
  )
  result_nuc <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc"
    )
  )
  aa_lineages <- result_aa$lineage[order(result_aa$strain)]
  nuc_lineages <- result_nuc$lineage[order(result_nuc$strain)]
  expect_equal(aa_lineages, nuc_lineages)
})

# ── ambiguous lineage tests ────────────────────────────────────────────────────

test_that("nuc pathway: N at diagnostic position returns 'ambiguous'", {
  alignment <- c(make_test_biostring(), make_test_ambiguous_seq())
  result <- suppressWarnings(suppressMessages(
    define_lineage(
      alignment,
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc",
      verbose = FALSE
    )
  ))
  expect_equal(
    result$lineage[result$strain == "WNV|2023|CO|AMB_001"],
    "ambiguous"
  )
})

test_that("nuc pathway: N at non-diagnostic position does not trigger 'ambiguous'", {
  # NST_NS2A_004 has N at pos3 (non-diagnostic) → NS4A, not ambiguous
  # EXTRA_NY10_005 has N at pos14 (non-diagnostic) → NY10, not ambiguous
  result <- suppressMessages(
    define_lineage(
      make_test_biostring(),
      make_test_ref(),
      make_test_muts(),
      mut_type = "nuc",
      verbose = FALSE
    )
  )
  expect_equal(
    result$lineage[result$strain == "WNV|2019|CO|NST_NS2A_004"],
    "NS4A"
  )
  expect_equal(
    result$lineage[result$strain == "WNV|2021|NM|EXTRA_NY10_005"],
    "NY10"
  )
})

test_that("aa pathway: N at diagnostic nuc position (fuzzy codon) returns 'ambiguous'", {
  alignment <- c(make_test_biostring(), make_test_ambiguous_seq())
  result <- suppressWarnings(suppressMessages(
    define_lineage(
      alignment,
      make_test_ref(),
      make_test_aa_muts(),
      mut_type = "aa",
      verbose = FALSE
    )
  ))
  expect_equal(
    result$lineage[result$strain == "WNV|2023|CO|AMB_001"],
    "ambiguous"
  )
})

test_that("ambiguous sequence triggers a warning naming the sequence", {
  alignment <- c(make_test_biostring(), make_test_ambiguous_seq())
  expect_warning(
    suppressMessages(
      define_lineage(
        alignment,
        make_test_ref(),
        make_test_muts(),
        mut_type = "nuc",
        verbose = FALSE
      )
    ),
    regexp = "AMB_001"
  )
})
