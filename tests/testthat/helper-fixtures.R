# helper-fixtures.R — shared in-memory test objects for all testthat tests
# Loaded automatically by testthat before any test file runs.

# ---------------------------------------------------------------------------
# Private reference sequence (dot-prefix keeps it out of test namespaces)
# ---------------------------------------------------------------------------

# 18nt synthetic reference: ATG start + 4 AAA codons + TAA stop.
# Layout (all positions 1-indexed):
#   pos  1-3  : ATG  — Met, start codon
#   pos  4-6  : AAA  — Lys (codon 2, gene_A) — DIAGNOSTIC: A→G gives GAA=Glu(E)
#   pos  7-9  : AAA  — Lys (codon 3, gene_A) — non-diagnostic
#   pos 10-12 : AAA  — Lys (codon 1, gene_B) — DIAGNOSTIC: A→C gives CAA=Gln(Q)
#   pos 13-15 : AAA  — Lys (codon 2, gene_B) — non-diagnostic
#   pos 16-18 : TAA  — stop (outside both genes in GFF)
.ref_seq <- paste0("ATG", "AAA", "AAA", "AAA", "AAA", "TAA")

# Sanity-check at load time so any accidental edit is caught immediately
stopifnot(nchar(.ref_seq) == 18L)

# ---------------------------------------------------------------------------
# Private helper: substitute bases at the two diagnostic positions
# ---------------------------------------------------------------------------

#' Build a variant of .ref_seq by substituting at positions 3, 4, and/or 10.
#'
#' @param base  Character scalar — template sequence (default: .ref_seq).
#' @param pos3  Single character — base at position 3 (wt = "G"; default keeps ATG start intact).
#' @param pos4  Single character — base at position 4 (wt = "A").
#' @param pos10 Single character — base at position 10 (wt = "A").
#' @return Character scalar of the same length as `base`.
.mutate_seq <- function(base = .ref_seq, pos4 = "A", pos10 = "A", pos3 = "G") {
  s <- base
  substr(s, 3, 3) <- pos3
  substr(s, 4, 4) <- pos4
  substr(s, 10, 10) <- pos10
  s
}

# ---------------------------------------------------------------------------
# Public factory functions — called by individual test files
# ---------------------------------------------------------------------------

#' In-memory DNAStringSet with 6 test sequences.
#'
#' Sequence layout (diagnostic positions 4 and 10):
#'   NY10_001      : pos4=G / pos10=C — both mutations; NY10 lineage
#'   NS2A_002      : pos4=G / pos10=A — pos4 mutation only; NS4A lineage
#'   NS4B_003      : pos4=A / pos10=C — pos10 mutation only; NS4B lineage
#'   OTHER_004     : pos4=A / pos10=A — no mutations; wildtype / unknown
#'   NST_NS2A_004  : pos3=N / pos4=G  — fuzzy start codon (ATN); NS4A lineage from AA2 downstream
#'   EXTRA_NY10_005: pos4=G / pos10=C — both mutations; N at pos14 (non-diagnostic); NY10 lineage
#'
#' NST_NS2A_004 simulates a sequencing error that corrupts the start codon.
#' The pipeline should still call NS4A from the downstream diagnostic position.
#' EXTRA_NY10_005 is intentionally absent from make_test_metadata().
#'
#' @return Biostrings::DNAStringSet with names in pipe-delimited format
#'   \code{VIRUS|YEAR|STATE|STRAIN_ID}
make_test_biostring <- function() {
  # NY10_001: both diagnostic mutations — NY10 lineage
  ny10 <- .mutate_seq(pos4 = "G", pos10 = "C")

  # NS2A_002: pos4=G only — NS4A lineage
  ns2a <- .mutate_seq(pos4 = "G")

  # NS4B_003: pos10=C only — NS4B lineage
  ns4b <- .mutate_seq(pos10 = "C")

  # OTHER_004: no mutations — wildtype, unknown
  other <- .mutate_seq()

  # NST_NS2A_004: fuzzy start codon (pos3=N → ATN) + pos4=G — tests that sequencing error at
  # start codon does not block lineage assignment from downstream diagnostic position
  nstart_ns2a <- .mutate_seq(pos3 = "N", pos4 = "G")

  # EXTRA_NY10_005: both mutations + N at pos14 (middle of codon 5, non-diagnostic) — NY10 lineage
  ns4b_extra <- .mutate_seq(pos4 = "G", pos10 = "C")
  substr(ns4b_extra, 14, 14) <- "N"

  Biostrings::DNAStringSet(c(
    `WNV|2021|CO|NY10_001` = ny10,
    `WNV|2020|CO|NS2A_002` = ns2a,
    `WNV|2022|WY|NS4B_003` = ns4b,
    `WNV|2019|CO|OTHER_004` = other,
    `WNV|2019|CO|NST_NS2A_004` = nstart_ns2a,
    `WNV|2021|NM|EXTRA_NY10_005` = ns4b_extra
  ))
}

#' In-memory reference DNAStringSet (length 1).
#'
#' @return Biostrings::DNAStringSet, single sequence named "WNV_REF_KU877344"
make_test_ref <- function() {
  Biostrings::DNAStringSet(c(WNV_REF_KU877344 = .ref_seq))
}

#' Metadata tibble — 4 rows, deliberately excludes NST_NS2A_004 and EXTRA_NY10_005.
#'
#' Columns: strain_id (chr), date (Date), location (chr), host (chr).
#' strain_id values match the 4th pipe-field of the DNAStringSet names
#' returned by \code{make_test_biostring()}.
#'
#' @return tibble::tibble with 4 rows
make_test_metadata <- function() {
  tibble::tibble(
    strain_id = c("NY10_001", "NS2A_002", "NS4B_003", "OTHER_004"),
    date = as.Date(c("2021-07-15", "2020-08-01", "2022-06-20", "2019-05-10")),
    location = c("Colorado", "Colorado", "Wyoming", "Colorado"),
    host = c("Culex tarsalis", "Culex pipiens", "Culex tarsalis", "Bird")
  )
}

#' Nucleotide mutation definitions tibble — 4 rows covering three lineages.
#'
#' NY10 requires pos4=G AND pos10=C (strict 2-mutation lineage).
#' NS4A requires pos4=G only         (single-position lineage).
#' NS4B requires pos10=C only        (single-position lineage).
#'
#' Wildtype residues (reference): pos4=A (Lys codon AAA), pos10=A (Lys codon AAA).
#'
#' Columns: lineage (chr), pos (int), residue (chr).
#'
#' @return tibble::tibble with 4 rows
make_test_muts <- function() {
  tibble::tibble(
    lineage = c("NY10", "NY10", "NS4A", "NS4B"),
    pos = c(4L, 10L, 4L, 10L),
    residue = c("G", "C", "G", "C")
  )
}

#' AA-position mutation definitions tibble — for the type = "aa" pathway.
#'
#' AA positions are counted from the CDS start (first ATG in ref = pos 1).
#'
#' Translation of .ref_seq:
#'   AA 1 : ATG (pos 1-3)  → Met  (M) — start codon
#'   AA 2 : pos 4-6        → Lys  (K) wt / Glu (E) when pos4=G — DIAGNOSTIC
#'   AA 3 : pos 7-9        → Lys  (K) — non-diagnostic
#'   AA 4 : pos 10-12      → Lys  (K) wt / Gln (Q) when pos10=C — DIAGNOSTIC
#'
#' NY10 requires AA2=E AND AA4=Q (both diagnostic positions).
#' NS4A requires AA2=E only.
#' NS4B requires AA4=Q only.
#'
#' @return tibble::tibble with columns: lineage, pos, residue
make_test_aa_muts <- function() {
  tibble::tibble(
    lineage = c("NY10", "NY10", "NS4A", "NS4B"),
    pos = c(2L, 4L, 2L, 4L),
    residue = c("E", "Q", "E", "Q")
  )
}

#' Gene-relative mutation definitions tibble — mirrors make_test_muts() at the AA level.
#'
#' Uses two synthetic genes defined in testdata/test_genome.gff:
#'   gene_A: nt 1–9  (3 codons: Met, Lys, Lys)
#'   gene_B: nt 10–15 (2 codons: Lys, Lys)
#'
#' Arithmetic:
#'   gene_A start=1;  AA 2 → codon_start = 1  + (2-1)*3 = 4   (== nuc pos=4  in make_test_muts)
#'   gene_B start=10; AA 1 → codon_start = 10 + (1-1)*3 = 10  (== nuc pos=10 in make_test_muts)
#'
#' Codons and residues in mutant sequences:
#'   pos4=G  → codon GAA = Glu (E)   pos4=A  → codon AAA = Lys (wildtype)
#'   pos10=C → codon CAA = Gln (Q)   pos10=A → codon AAA = Lys (wildtype)
#'
#' NY10 requires gene_A AA2=E AND gene_B AA1=Q.
#' NS4A requires gene_A AA2=E only.
#' NS4B requires gene_B AA1=Q only.
#'
#' @return tibble::tibble with columns: lineage, gene, pos, residue
make_test_gene_muts <- function() {
  tibble::tibble(
    lineage = c("NY10", "NY10", "NS4A", "NS4B"),
    gene = c("gene_A", "gene_B", "gene_A", "gene_B"),
    pos = c(2L, 1L, 2L, 1L),
    residue = c("E", "Q", "E", "Q")
  )
}

#' Sequence with an ambiguous base AT a diagnostic position (pos4 = N).
#'
#' pos4 is diagnostic for NY10 and NS4A; N there means the lineage cannot be
#' determined. pos10 = A (wildtype) so NS4B is definitively ruled out.
#' Expected output from define_lineage(): "ambiguous".
#'
#' @return Biostrings::DNAStringSet, single sequence named "WNV|2023|CO|AMB_001"
make_test_ambiguous_seq <- function() {
  seq <- .mutate_seq() # wildtype: all A at diagnostic positions
  substr(seq, 4, 4) <- "N" # N at diagnostic pos4
  Biostrings::DNAStringSet(c(`WNV|2023|CO|AMB_001` = seq))
}

#' Minimal in-memory phylo object using the same tip labels as make_test_biostring().
#'
#' Tips are the 4 fully-annotated WNV pipe-delimited headers.
#' Topology is arbitrary — tests care about label identity, not branch lengths.
#'
#' @return ape::phylo
make_test_tree <- function() {
  ape::read.tree(
    text = "((WNV|2021|CO|NY10_001,WNV|2020|CO|NS2A_002),(WNV|2022|WY|NS4B_003,WNV|2019|CO|OTHER_004));"
  )
}

# ---------------------------------------------------------------------------
# Introduction-detection fixtures (explode_tree and helpers)
# ---------------------------------------------------------------------------

#' Small dated tree with a known introduction structure.
#'
#' Topology: ((tB1,tB2),(tA1,tC1)). Tip labels are pipe-delimited
#' \code{strain|date|deme}. Branch lengths increase toward the tips so internal
#' nodes are dated earlier than tips (a usable molecular clock).
#'
#' Intended deme states (set in make_intro_node_probs):
#'   root            -> regionA
#'   MRCA(tB1,tB2)   -> regionB  (an ESTABLISHED introduction, clade size 2)
#'   MRCA(tA1,tC1)   -> regionA
#'   tC1             -> regionC  (a SINGLETON introduction, clade size 1)
#'   tA1             -> regionA  (no transition)
#'
#' @return ape::phylo with branch lengths.
make_intro_tree <- function() {
  # Branch lengths chosen so root-to-tip distance increases with sampling date
  # (a valid molecular clock): depth 4->2019 ... depth 7->2022, +1yr per unit.
  ape::read.tree(
    text = paste0(
      "((tB1|2020-01-01|regionB:2,tB2|2021-01-01|regionB:3):3,",
      "(tA1|2019-01-01|regionA:2,tC1|2022-01-01|regionC:5):2);"
    )
  )
}

#' Node x deme probability matrix matching make_intro_tree().
#'
#' Tips are one-hot on their own deme (probability 1). Internal node states are
#' assigned by MRCA so the fixture is robust to ape's node numbering.
#'
#' @param tree The phylo from make_intro_tree().
#' @return data.frame with a node column + one column per deme.
make_intro_node_probs <- function(tree = make_intro_tree()) {
  ntip <- ape::Ntip(tree)
  demes <- c("regionA", "regionB", "regionC")

  # Resolve internal node ids by MRCA (independent of ape numbering).
  lbl <- function(p) grep(p, tree$tip.label, value = TRUE)
  root <- ntip + 1L
  node_B <- ape::getMRCA(tree, c(lbl("^tB1"), lbl("^tB2")))
  node_AC <- ape::getMRCA(tree, c(lbl("^tA1"), lbl("^tC1")))

  # State per node id: tips from their label deme, internals as documented.
  tip_deme <- vapply(
    strsplit(tree$tip.label, "|", fixed = TRUE),
    `[`,
    character(1),
    3L
  )
  state <- character(ntip + tree$Nnode)
  state[seq_len(ntip)] <- tip_deme
  state[root] <- "regionA"
  state[node_B] <- "regionB"
  state[node_AC] <- "regionA"

  # One-hot probability matrix.
  mat <- matrix(
    0,
    nrow = length(state),
    ncol = length(demes),
    dimnames = list(NULL, demes)
  )
  mat[cbind(seq_along(state), match(state, demes))] <- 1

  data.frame(node = seq_along(state), mat, check.names = FALSE)
}

#' Tip-date tibble matching make_intro_tree() (parsed from the tip labels).
#'
#' @param tree The phylo from make_intro_tree().
#' @return tibble with tip_label and date (Date).
make_intro_tip_dates <- function(tree = make_intro_tree()) {
  tibble::tibble(
    tip_label = tree$tip.label,
    date = as.Date(vapply(
      strsplit(tree$tip.label, "|", fixed = TRUE),
      `[`,
      character(1),
      2L
    ))
  )
}

#' Path to a treeio-readable native-annotated NEXUS file matching make_intro_tree().
#'
#' Per-node dates follow the same "+1 year per depth unit" scheme documented in
#' make_intro_tree()'s own comment (depth 4 -> 2019 ... depth 7 -> 2022), just
#' extended to the 3 internal nodes: root = 2015-01-01, MRCA(tA1,tC1) =
#' 2017-01-01, MRCA(tB1,tB2) = 2018-01-01. Tip dates match
#' make_intro_tip_dates() exactly (hand-verified against tree depth) unless
#' overridden. Written via treeio::write.beast() -- round-trip verified via
#' treeio::read.beast() to produce exactly the NHX `[&date=...]` shape
#' node_dates_from_timetree() parses (including a BEGIN TAXA; block, like
#' TreeTime's real output, so this fixture exercises the same file shape
#' `treeio::read.beast()` sees in production, not a simplified stand-in).
#'
#' @param tree The phylo from make_intro_tree().
#' @param tip_dates_override Optional named character vector (tip label ->
#'   raw annotation string) overriding specific tip dates -- used to build a
#'   deliberately non-clocklike or non-ISO (bare-year / decimal-year) fixture.
#' @return Path to a temp NEXUS file.
make_intro_annotated_tree_path <- function(tree = make_intro_tree(), tip_dates_override = NULL) {
  ntip <- ape::Ntip(tree)
  lbl <- function(p) grep(p, tree$tip.label, value = TRUE)
  root <- ntip + 1L
  node_B <- ape::getMRCA(tree, c(lbl("^tB1"), lbl("^tB2")))
  node_AC <- ape::getMRCA(tree, c(lbl("^tA1"), lbl("^tC1")))

  tip_date_of <- function(label) {
    if (!is.null(tip_dates_override) && label %in% names(tip_dates_override)) {
      return(tip_dates_override[[label]])
    }
    strsplit(label, "|", fixed = TRUE)[[1]][2]
  }
  date_str <- character(ntip + tree$Nnode)
  for (i in seq_len(ntip)) date_str[i] <- tip_date_of(tree$tip.label[i])
  date_str[root] <- "2015-01-01"
  date_str[node_AC] <- "2017-01-01"
  date_str[node_B] <- "2018-01-01"

  td <- treeio::treedata(
    phylo = tree,
    data = tibble::tibble(node = seq_along(date_str), date = date_str)
  )
  tf <- tempfile(fileext = ".nex")
  treeio::write.beast(td, file = tf)
  tf
}

# ---------------------------------------------------------------------------
# File writer — single source of truth for on-disk fixtures
# ---------------------------------------------------------------------------

#' Write all test fixtures to disk.
#'
#' Called manually via tests/write_fixtures.R — never called during test runs.
#' Outputs mirror the in-memory make_*() functions defined above so that I/O
#' tests can reference real files without diverging from unit-test data.
#'
#' @param out_dir Path to output directory (created if absent).
#' @return out_dir invisibly.
write_test_fixtures <- function(out_dir = "tests/testthat/testdata") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # FASTA files — written via Biostrings
  Biostrings::writeXStringSet(
    make_test_ref(),
    file.path(out_dir, "wnv_ref.fasta")
  )
  Biostrings::writeXStringSet(
    make_test_biostring(),
    file.path(out_dir, "wnv_seqs.fasta")
  )

  # CSV files — metadata and mutation definitions
  utils::write.csv(
    make_test_metadata(),
    file.path(out_dir, "wnv_metadata.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    make_test_muts(),
    file.path(out_dir, "wnv_muts.csv"),
    row.names = FALSE
  )

  # Minimal GFF3 for gene-relative mutation tests.
  # Two genes covering the 18nt reference; stop codon (pos 16-18) is outside both.
  #   gene_A: pos  1-9  (Met + 2 Lys codons — AA2 is the diagnostic Glu/Lys position)
  #   gene_B: pos 10-15 (2 Lys codons     — AA1 is the diagnostic Gln/Lys position)
  make_gff_row <- function(start, end, product) {
    paste(
      c(
        "test_seq",
        ".",
        "mature_protein_region_of_CDS",
        start,
        end,
        ".",
        "+",
        ".",
        paste0("product=", product)
      ),
      collapse = "\t"
    )
  }
  gff_lines <- c(
    "##gff-version 3",
    make_gff_row(1, 9, "gene_A"),
    make_gff_row(10, 15, "gene_B")
  )
  writeLines(gff_lines, file.path(out_dir, "test_genome.gff"))

  message("Wrote 5 fixture files to ", out_dir)
  invisible(out_dir)
}
