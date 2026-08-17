# Generates data/wnv_mut_gff.rda and data/wnv_mut_orf.rda.
#
# These are diagnostic-mutation lookup tables for define_lineage(). Each
# lineage is an explicit, literal set of required mutations (all rows must
# match) — including cross-family combinations such as "NY10-SW03_NS5" for
# sequences carrying markers from both the NY10 and SW03 families at once.
# See docs/superpowers/plans/staged-wondering-meerkat.md for the design.
#
# Markers are defined once (position + residue, in both gff gene-relative
# and orf polyprotein-relative coordinates) and lineages are built by
# selecting which markers they require, so the combinatorial expansion below
# is composed from these building blocks rather than retyping raw numbers.

# ── Marker definitions ───────────────────────────────────────────────────────
# env_A: the WN02-clade marker, shared by WN02/SW03/NY10 and all their
#   sub-lineages. env_V: the ancestral NY99 marker (mutually exclusive with
#   env_A — a sequence can only carry one).
markers <- tibble::tribble(
  ~marker , ~gene  , ~gff_pos , ~orf_pos , ~residue ,
  "env_A" , "env"  , 159L     ,  449L    , "A"      ,
  "env_V" , "env"  , 159L     ,  449L    , "V"      ,
  "NS4A"  , "NS4A" ,  85L     , 2209L    , "T"      ,
  "NS5"   , "NS5"  , 314L     , 2842L    , "R"      ,
  "NS2A"  , "NS2A" , 188L     , 1331L    , "K"      ,
  "NS4B"  , "NS4B" , 240L     , 2513L    , "M"
)

# ── Lineage definitions ──────────────────────────────────────────────────────
# Each lineage lists the markers it requires. Single-family lineages/partials
# (SW03, SW03_NS4A, ...) are the pre-existing set; the 9 cross-family entries
# are new, one per presence/absence state of (NS4A, NS5, NS2A, NS4B) not
# already covered by a single-family lineage.
#
# Naming: cross-family names join the two family pieces with "-", ordering
# by number of required mutations DESCENDING (the more specific piece first),
# falling back to alphabetical when the two pieces have an equal count — e.g.
# SW03 (3 muts) + NY10_NS2A (2 muts) -> "SW03-NY10_NS2A", but NY10_NS4B (2) +
# SW03_NS4A (2) -> "NY10_NS4B-SW03_NS4A" (tied, alphabetical).
lineages <- list(
  # ── base lineages (unchanged) ──
  NY99 = c("env_V"),
  WN02 = c("env_A"),
  SW03 = c("env_A", "NS4A", "NS5"),
  SW03_NS4A = c("env_A", "NS4A"),
  SW03_NS5 = c("env_A", "NS5"),
  NY10 = c("env_A", "NS2A", "NS4B"),
  NY10_NS2A = c("env_A", "NS2A"),
  NY10_NS4B = c("env_A", "NS4B"),

  # ── cross-family combinations (new) ──
  "NY10_NS4B-SW03_NS4A" = c("env_A", "NS4A", "NS4B"), # 2 vs 2: tied -> alphabetical
  "NY10_NS2A-SW03_NS4A" = c("env_A", "NS4A", "NS2A"), # 2 vs 2: tied -> alphabetical
  "NY10-SW03_NS4A" = c("env_A", "NS4A", "NS2A", "NS4B"), # 3 vs 2: NY10 more specific
  "NY10_NS4B-SW03_NS5" = c("env_A", "NS5", "NS4B"), # 2 vs 2: tied -> alphabetical
  "NY10_NS2A-SW03_NS5" = c("env_A", "NS5", "NS2A"), # 2 vs 2: tied -> alphabetical
  "NY10-SW03_NS5" = c("env_A", "NS5", "NS2A", "NS4B"), # 3 vs 2: NY10 more specific
  "SW03-NY10_NS4B" = c("env_A", "NS4A", "NS5", "NS4B"), # 3 vs 2: SW03 more specific
  "SW03-NY10_NS2A" = c("env_A", "NS4A", "NS5", "NS2A"), # 3 vs 2: SW03 more specific
  "NY10-SW03" = c("env_A", "NS4A", "NS5", "NS2A", "NS4B") # 3 vs 3: tied -> alphabetical
)

# ── Major lineage collapse ───────────────────────────────────────────────────
# SW03 = NS4A/NS5, NY10 = NS2A/NS4B are the family-defining markers; env_A/
# env_V are shared (WN02 clade) or ancestral (NY99) and don't count toward the
# SW03-vs-NY10 tally. A lineage collapses to whichever family has more
# markers present. Genuine ties (equal nonzero counts — the NY10/SW03 hybrid
# lineages split 1-1 or 2-2) are left as their own lineage name rather than
# forced to a side, since those cases are ambiguous by construction.
classify_major_lineage <- function(marker_names, lineage_name) {
  sw03_n <- sum(marker_names %in% c("NS4A", "NS5"))
  ny10_n <- sum(marker_names %in% c("NS2A", "NS4B"))

  dplyr::case_when(
    sw03_n == 0 & ny10_n == 0 & "env_V" %in% marker_names ~ "NY99",
    sw03_n == 0 & ny10_n == 0 ~ "WN02",
    sw03_n > ny10_n ~ "SW03",
    ny10_n > sw03_n ~ "NY10",
    TRUE ~ lineage_name
  )
}

major_lineages <- purrr::imap_chr(lineages, classify_major_lineage)

# ── Expand into the two shipped tables ───────────────────────────────────────

lineage_rows <- purrr::imap_dfr(lineages, function(marker_names, lineage_name) {
  dplyr::mutate(
    markers[match(marker_names, markers$marker), ],
    lineage = lineage_name,
    major_lineage = major_lineages[[lineage_name]],
    .before = 1
  )
})

wnv_mut_gff <- dplyr::select(
  lineage_rows,
  lineage,
  major_lineage,
  gene,
  pos = gff_pos,
  residue
)
wnv_mut_orf <- dplyr::select(
  lineage_rows,
  lineage,
  major_lineage,
  pos = orf_pos,
  residue
)

usethis::use_data(wnv_mut_gff, wnv_mut_orf, overwrite = TRUE)
