#' WNV lineage-defining mutations (polyprotein/ORF coordinates)
#'
#' Diagnostic amino acid mutations used to assign West Nile virus lineages
#' with [define_lineage()] (`mut_type = "aa"`). `pos` is an amino acid
#' position counted from the polyprotein CDS start (position 1 = the first
#' Met), so this table can be passed to `define_lineage()` directly — no
#' [resolve_mut_positions()] step needed.
#'
#' Includes the ancestral WN02 marker plus lineage-specific markers for NY10
#' and SW03, and single-marker "partial" lineages (e.g. `NY10_NS2A`) that
#' catch sequences carrying only one of a lineage's markers.
#'
#' @format A tibble with 16 rows and 3 columns:
#' \describe{
#'   \item{lineage}{Lineage name (chr).}
#'   \item{pos}{Amino acid position from the polyprotein CDS start (int).}
#'   \item{residue}{Required amino acid residue at `pos` (chr).}
#' }
#' @seealso [wnv_mut_gff] for the same mutations in gene-relative coordinates.
"wnv_mut_orf"

#' WNV lineage-defining mutations (gene-relative coordinates)
#'
#' The same diagnostic mutations as [wnv_mut_orf], expressed as gene-relative
#' amino acid positions (e.g. NS2A position 188) rather than polyprotein
#' coordinates. Resolve to alignment coordinates first with
#' [resolve_mut_positions()] (supplying the genome GFF) before passing to
#' [define_lineage()] — see that function's documentation for the workflow.
#'
#' @format A tibble with 16 rows and 4 columns:
#' \describe{
#'   \item{lineage}{Lineage name (chr).}
#'   \item{gene}{Gene the mutation falls in, matching GFF product names (chr).}
#'   \item{pos}{Amino acid position from the gene's own start (int).}
#'   \item{residue}{Required amino acid residue at `pos` (chr).}
#' }
#' @seealso [wnv_mut_orf] for the same mutations in polyprotein coordinates.
"wnv_mut_gff"
