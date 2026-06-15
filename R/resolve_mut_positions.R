#' Resolve gene-relative amino acid positions to absolute genomic codon starts
#'
#' Converts a mutation table that uses gene names and amino acid positions into
#' one that carries the absolute genomic nucleotide position of each codon. This
#' is a preprocessing step for [define_lineage()]: supply `cds_start` to compute
#' the `aa_start` column, then pass it as `pos` with `mut_type = "aa"`, e.g.
#' `resolve_mut_positions(muts, gff, cds_start = 97L) |>`
#' `dplyr::mutate(pos = aa_start) |> define_lineage(aln, ref, mut_type = "aa")`.
#'
#' Rows where `gene == "nuc"` are passed through unchanged with
#' `comparison_type = "nuc"` and `codon_start = NA`. All other rows are treated
#' as amino acid positions within the named gene, resolved via the GFF.
#'
#' Gene matching is case-insensitive substring search against the `product=`
#' attribute of `mature_protein_region_of_CDS` features in the GFF. Products
#' whose names contain `"/"` (frameshift/fusion products such as `N-NS4B/WARF4`)
#' are deprioritised so that queries like `"NS4B"` resolve to the canonical
#' protein. An informative error is raised if a query matches zero candidates or
#' more than one candidate after deprioritisation.
#'
#' @param muts data frame with columns: `lineage` (chr), `gene` (chr),
#'   `pos` (int), `residue` (chr). `pos` is an amino acid position within `gene`
#'   (1-indexed) for non-nuc rows, or an absolute nucleotide position for
#'   `gene == "nuc"` rows.
#' @param gff_path length-1 character; path to a GFF3 file.
#' @param cds_start integer (optional); the genomic nucleotide position of the
#'   first base of the CDS (i.e. the A in the polyprotein start ATG), in the
#'   same coordinate space as the GFF. When supplied, an `aa_start` column is
#'   computed: the position of the codon in the translated sequence, allowing
#'   direct comparison with [find_residue()] output.
#'   Formula: `ceiling((codon_start - cds_start) / 3 + 1)`.
#'   For WNV genomes with a 96-nt 5' UTR, use `cds_start = 97L`.
#'
#' @return tibble with the original four columns plus:
#'   \describe{
#'     \item{comparison_type}{`"nuc"` or `"aa"`}
#'     \item{codon_start}{integer genomic position of the first nucleotide of
#'       the codon; `NA_integer_` for `"nuc"` rows}
#'     \item{aa_start}{integer position of the codon in the translated sequence
#'       (counts from 1 at the CDS start); `NA_integer_` for `"nuc"` rows or
#'       when `cds_start` is not supplied}
#'   }
#'
#' @export
#' @importFrom stringr str_match str_detect fixed
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
resolve_mut_positions <- function(muts, gff_path, cds_start = NULL) {

  # ── Input validation ─────────────────────────────────────────────────────────

  if (!is.null(cds_start)) {
    if (!is.numeric(cds_start) || length(cds_start) != 1L || is.na(cds_start)) {
      stop("cds_start must be a single integer: the first nucleotide position of the CDS in GFF coordinates.")
    }
    cds_start <- as.integer(cds_start)
  }

  required_cols <- c("lineage", "gene", "pos", "residue")
  if (!all(required_cols %in% names(muts))) {
    stop("muts must have columns: ", paste(required_cols, collapse = ", "))
  }

  if (!file.exists(gff_path)) {
    stop("GFF file not found: ", gff_path)
  }

  # ── Parse GFF ────────────────────────────────────────────────────────────────
  # Read raw lines, drop comment/blank lines, filter to mature protein regions,
  # extract product name and nucleotide start position.

  raw_lines  <- readLines(gff_path, warn = FALSE)
  data_lines <- raw_lines[!startsWith(raw_lines, "#") & nchar(trimws(raw_lines)) > 0]

  if (length(data_lines) == 0L) {
    stop("GFF file contains no data rows (all lines are comments or blank): ", gff_path)
  }

  fields <- strsplit(data_lines, "\t", fixed = TRUE)

  # Keep only mature_protein_region_of_CDS rows with 9 fields
  is_mature <- vapply(fields, function(f) {
    length(f) >= 9L && f[[3L]] == "mature_protein_region_of_CDS"
  }, logical(1L))

  if (!any(is_mature)) {
    stop("No mature_protein_region_of_CDS features found in: ", gff_path)
  }

  mature_fields <- fields[is_mature]

  # Extract product= value from the attributes column (field 9)
  # using a capture group to avoid lookbehind portability issues
  attributes <- vapply(mature_fields, `[[`, character(1L), 9L)
  products   <- stringr::str_match(attributes, "product=([^;]+)")[, 2L]
  nt_starts  <- as.integer(vapply(mature_fields, `[[`, character(1L), 4L))

  gff_tbl <- tibble::tibble(product = products, nt_start = nt_starts)
  # Drop any rows where product could not be extracted
  gff_tbl <- gff_tbl[!is.na(gff_tbl$product), ]

  # ── Per-row resolution ───────────────────────────────────────────────────────

  resolve_one <- function(gene_val, pos_val) {

    # Rows explicitly marked as absolute nucleotide positions
    if (tolower(gene_val) == "nuc") {
      return(list(comparison_type = "nuc", codon_start = NA_integer_))
    }

    # Case-insensitive fixed substring match against all product names
    pattern_lower  <- tolower(gene_val)
    products_lower <- tolower(gff_tbl$product)
    match_idx      <- which(stringr::str_detect(products_lower, stringr::fixed(pattern_lower)))

    if (length(match_idx) == 0L) {
      stop(
        "Gene '", gene_val, "' not found in GFF. ",
        "Available products: ", paste(gff_tbl$product, collapse = ", ")
      )
    }

    # Deprioritise frameshift/fusion products (names containing "/")
    canonical_idx <- match_idx[!stringr::str_detect(gff_tbl$product[match_idx], "/")]
    candidate_idx <- if (length(canonical_idx) > 0L) canonical_idx else match_idx

    if (length(candidate_idx) > 1L) {
      stop(
        "Gene '", gene_val, "' matches multiple products: ",
        paste(gff_tbl$product[candidate_idx], collapse = ", "),
        ". Use a more specific name."
      )
    }

    # codon_start: first nucleotide of the codon for this amino acid position
    gene_nt_start <- gff_tbl$nt_start[candidate_idx[1L]]
    codon_start   <- as.integer(gene_nt_start + (pos_val - 1L) * 3L)

    list(comparison_type = "aa", codon_start = codon_start)
  }

  resolved <- purrr::map_dfr(seq_len(nrow(muts)), function(i) {
    r <- resolve_one(muts$gene[i], muts$pos[i])
    tibble::tibble(comparison_type = r$comparison_type, codon_start = r$codon_start)
  })

  # ── Compute aa_start ─────────────────────────────────────────────────────────
  # Position of the codon in the translated sequence (1-indexed from cds_start).
  # Only defined for AA rows when cds_start is supplied; NA otherwise.

  aa_start <- if (!is.null(cds_start)) {
    ifelse(
      resolved$comparison_type == "aa",
      as.integer(ceiling((resolved$codon_start - cds_start) / 3 + 1L)),
      NA_integer_
    )
  } else {
    rep(NA_integer_, nrow(resolved))
  }

  # ── Assemble output ──────────────────────────────────────────────────────────

  tibble::tibble(
    lineage         = muts$lineage,
    gene            = muts$gene,
    pos             = muts$pos,
    residue         = muts$residue,
    comparison_type = resolved$comparison_type,
    codon_start     = resolved$codon_start,
    aa_start        = aa_start
  )
}
