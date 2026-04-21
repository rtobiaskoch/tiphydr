#' Assign lineage labels to sequences based on diagnostic mutations
#'
#' Requires the alignment to have been produced by fasta_trim_ref() so that
#' alignment position N equals ungapped reference position N. Positions in
#' `muts$pos` must use this same reference coordinate system.
#'
#' Assignment is strict: all mutations for a lineage must be present. If a
#' sequence matches multiple lineages, the one with the most required mutations
#' (most specific) is assigned. Ties are broken alphabetically with a warning.
#'
#' @param alignment DNAStringSet; pre-aligned sequences, all equal length.
#'   Must be output of fasta_trim_ref() to satisfy coordinate contract.
#' @param ref length-1 DNAStringSet; reference sequence. Used to locate the
#'   first ATG start codon when type = "aa".
#' @param muts dataframe with columns: lineage (chr), pos (int), residue (chr).
#'   pos is in reference coordinates (ungapped).
#' @param type "nuc" (default) or "aa"
#' @param verbose logical; if TRUE prints lineage count table (default TRUE)
#'
#' @return dataframe with columns strain (chr) and lineage (chr)
#' @export
define_lineage <- function(alignment, ref, muts, type = "nuc", verbose = TRUE) {

  # ── Input validation ─────────────────────────────────────────────────────────

  # Ensure muts has the expected columns for lineage assignment
  required_cols <- c("lineage", "pos", "residue")
  if (!all(required_cols %in% names(muts))) {
    stop("muts must have columns: ", paste(required_cols, collapse = ", "))
  }

  # Only nucleotide and amino acid modes are supported
  if (!type %in% c("nuc", "aa")) stop("type must be 'nuc' or 'aa'")

  # Sequence names become the `strain` column in the output
  seq_names <- names(alignment)

  # ── Build residue lookup function ────────────────────────────────────────────
  # Returns the base/residue at a given position for a given sequence index.
  # Two modes:
  #   "nuc" — alignment position == reference position (coordinate contract)
  #   "aa"  — translate CDS (starting at first ATG in ref) then index by codon

  if (type == "aa") {
    # Find first ATG in reference to establish CDS start
    atg_matches <- Biostrings::matchPattern("ATG", ref[[1]])
    if (length(atg_matches) == 0L) stop("No ATG start codon found in ref")
    cds_start <- Biostrings::start(atg_matches)[1]

    # Translate all sequences from CDS start; fuzzy codons become "X"
    cds_seqs   <- Biostrings::subseq(alignment, start = cds_start)
    translated <- Biostrings::translate(cds_seqs, if.fuzzy.codon = "X")

    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(translated[seq_idx])) return(NA_character_)
      as.character(Biostrings::subseq(translated[seq_idx], pos, pos))
    }

  } else {
    # nuc: alignment position = reference position (coordinate contract from fasta_trim_ref)
    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(alignment[seq_idx])) return(NA_character_)
      toupper(as.character(Biostrings::subseq(alignment[seq_idx], pos, pos)))
    }
  }

  # ── Group mutations by lineage ───────────────────────────────────────────────
  # Split muts into a named list: lineage name → dataframe of (pos, residue) rows
  lineage_muts <- split(
    muts[, c("pos", "residue"), drop = FALSE],
    muts$lineage
  )

  # ── Assign lineage to each sequence ─────────────────────────────────────────
  # For each sequence, test all lineages and keep only those where every
  # required mutation is present. If multiple match, prefer the one with the
  # most mutations (most specific). Ties → alphabetical with warning.

  assign_one <- function(seq_idx) {

    # Collect lineages where ALL required mutations are satisfied
    matched <- purrr::keep(names(lineage_muts), function(lin) {
      lin_muts <- lineage_muts[[lin]]
      all(purrr::map2_lgl(
        lin_muts$pos, lin_muts$residue,
        function(p, r) {
          residue_val <- get_residue(seq_idx, p)
          !is.na(residue_val) && residue_val == toupper(r)
        }
      ))
    })

    if (length(matched) == 0L) return("unknown")
    if (length(matched) == 1L) return(matched)

    # Multiple matches: pick the lineage with the most required mutations (most specific)
    n_muts      <- purrr::map_int(matched, ~ nrow(lineage_muts[[.x]]))
    max_n       <- max(n_muts)
    candidates  <- matched[n_muts == max_n]

    if (length(candidates) > 1L) {
      warning(
        "Sequence '", seq_names[seq_idx],
        "' matched lineages with equal specificity: ",
        paste(sort(candidates), collapse = ", "),
        ". Assigning first alphabetically."
      )
      candidates <- sort(candidates)[1]
    }

    candidates
  }

  # Apply assignment across all sequences; returns character vector of lineage names
  result_lineages <- purrr::map_chr(seq_along(alignment), assign_one)

  # ── Assemble output dataframe ────────────────────────────────────────────────
  result <- data.frame(
    strain  = seq_names,
    lineage = result_lineages,
    stringsAsFactors = FALSE
  )

  # ── Optional verbose summary ─────────────────────────────────────────────────
  if (verbose) {
    message("\nLineage assignments:")
    print(table(result$lineage))

    n_unknown <- sum(result$lineage == "unknown")
    if (n_unknown > 0L) {
      message("\n", n_unknown, " unknown sequences:")
      message(paste(result$strain[result$lineage == "unknown"], collapse = "\n"))
    }
  }

  result
}
