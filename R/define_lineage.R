# IUPAC ambiguity codes: everything except the four unambiguous bases and gap chars
.iupac_ambiguous <- c("N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V")

#' Assign lineage labels to sequences based on diagnostic mutations
#'
#' Requires the alignment to have been produced by fasta_trim_ref() so that
#' alignment position N equals ungapped reference position N. Positions in
#' `muts$pos` must use this same reference coordinate system.
#'
#' Assignment is strict: all mutations for a lineage must be present. If a
#' sequence matches multiple lineages, the one with the most required mutations
#' (most specific) is assigned. A tie between two or more equally-specific
#' lineages has no principled winner, so it is not broken arbitrarily: the
#' tied lineage names are sorted and joined with "-" to form a composite
#' label (e.g. "NY10_NS2A-SW03_NS5"), with a warning naming the tie.
#'
#' For gene-relative amino acid mutations, resolve them to alignment coordinates
#' first with resolve_mut_positions(), then pass its `aa_start` column as `pos`
#' with mut_type = "aa". See that function's documentation for the workflow.
#'
#' @param alignment DNAStringSet (nucleotides) or AAStringSet (amino acids);
#'   pre-aligned sequences, all equal length. A nucleotide alignment must be the
#'   output of fasta_trim_ref() to satisfy the coordinate contract. An AAStringSet
#'   is only valid with mut_type = "aa" and is indexed directly without translation.
#' @param ref length-1 DNAStringSet; reference sequence. Required for a
#'   nucleotide alignment: used to locate the first ATG start codon when
#'   mut_type = "aa", and to sanity-check coordinates. Ignored (and may be NULL)
#'   when alignment is an AAStringSet.
#' @param muts dataframe with columns: lineage (chr), pos (int), residue (chr).
#'   pos is in reference coordinates (ungapped) for mut_type = "nuc", or an
#'   amino acid position from the CDS start for mut_type = "aa".
#' @param mut_type "nuc" or "aa" — what the mutations are defined as. Required;
#'   there is no default, so the caller must always state mutation type.
#' @param verbose logical; if TRUE prints lineage count table (default TRUE)
#'
#' @return dataframe with columns strain (chr) and lineage (chr)
#' @export
define_lineage <- function(
  alignment,
  ref = NULL,
  muts,
  mut_type,
  verbose = TRUE
) {
  # ── Input validation ─────────────────────────────────────────────────────────

  # mut_type is required: there is no safe default. Nucleotide and amino acid
  # mutation tables can look identical (e.g. residues A/K/M/T/V are valid in
  # both IUPAC and amino acid alphabets), so the caller must declare intent.
  if (missing(mut_type)) {
    stop("mut_type is required: specify 'nuc' or 'aa' for the mutation table.")
  }
  if (!mut_type %in% c("nuc", "aa")) {
    stop("mut_type must be 'nuc' or 'aa'")
  }

  required_cols <- c("lineage", "pos", "residue")
  if (!all(required_cols %in% names(muts))) {
    stop("muts must have columns: ", paste(required_cols, collapse = ", "))
  }

  # The alignment may be nucleotides (DNAStringSet) or already-translated amino
  # acids (AAStringSet). An AA alignment only makes sense with mut_type = "aa":
  # nucleotides cannot be recovered from amino acids.
  is_aa_input <- methods::is(alignment, "AAStringSet")
  if (mut_type == "nuc" && is_aa_input) {
    stop("mut_type = 'nuc' needs a nucleotide alignment (DNAStringSet), but an AAStringSet was supplied.")
  }

  # ref is only needed for a nucleotide alignment: to sanity-check coordinates and
  # (for mut_type = "aa") to locate the CDS start. An AA alignment ignores it.
  if (!is_aa_input && is.null(ref)) {
    stop("ref is required for a nucleotide alignment.")
  }

  # Coordinate sanity check only applies when comparing a nucleotide alignment to
  # a nucleotide reference.
  if (!is_aa_input && any(Biostrings::width(alignment) != Biostrings::width(ref))) {
    warning(
      "Genome lengths of input alignment do not match ref.
            run fasta_trim_ref(fasta, ref) Ensure genomes are aligned and trimmed."
    )
  }

  # Sequence names become the `strain` column in the output
  seq_names <- names(alignment)

  # ── Build residue lookup function ────────────────────────────────────────────
  # Returns the base/residue at a given position for a given sequence index.
  # Two modes:
  #   "nuc" — alignment position == reference position (coordinate contract)
  #   "aa"  — translate CDS (starting at first ATG in ref) then index by codon

  if (mut_type == "aa") {
    # DETERMINE STARTING POINT, then index amino acids -------------------------
    # The start codon establishes amino acid position 1 (the coordinate that
    # muts$pos is defined against). detect_sequence_start() finds its modal
    # position so leading UTR / pre-start residues are trimmed off.

    if (is_aa_input) {
      # Alignment is already amino acids — locate the start residue M (Met) and
      # index from there; no translation needed.
      start_info <- detect_sequence_start(alignment, pattern = "M", verbose = verbose)
      aa_start <- start_info$position

      if (aa_start != 1) {
        warning("Your amino acid sequences do not start with M (Met)")
      }

      translated <- Biostrings::subseq(alignment, start = aa_start)
    } else {
      # Nucleotide alignment — locate the ATG start codon in the reference, then
      # translate every sequence from that position.
      start_info <- detect_sequence_start(ref, pattern = "ATG", verbose = verbose)
      cds_start <- start_info$position

      if (cds_start != 1) {
        warning("Your reference sequences does not start with a start codon")
      }

      # Translate all sequences from CDS start; fuzzy codons become "X"
      cds_seqs <- Biostrings::subseq(alignment, start = cds_start)
      # Replace gap/ambiguous characters with N so translation treats them as fuzzy codons ("X")
      cds_seqs <- Biostrings::chartr(".-", "NN", cds_seqs)

      translated <- Biostrings::translate(cds_seqs, if.fuzzy.codon = "X")
    }

    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(translated[seq_idx])) {
        return(NA_character_)
      }
      as.character(Biostrings::subseq(translated[seq_idx], pos, pos))
    }
  } else {
    # nuc: alignment position = reference position (coordinate contract from fasta_trim_ref)
    get_residue <- function(seq_idx, pos) {
      if (pos > Biostrings::width(alignment[seq_idx])) {
        return(NA_character_)
      }
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
  assign_one <- function(seq_idx) {
    # Collect lineages where ALL required mutations are satisfied
    matched <- purrr::keep(names(lineage_muts), function(lin) {
      lin_muts <- lineage_muts[[lin]]
      all(purrr::map2_lgl(
        lin_muts$pos,
        lin_muts$residue,
        function(p, r) {
          residue_val <- get_residue(seq_idx, p)
          !is.na(residue_val) && residue_val == toupper(r)
        }
      ))
    })

    if (length(matched) == 0L) {
      # Return "ambiguous" when a diagnostic position has an unresolvable base.
      # nuc mode: check for IUPAC ambiguity codes; aa mode: check for "X" (fuzzy codon).
      any_ambiguous <- any(purrr::map_lgl(unique(muts$pos), function(p) {
        val <- get_residue(seq_idx, p)
        if (mut_type == "nuc") {
          !is.na(val) && val %in% .iupac_ambiguous
        } else {
          !is.na(val) && val == "X"
        }
      }))
      if (any_ambiguous) {
        return("ambiguous")
      }
      return("unknown")
    }
    if (length(matched) == 1L) {
      return(matched)
    }

    # Multiple matches: pick the lineage with the most required mutations (most specific)
    n_muts <- purrr::map_int(matched, ~ nrow(lineage_muts[[.x]]))
    max_n <- max(n_muts)
    candidates <- matched[n_muts == max_n]

    if (length(candidates) > 1L) {
      # A tie between equally-specific lineages has no principled winner, so it
      # is not resolved arbitrarily: the tied names are combined into a
      # composite label instead.
      hybrid_label <- paste(sort(candidates), collapse = "-")
      warning(
        "Sequence '",
        seq_names[seq_idx],
        "' matched lineages with equal specificity: ",
        paste(sort(candidates), collapse = ", "),
        ". Assigning composite lineage '", hybrid_label, "'."
      )
      return(hybrid_label)
    }
    candidates
  }

  # Apply assignment across all sequences; returns character vector of lineage names
  result_lineages <- purrr::map_chr(seq_along(alignment), assign_one)

  # ── Assemble output dataframe ────────────────────────────────────────────────
  result <- data.frame(
    strain = seq_names,
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
      message(paste(
        result$strain[result$lineage == "unknown"],
        collapse = "\n"
      ))
    }
  }

  # Always warn — ambiguous bases at diagnostic positions flag a data quality issue
  n_ambiguous <- sum(result$lineage == "ambiguous")
  if (n_ambiguous > 0L) {
    warning(
      n_ambiguous,
      " sequence(s) had ambiguous bases at lineage-defining positions:\n",
      paste(result$strain[result$lineage == "ambiguous"], collapse = "\n"),
      call. = FALSE
    )
  }

  result
}
