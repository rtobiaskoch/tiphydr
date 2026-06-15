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
#' (most specific) is assigned. Ties are broken alphabetically with a warning.
#'
#' @param alignment DNAStringSet; pre-aligned sequences, all equal length.
#'   Must be output of fasta_trim_ref() to satisfy coordinate contract.
#' @param ref length-1 DNAStringSet; reference sequence. Used to locate the
#'   first ATG start codon when mut_type = "aa".
#' @param muts dataframe with columns: lineage (chr), pos (int), residue (chr).
#'   pos is in reference coordinates (ungapped). Alternatively, a dataframe with
#'   columns lineage, gene, pos, residue for gene-relative amino acid mutations
#'   (requires gff to be supplied).
#' @param mut_type "nuc" (default) or "aa" for what the mutations are defined as.
#'   Ignored when muts has a gene column.
#' @param verbose logical; if TRUE prints lineage count table (default TRUE)
#' @param gff length-1 character path to a GFF3 file, or NULL (default). Required
#'   when muts has a gene column with non-"nuc" gene names.
#'
#' @return dataframe with columns strain (chr) and lineage (chr)
#' @export
define_lineage <- function(
  alignment,
  ref,
  muts,
  mut_type = "aa",
  verbose = TRUE,
  gff = NULL
) {
  # ── Input validation ─────────────────────────────────────────────────────────

  has_gene_col <- "gene" %in% names(muts)

  if (has_gene_col) {
    # Gene-relative pathway: muts must have lineage, gene, pos, residue
    required_cols <- c("lineage", "gene", "pos", "residue")
    if (!all(required_cols %in% names(muts))) {
      stop("muts must have columns: ", paste(required_cols, collapse = ", "))
    }
    non_nuc_genes <- unique(muts$gene[tolower(muts$gene) != "nuc"])
    if (length(non_nuc_genes) > 0L && is.null(gff)) {
      stop(
        "muts has gene-relative positions (",
        paste(non_nuc_genes, collapse = ", "),
        ") but no 'gff' path was supplied. Provide gff = path_to_gff."
      )
    }
    if (!is.null(gff)) {
      muts_resolved <- resolve_mut_positions(muts, gff)
    } else {
      # All rows are gene == "nuc" — fall through to existing nuc path
      has_gene_col <- FALSE
    }
  } else {
    # Ensure muts has the expected columns for the existing pathway
    required_cols <- c("lineage", "pos", "residue")
    if (!all(required_cols %in% names(muts))) {
      stop("muts must have columns: ", paste(required_cols, collapse = ", "))
    }
  } #end if has_gene_col

  # Only nucleotide and amino acid modes are supported (existing pathway)
  if (!mut_type %in% c("nuc", "aa")) {
    stop("mut_type must be 'nuc' or 'aa'")
  }

  if (any(Biostrings::width(alignment) != Biostrings::width(ref))) {
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
    # Find first ATG in reference to establish CDS start
    # Use detect_sequence_start() to find modal start position

    #DETERMINE STARTING POINT --------------------------------------------------

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

  # ── Assign lineage to each sequence ─────────────────────────────────────────
  # Two paths depending on whether muts carries a gene column:
  #
  #   has_gene_col == TRUE  — per-row dispatch: nucleotide or codon-AA comparison
  #                           chosen per row from muts_resolved$comparison_type
  #   has_gene_col == FALSE — existing path: one get_residue closure for the
  #                           whole call, governed by the `type` parameter

  if (has_gene_col) {
    # Per-row residue extractor: returns the nucleotide or translated amino acid
    # at the position specified by this mutation row for a given sequence.
    get_residue_row <- function(seq_idx, row) {
      if (row$comparison_type == "nuc") {
        pos <- row$pos
        if (pos > Biostrings::width(alignment[seq_idx])) {
          return(NA_character_)
        }
        toupper(as.character(Biostrings::subseq(alignment[seq_idx], pos, pos)))
      } else {
        # "aa": extract 3-nt codon, replace gap chars with N, translate
        cs <- row$codon_start
        if (cs + 2L > Biostrings::width(alignment[seq_idx])) {
          return(NA_character_)
        }
        codon_str <- gsub(
          "[.-]",
          "N",
          as.character(Biostrings::subseq(alignment[seq_idx], cs, cs + 2L))
        )
        unname(as.character(
          Biostrings::translate(
            Biostrings::DNAStringSet(codon_str),
            if.fuzzy.codon = "X"
          )
        ))
      }
    }

    lineage_muts_r <- split(muts_resolved, muts_resolved$lineage)

    assign_one <- function(seq_idx) {
      matched <- purrr::keep(names(lineage_muts_r), function(lin) {
        df <- lineage_muts_r[[lin]]
        all(purrr::map_lgl(seq_len(nrow(df)), function(ri) {
          val <- get_residue_row(seq_idx, df[ri, , drop = FALSE])
          !is.na(val) && val == toupper(df$residue[ri])
        }))
      })

      if (length(matched) == 0L) {
        # Return "ambiguous" when any diagnostic position has an unresolvable base.
        # comparison_type drives whether to check IUPAC codes (nuc) or "X" (aa).
        any_ambiguous <- any(purrr::map_lgl(
          seq_len(nrow(muts_resolved)),
          function(ri) {
            val <- get_residue_row(seq_idx, muts_resolved[ri, , drop = FALSE])
            comp <- muts_resolved$comparison_type[ri]
            if (comp == "nuc") {
              !is.na(val) && val %in% .iupac_ambiguous
            } else {
              !is.na(val) && val == "X"
            }
          }
        ))
        if (any_ambiguous) {
          return("ambiguous")
        }
        return("unknown")
      }
      if (length(matched) == 1L) {
        return(matched)
      }

      n_muts <- purrr::map_int(matched, ~ nrow(lineage_muts_r[[.x]]))
      candidates <- matched[n_muts == max(n_muts)]

      if (length(candidates) > 1L) {
        warning(
          "Sequence '",
          seq_names[seq_idx],
          "' matched lineages with equal specificity: ",
          paste(sort(candidates), collapse = ", "),
          ". Assigning first alphabetically."
        )
        candidates <- sort(candidates)[1L]
      }
      candidates
    }
  } else {
    #end has gene_col
    # ── Group mutations by lineage (existing pathway) ──────────────────────────
    # Split muts into a named list: lineage name → dataframe of (pos, residue) rows
    lineage_muts <- split(
      muts[, c("pos", "residue"), drop = FALSE],
      muts$lineage
    )

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
        warning(
          "Sequence '",
          seq_names[seq_idx],
          "' matched lineages with equal specificity: ",
          paste(sort(candidates), collapse = ", "),
          ". Assigning first alphabetically."
        )
        candidates <- sort(candidates)[1L]
      }
      candidates
    }
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
