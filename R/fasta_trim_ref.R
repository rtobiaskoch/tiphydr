# Internal helper: coerce a DNAbin object to a named character matrix.
# ape::as.character() returns a list when given a single-sequence DNAbin, and a
# matrix when given multiple sequences. This wrapper normalises both cases to a
# 1-row (or N-row) character matrix with the sequence names as rownames, and
# all characters uppercased.
.dnabin_to_matrix <- function(dnabin) {
  raw_chars <- as.character(dnabin)

  # ape::as.character() on a DNAbin list (unaligned or short sequences) returns
  # a named list where each element is a character vector of individual bases.
  # Do NOT apply toupper() before rbind — toupper() on a list collapses it into
  # a flat character vector, losing names and structure.
  if (is.list(raw_chars)) {
    seq_names  <- names(raw_chars)
    seq_matrix <- do.call(rbind, raw_chars)
    rownames(seq_matrix) <- seq_names
    return(toupper(seq_matrix))
  }

  # Matrix case (aligned sequences with gap characters): toupper() works fine.
  toupper(raw_chars)
}

#' Align sequences to a reference and trim to reference coordinates
#'
#' Uses MAFFT for alignment. When fasta is the first alignment, all sequences
#' plus the reference are aligned together with ips::mafft(). When an existing
#' alignment is provided, new sequences are inserted via MAFFT --add profile
#' mode, which preserves existing alignment columns exactly.
#'
#' Coordinate invariant: after trimming, alignment position N equals ungapped
#' reference position N. This is the contract required by define_lineage().
#'
#' @param fasta DNAStringSet; sequences to align and trim
#' @param ref length-1 DNAStringSet; reference sequence defining trim coordinates
#' @param alignment optional DNAStringSet; existing alignment to add sequences
#'   into via MAFFT --add profile mode (NULL = build from scratch)
#' @param drop logical; if TRUE and alignment provided, remove sequences present
#'   in alignment but absent from fasta (default FALSE)
#' @param verbose logical; if TRUE prints progress messages (default TRUE)
#'
#' @return DNAStringSet of trimmed sequences (reference excluded)
#' @export
fasta_trim_ref <- function(fasta, ref, alignment = NULL, drop = FALSE,
                           verbose = TRUE) {

  if (length(ref) != 1L) stop("ref must be a single-sequence DNAStringSet (length 1)")

  if (nchar(Sys.which("mafft")) == 0L) {
    stop(
      "MAFFT not found on PATH.\n",
      "Install from: https://mafft.cbrc.jp/alignment/software/\n",
      "Alternative for basic alignment (no --add): BiocManager::install('msa')\n",
      "Profile alignment via --add requires MAFFT."
    )
  }

  ref_name <- names(ref)

  if (verbose) {
    message("Input: ", length(fasta), " sequences, avg length ",
            round(mean(Biostrings::width(fasta)), 1), " nt")
    message("Reference: ", ref_name, " (", Biostrings::width(ref), " nt)")
  }

  # ── Build alignment ──────────────────────────────────────────────────────────
  if (is.null(alignment)) {
    # Full alignment from scratch
    combined     <- c(ref, fasta)
    combined_bin <- ape::as.DNAbin(combined)
    aligned_bin  <- ips::mafft(combined_bin)
    aligned_mat  <- .dnabin_to_matrix(aligned_bin)

  } else {
    # Profile alignment: add new sequences via MAFFT --add
    existing_names <- names(alignment)
    new_seqs       <- fasta[!names(fasta) %in% existing_names]
    n_new          <- length(new_seqs)

    if (verbose) message("Adding ", n_new, " new sequences to existing alignment")

    if (n_new > 0L) {
      tmp_new   <- tempfile(fileext = ".fasta")
      tmp_exist <- tempfile(fileext = ".fasta")
      tmp_out   <- tempfile(fileext = ".fasta")
      on.exit(unlink(c(tmp_new, tmp_exist, tmp_out)), add = TRUE)

      # Include ref in the existing alignment so we can trim from it
      aln_with_ref <- if (ref_name %in% existing_names) {
        alignment
      } else {
        c(ref, alignment)
      }

      Biostrings::writeXStringSet(new_seqs,     tmp_new)
      Biostrings::writeXStringSet(aln_with_ref, tmp_exist)

      exit_code <- system2(
        "mafft",
        args   = c("--quiet", "--add", tmp_new, tmp_exist),
        stdout = tmp_out
      )
      if (exit_code != 0L) stop("MAFFT --add alignment failed (exit code ", exit_code, ")")

      added_aln   <- Biostrings::readDNAStringSet(tmp_out)
      aligned_mat <- .dnabin_to_matrix(ape::as.DNAbin(added_aln))

    } else {
      # No new sequences — use existing alignment directly (ensure ref is present)
      aln_with_ref <- if (ref_name %in% existing_names) {
        alignment
      } else {
        c(ref, alignment)
      }
      aligned_mat <- .dnabin_to_matrix(ape::as.DNAbin(aln_with_ref))
    }
  }

  if (verbose) message("Alignment length: ", ncol(aligned_mat), " nt")

  # ── Trim to reference coordinates ───────────────────────────────────────────
  # The coordinate invariant: only keep columns where the reference has a real
  # nucleotide (A/T/C/G). Gap columns inserted to accommodate other sequences
  # are discarded, so output position N == ungapped reference position N.
  ref_row         <- aligned_mat[ref_name, ]
  valid_positions <- which(ref_row %in% c("A", "T", "C", "G"))

  if (verbose) {
    message("Trimming columns ", valid_positions[1], " to ",
            tail(valid_positions, 1), " (", length(valid_positions), " reference nt)")
  }

  aligned_trimmed <- aligned_mat[, valid_positions, drop = FALSE]

  # Remove reference row — caller only needs the query sequences
  aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) != ref_name, , drop = FALSE]

  # Optionally drop sequences not in fasta (used when working from an existing alignment)
  if (drop && !is.null(alignment)) {
    keep_names      <- names(fasta)
    n_dropped       <- sum(!rownames(aligned_trimmed) %in% keep_names)
    if (verbose) message("Dropping ", n_dropped, " sequences not in fasta")
    aligned_trimmed <- aligned_trimmed[rownames(aligned_trimmed) %in% keep_names, , drop = FALSE]
  }

  # Collapse each row character vector back to a single sequence string
  seqs       <- apply(aligned_trimmed, 1, paste0, collapse = "")
  trimmed_set <- Biostrings::DNAStringSet(seqs)
  names(trimmed_set) <- rownames(aligned_trimmed)

  # Validate the coordinate invariant — warn rather than error so partial
  # results are still returned (caller can inspect).
  expected_len <- Biostrings::width(ref)
  actual_lens  <- unique(Biostrings::width(trimmed_set))
  if (!all(actual_lens == expected_len)) {
    warning("Trimmed length (", paste(actual_lens, collapse = ","),
            " nt) does not match reference length (", expected_len, " nt)")
  }

  if (verbose) {
    message("Output: ", length(trimmed_set), " sequences, ",
            expected_len, " nt (= reference length)")
  }

  trimmed_set
}
