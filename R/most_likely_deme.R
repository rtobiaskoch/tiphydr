#' Reduce a node x deme probability matrix to the most likely deme per node
#'
#' Discrete trait ancestral reconstruction (e.g. TreeTime DTA / mugration) returns
#' a probability for every deme at every node. This collapses that matrix to the
#' single most probable deme per node and its probability, the form needed to
#' detect state transitions along the tree.
#'
#' Input rows may be duplicated (some tools emit tip rows twice); identical rows
#' are de-duplicated first. Ties are broken by column order (first deme wins).
#'
#' @param node_probs A data frame / tibble with an integer \code{node} column and
#'   one numeric column per deme holding probabilities. Column names are the deme
#'   names.
#'
#' @return A tibble with columns \code{node}, \code{inferred_state} (deme name),
#'   and \code{confidence_state} (the probability of that deme).
#' @export
most_likely_deme <- function(node_probs) {

  # ── Validation ───────────────────────────────────────────────────────────────
  np <- tibble::as_tibble(node_probs)

  if (!"node" %in% names(np)) {
    stop("node_probs must contain a 'node' column")
  }

  deme_cols <- setdiff(names(np), "node")
  if (length(deme_cols) == 0L) {
    stop("node_probs must contain at least one deme (probability) column")
  }

  # ── Core ─────────────────────────────────────────────────────────────────────

  # Drop redundant duplicate rows (e.g. tip rows emitted twice by TreeTime).
  np <- dplyr::distinct(np)

  # Probability matrix; max.col is a vectorised "which is the max per row".
  mat <- as.matrix(np[deme_cols])
  best_idx <- max.col(mat, ties.method = "first")

  tibble::tibble(
    node             = np$node,
    inferred_state   = deme_cols[best_idx],
    # value at [row, winning column] via matrix indexing
    confidence_state = mat[cbind(seq_len(nrow(mat)), best_idx)]
  )
}
