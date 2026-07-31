#' Read calendar dates off every node, as computed by the dating tool itself
#'
#' Dating tools like LSD2 already compute a calendar date for every node
#' (tip and internal) and embed it as a NHX-style `[&date=...]` comment in
#' their native NEXUS output. This reads those dates directly rather than
#' re-deriving anything from branch lengths -- no regression, no model, no
#' assumption about the tree's time unit. Because the dating tool enforces
#' chronological node ordering internally (an ancestor is never dated later
#' than its descendants), reading its own dates guarantees the same holds
#' here, regardless of clock signal strength or rate heterogeneity.
#'
#' `annotated_tree_path` must point to the dating tool's native output --
#' e.g. `{prefix}.timetree.nex`, before any later step (such as a
#' zero-branch adjustment) round-trips the tree through `ape::read.nexus()`
#' / `ape::write.nexus()`, which silently drops NHX comments. `tree` is
#' whatever downstream working tree (e.g. after that adjustment, or with
#' tips relabeled) the dates need to be matched onto; it must share the
#' exact same topology and tip order as the annotated tree, since node
#' identity is matched by ape node id.
#'
#' @param tree An `ape::phylo` tree to attach dates to. Must have the same
#'   topology (tip order, edge structure) as the tree in
#'   `annotated_tree_path`.
#' @param annotated_tree_path Path to a NEXUS/BEAST-style tree file (LSD2's
#'   native output) with a `date` annotation on every node, readable by
#'   `treeio::read.beast()`.
#'
#' @return A tibble with columns `node` (ape node id) and `inferred_date`
#'   (`Date`) for every tip and internal node of `tree`.
#' @export
node_dates_from_timetree <- function(tree, annotated_tree_path) {
  annotated <- treeio::read.beast(annotated_tree_path)
  atree <- annotated@phylo

  if (
    !identical(ape::Ntip(atree), ape::Ntip(tree)) ||
      !identical(atree$edge[, 1:2], tree$edge[, 1:2])
  ) {
    stop(
      "tree and annotated_tree_path do not share the same topology -- ",
      "node ids would not refer to the same nodes"
    )
  }

  node_dates <- tibble::as_tibble(annotated@data) |>
    dplyr::transmute(
      node = as.integer(.data$node),
      inferred_date = parse_annotated_date(.data$date)
    )

  all_nodes <- tibble::tibble(node = seq_len(ape::Ntip(tree) + ape::Nnode(tree)))
  result <- dplyr::left_join(all_nodes, node_dates, by = "node")

  if (anyNA(result$inferred_date)) {
    missing <- result$node[is.na(result$inferred_date)]
    stop(
      "no annotated date found for node(s): ",
      paste(missing, collapse = ", ")
    )
  }

  result
}

#' Parse a dating tool's native per-node date annotation
#'
#' Handles the three date shapes actually observed across this pipeline's
#' two dating engines: LSD2/IQ-TREE3 writes a bare 4-digit year when the
#' fitted date has no finer resolution, or a full `YYYY-MM-DD` otherwise.
#' TreeTime (after `treetime_to_timetree.R`'s NEXUS normalization) writes a
#' decimal year (e.g. `"2019.40"`) -- normalizing the annotation *shape* to
#' LSD2's does not convert the date *value* itself, so this still has to be
#' handled separately here.
#'
#' @param x Character vector of raw date strings from the annotated tree.
#' @return A `Date` vector, `NA` for any string matching none of the three
#'   known shapes.
#' @noRd
parse_annotated_date <- function(x) {
  is_bare_year <- grepl("^[0-9]{4}$", x)
  is_iso_date <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  is_decimal_year <- !is_bare_year & !is_iso_date & grepl("^[0-9]{4}[.][0-9]+$", x)

  # Assign by subset (not dplyr::case_when()) for the same reason
  # extract_uncertain_dates.R does: case_when() evaluates every branch's RHS
  # across the whole vector before masking, so as.Date() on a decimal-year
  # string (or decimal_year_to_date() on a bare-year/ISO string) would error
  # even on rows that branch never actually applies to. Each branch is also
  # guarded by any(): paste0(character(0), "-07-02") does NOT return
  # character(0) as a naive reading would suggest -- paste0() treats a
  # zero-length argument as contributing nothing rather than collapsing the
  # whole result to zero length, so it silently returns "-07-02" (length 1)
  # when is_bare_year is all-FALSE, which as.Date() then fails to parse.
  out <- as.Date(rep(NA_character_, length(x)))
  # LSD2 emits a bare "YYYY" when the fitted date has no finer resolution;
  # treat it as the midpoint of that year.
  if (any(is_bare_year)) {
    out[is_bare_year] <- as.Date(paste0(x[is_bare_year], "-07-02"))
  }
  if (any(is_iso_date)) {
    out[is_iso_date] <- as.Date(x[is_iso_date])
  }
  if (any(is_decimal_year)) {
    out[is_decimal_year] <- decimal_year_to_date(x[is_decimal_year])
  }
  out
}

#' Convert a TreeTime decimal year (e.g. "2021.31") to a calendar `Date`
#'
#' Leap-year aware: uses the actual number of days in the specific calendar
#' year the decimal falls in, not a fixed 365/365.25. Mirrors
#' `wnv-ge_dta_stan`'s `scripts/extract_uncertain_dates.R::decimal_year_to_date()`
#' exactly, so both pipelines convert the same TreeTime output to the same
#' calendar date.
#'
#' @param x Character vector of decimal-year strings.
#' @return A `Date` vector.
#' @noRd
decimal_year_to_date <- function(x) {
  dy <- as.numeric(x)
  year <- floor(dy)
  start_of_year <- as.Date(sprintf("%d-01-01", year))
  start_of_next <- as.Date(sprintf("%d-01-01", year + 1L))
  days_in_year <- as.numeric(start_of_next - start_of_year)
  start_of_year + round((dy - year) * days_in_year)
}
