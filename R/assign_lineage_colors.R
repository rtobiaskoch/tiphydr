#' Assign colors to lineage names by family, specificity, and hybrid mixture
#'
#' Generates a color for every lineage name in `muts`, plus dedicated
#' `"unknown"`/`"ambiguous"` colors (the sentinel labels [define_lineage()]
#' produces when no or an ambiguous match is found). Colors follow the
#' lineage-naming convention used throughout this package (see
#' [define_lineage()]): a bare family name (e.g. `"NY10"`) gets its base
#' color at full saturation; a single-family partial (e.g. `"NY10_NS2A"`)
#' gets a lighter tint of that family's color, proportional to how many of
#' the family's required mutations it carries; a hybrid name joining two
#' pieces with `"-"` (e.g. `"NY10-SW03_NS5"`) gets a blend of its two
#' pieces' colors, weighted toward whichever piece requires more mutations.
#'
#' @param muts dataframe with columns `lineage` (chr), `pos` (int),
#'   `residue` (chr) — the same lookup table shape passed to
#'   [define_lineage()] (e.g. [wnv_mut_gff] or [wnv_mut_orf]). Row counts per
#'   lineage drive both the partial-lightening fraction and the hybrid blend
#'   weighting.
#' @param base_colors dataframe with columns `lineage` (chr, the family
#'   name) and `color` (chr, a hex color) giving each family's base color at
#'   full saturation. Defaults to [wnv_lineage_colors]. Every lineage in
#'   `muts` must resolve to exactly one `base_colors` family, either by exact
#'   match (a bare family name) or by `"<family>_"` prefix (a partial).
#' @param unknown_color,ambiguous_color length-1 character; hex colors for
#'   the `"unknown"` and `"ambiguous"` sentinel labels [define_lineage()]
#'   can produce. These are not mutation-based lineages (they have no rows
#'   in `muts`), so they're always appended to the output rather than driven
#'   by `muts` row counts.
#'
#' @return tibble with columns `lineage` (chr) and `color` (chr, hex).
#' @export
#' @importFrom purrr map_chr map_dbl map_lgl
#' @importFrom tibble tibble
assign_lineage_colors <- function(
  muts,
  base_colors = NULL,
  unknown_color = "#BBBBBB",
  ambiguous_color = "#888888"
) {

  # ── Input validation ─────────────────────────────────────────────────────────

  # base_colors defaults to the package's own wnv_lineage_colors data; NULL
  # is used as the formal default (rather than the data object directly) so
  # R CMD check doesn't flag it as an unresolved global variable.
  if (is.null(base_colors)) {
    base_colors <- wnv_lineage_colors
  }

  required_cols <- c("lineage", "pos", "residue")
  if (!all(required_cols %in% names(muts))) {
    stop("muts must have columns: ", paste(required_cols, collapse = ", "))
  }
  if (!all(c("lineage", "color") %in% names(base_colors))) {
    stop("base_colors must have columns: lineage, color")
  }

  # Row count per lineage: how many mutations each lineage requires. Drives
  # both the partial-lightening fraction and the hybrid blend weighting.
  muts_counts <- as.data.frame(table(lineage = muts$lineage), stringsAsFactors = FALSE)
  names(muts_counts) <- c("lineage", "n")

  # ── Per-lineage color resolution ─────────────────────────────────────────────

  resolve_color <- function(lineage_name) {

    # Exact match to a known base family: full saturation, no lightening.
    exact_hit <- base_colors$color[base_colors$lineage == lineage_name]
    if (length(exact_hit) == 1L) {
      return(exact_hit)
    }

    # Hybrid: two pieces joined by "-", each itself a base family or a
    # single-family partial. Resolve each piece's color and mutation-count
    # weight independently, then blend toward the more specific (higher
    # mutation count) piece.
    if (grepl("-", lineage_name, fixed = TRUE)) {
      pieces <- strsplit(lineage_name, "-", fixed = TRUE)[[1]]
      if (length(pieces) != 2L) {
        stop(
          "Cannot resolve hybrid lineage name with more than two '-'-joined ",
          "pieces: '", lineage_name, "'"
        )
      }
      piece_colors  <- purrr::map_chr(pieces, resolve_color)
      piece_weights <- purrr::map_dbl(pieces, function(p) muts_counts$n[muts_counts$lineage == p])

      # colorspace::mixcolor(alpha, color1, color2): alpha = 0 -> pure color1,
      # alpha = 1 -> pure color2. To favor whichever piece has MORE required
      # mutations, alpha must be the OTHER piece's weight fraction.
      alpha <- piece_weights[2] / sum(piece_weights)
      mixed <- colorspace::mixcolor(
        alpha,
        colorspace::hex2RGB(piece_colors[1]),
        colorspace::hex2RGB(piece_colors[2])
      )
      return(colorspace::hex(mixed))
    }

    # Single-family partial (e.g. "NY10_NS2A"): find the base family it
    # belongs to (the family name that is a "<family>_" prefix of this
    # lineage), then lighten that family's color in proportion to how much
    # of the family's full mutation set is present.
    family <- base_colors$lineage[
      purrr::map_lgl(base_colors$lineage, function(fam) startsWith(lineage_name, paste0(fam, "_")))
    ]
    if (length(family) != 1L) {
      stop(
        "Could not resolve a unique base family for lineage '", lineage_name,
        "'. Check that base_colors contains exactly one matching '<family>_' prefix."
      )
    }

    family_color <- base_colors$color[base_colors$lineage == family]
    family_n     <- muts_counts$n[muts_counts$lineage == family]
    own_n        <- muts_counts$n[muts_counts$lineage == lineage_name]
    fraction     <- own_n / family_n

    colorspace::lighten(family_color, amount = 1 - fraction)
  }

  # ── Assemble output ──────────────────────────────────────────────────────────
  # One row per lineage present in muts, plus the two sentinel labels
  # define_lineage() can produce even though they have no rows in muts.

  mut_lineages <- unique(muts$lineage)
  result <- tibble::tibble(
    lineage = c(mut_lineages, "unknown", "ambiguous"),
    color = c(
      purrr::map_chr(mut_lineages, resolve_color),
      unknown_color,
      ambiguous_color
    )
  )

  result
}

# wnv_lineage_colors is package data (data/wnv_lineage_colors.rda), lazy-loaded
# at runtime but invisible to R CMD check's static analysis of the function
# body above; this silences the resulting "no visible binding" NOTE.
utils::globalVariables("wnv_lineage_colors")
