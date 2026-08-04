#' Filter, order, and number clades within each deme for plotting
#'
#' Shared row-ordering logic between plot_explode_tree() and
#' plot_clade_persistence() -- keeping both in one place guarantees their
#' rows line up when the two plots are viewed side by side.
#'
#' @param clade_tbl A tibble with (at least) clade_size, deme,
#'   inferred_intro_date columns, one row per clade.
#' @param min_clade Minimum clade_size to keep.
#' @return clade_tbl, filtered, arranged by (deme, inferred_intro_date), with
#'   a local_clade_num column added (1..N within each deme).
.order_clades <- function(clade_tbl, min_clade = 1) {
  clade_tbl |>
    dplyr::filter(.data$clade_size >= min_clade) |>
    dplyr::arrange(.data$deme, .data$inferred_intro_date) |>
    dplyr::mutate(
      local_clade_num = dplyr::row_number(),
      .by = "deme"
    )
}

#' Plot introduction clade timelines from explode_tree output
#'
#' Visualises each introduction clade as a horizontal segment: from
#' \code{inferred_intro_date} to \code{last_sample_date}, with the
#' introduction point (filled diamond) coloured by \code{inferred_intro_source}
#' and the last-sample point (open circle) and segment coloured by persistence
#' (\code{last_sample_date - inferred_intro_date > persistence_days}).
#' One facet per destination deme.
#'
#' @param introductions A tibble from \code{\link{explode_tree}()$introductions}.
#'   Required columns: \code{intro_clade_id}, \code{deme},
#'   \code{inferred_intro_date}, \code{last_sample_date},
#'   \code{inferred_intro_source}, \code{clade_size}.
#' @param persistence_days Numeric scalar. Days below which a clade is labelled
#'   \code{"transient"}. Default 365 (one year). When \code{persistent_days} is
#'   also supplied, clades between the two thresholds are labelled
#'   \code{"semi-persistent"}.
#' @param persistent_days Numeric scalar or \code{NULL}. Days above which a clade
#'   is labelled \code{"persistent"}. If \code{NULL} (default), the old two-class
#'   behaviour is used (\code{"persistent"} / \code{"transient"} only).
#' @param palette_source RColorBrewer palette for the source-deme fill scale
#'   (introduction-date point). Default \code{"Set2"}.
#'
#' @return A \code{ggplot} object. Additional \pkg{ggplot2} layers can be
#'   added with \code{+}.
#' @export
plot_explode_tree <- function(
  introductions,
  persistence_days = 365,
  persistent_days  = NULL,
  min_clade = 1,
  palette_source = "Set2",
  named_palette = NULL
) {
  # --- Validate required columns -------------------------------------------
  required_cols <- c(
    "intro_clade_id",
    "deme",
    "inferred_intro_date",
    "last_sample_date",
    "inferred_intro_source",
    "clade_size"
  )
  missing_cols <- setdiff(required_cols, names(introductions))
  if (length(missing_cols) > 0L) {
    stop(
      "introductions is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ".\nPass the $introductions element from explode_tree()."
    )
  }

  # Introductions tibble has one row per tip; collapse to one row per clade.
  # Three-class when persistent_days is supplied; two-class otherwise (backward compat).
  clade_tbl <- introductions |>
    dplyr::distinct(
      .data$intro_clade_id,
      .data$deme,
      .data$inferred_intro_date,
      .data$last_sample_date,
      .data$inferred_intro_source,
      .data$clade_size
    ) |>
    dplyr::mutate(
      .persist_num = as.numeric(.data$last_sample_date - .data$inferred_intro_date),
      persistence = if (is.null(persistent_days)) {
        dplyr::if_else(.data$.persist_num > persistence_days, "persistent", "transient")
      } else {
        dplyr::case_when(
          .data$.persist_num <  persistence_days ~ "transient",
          .data$.persist_num >  persistent_days  ~ "persistent",
          .default                               = "semi-persistent"
        )
      }
    ) |>
    .order_clades(min_clade)

  # Guard: RColorBrewer palettes have a fixed maximum. Skip when named_palette
  # is supplied (no maxcolors constraint for explicit named vectors).
  if (is.null(named_palette)) {
    n_sources <- length(unique(clade_tbl$inferred_intro_source))
    max_source_colors <- RColorBrewer::brewer.pal.info[
      palette_source,
      "maxcolors"
    ]
    if (n_sources > max_source_colors) {
      stop(
        "Number of source demes (",
        n_sources,
        ") exceeds the maximum colours in ",
        "palette '",
        palette_source,
        "' (",
        max_source_colors,
        "). ",
        "Supply a palette with more colours via palette_source or named_palette."
      )
    }
  }

  # Long format: two rows per clade (one per date type).
  # Using bind_rows avoids a tidyr dependency.
  clade_long <- dplyr::bind_rows(
    dplyr::mutate(
      clade_tbl,
      date = .data$inferred_intro_date,
      date_type = "inferred_intro_date"
    ),
    dplyr::mutate(
      clade_tbl,
      date = .data$last_sample_date,
      date_type = "last_sample_date"
    )
  )

  # Separate layers: intro date (fill = source) vs last sample (color = persistence).
  # color and fill are independent scales in ggplot2 — no extra package needed.
  intro_pts <- dplyr::filter(
    clade_long,
    .data$date_type == "inferred_intro_date"
  )
  sample_pts <- dplyr::filter(clade_long, .data$date_type == "last_sample_date")

  ggplot2::ggplot(clade_long, ggplot2::aes(y = .data$local_clade_num)) +
    # Segment colored by persistence (was the clade transient or long-lived?)
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data$date,
        group = .data$intro_clade_id,
        color = .data$persistence
      ),
      linewidth = 1,
      alpha = 0.8
    ) +
    # Last-sample point: open circle, outline colored by persistence
    ggplot2::geom_point(
      data = sample_pts,
      ggplot2::aes(
        x = .data$date,
        color = .data$persistence,
        size = .data$clade_size
      ),
      shape = 21,
      fill = "white"
    ) +
    # Introduction-date point: filled diamond, fill colored by source deme
    ggplot2::geom_point(
      data = intro_pts,
      ggplot2::aes(
        x = .data$date,
        fill = .data$inferred_intro_source,
        size = .data$clade_size
      ),
      shape = 22,
      color = "white"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "persistent"      = "grey25",
        "semi-persistent" = "grey45",
        "transient"       = "grey65"
      ),
      name = "Persistence"
    ) +
    (if (!is.null(named_palette)) {
      ggplot2::scale_fill_manual(values = named_palette, name = "Source deme")
    } else {
      ggplot2::scale_fill_brewer(palette = palette_source, name = "Source deme")
    }) +
    ggplot2::scale_size_continuous(range = c(2.5, 8), name = "Clade size") +
    ggplot2::scale_y_continuous(breaks = function(x) {
      seq(ceiling(x[1]), floor(x[2]))
    }) +
    ggplot2::facet_wrap(ggplot2::vars(.data$deme), scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Introduction") +
    ggplot2::theme_classic()
}
