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
#' @param persistence_days Numeric scalar. Days a clade must span (intro date to
#'   last sample) to be labelled \code{"persistent"}. Default 365 (one year).
#' @param palette_persistence RColorBrewer palette for the persistence colour
#'   scale (line + last-sample point). Default \code{"Dark2"}.
#' @param palette_source RColorBrewer palette for the source-deme fill scale
#'   (introduction-date point). Default \code{"Set2"}.
#'
#' @return A \code{ggplot} object. Additional \pkg{ggplot2} layers can be
#'   added with \code{+}.
#' @export
plot_explode_tree <- function(introductions,
                              persistence_days    = 365,
                              palette_persistence = "Dark2",
                              palette_source      = "Set2") {

  # --- Validate required columns -------------------------------------------
  required_cols <- c(
    "intro_clade_id", "deme", "inferred_intro_date",
    "last_sample_date", "inferred_intro_source", "clade_size"
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
  # persistence: did the clade persist for longer than persistence_days?
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
      persistence = dplyr::if_else(
        as.numeric(.data$last_sample_date - .data$inferred_intro_date) > persistence_days,
        "persistent",
        "transient"
      )
    )

  # Long format: two rows per clade (one per date type).
  # Using bind_rows avoids a tidyr dependency.
  clade_long <- dplyr::bind_rows(
    dplyr::mutate(clade_tbl,
                  date      = .data$inferred_intro_date,
                  date_type = "inferred_intro_date"),
    dplyr::mutate(clade_tbl,
                  date      = .data$last_sample_date,
                  date_type = "last_sample_date")
  )

  # Separate layers: intro date (fill = source) vs last sample (color = persistence).
  # color and fill are independent scales in ggplot2 — no extra package needed.
  intro_pts  <- dplyr::filter(clade_long, .data$date_type == "inferred_intro_date")
  sample_pts <- dplyr::filter(clade_long, .data$date_type == "last_sample_date")

  ggplot2::ggplot(clade_long,
                  ggplot2::aes(y = .data$intro_clade_id)) +
    # Segment colored by persistence (was the clade transient or long-lived?)
    ggplot2::geom_line(
      ggplot2::aes(x     = .data$date,
                   group = .data$intro_clade_id,
                   color = .data$persistence),
      linewidth = 1,
      alpha     = 0.8
    ) +
    # Last-sample point: open circle, outline colored by persistence
    ggplot2::geom_point(
      data = sample_pts,
      ggplot2::aes(x     = .data$date,
                   color = .data$persistence,
                   size  = .data$clade_size),
      shape = 21,
      fill  = "white"
    ) +
    # Introduction-date point: filled diamond, fill colored by source deme
    ggplot2::geom_point(
      data = intro_pts,
      ggplot2::aes(x    = .data$date,
                   fill = .data$inferred_intro_source,
                   size = .data$clade_size),
      shape = 22,
      color = "white"
    ) +
    ggplot2::scale_color_brewer(palette = palette_persistence, name = "Persistence") +
    ggplot2::scale_fill_brewer(palette = palette_source,       name = "Source deme") +
    ggplot2::scale_size_continuous(range = c(2.5, 8),          guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(.data$deme), scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Introduction clade") +
    ggplot2::theme_classic()
}
