#' Plot per-clade persistence_days distributions from simmap replicates
#'
#' Companion to plot_explode_tree(): one horizontal box-and-whisker per
#' introduction clade, showing the distribution of persistence_days across
#' the simmap posterior (the n_sim make.simmap() replicates). persistence_days
#' is a duration, not a calendar date, so it can't share plot_explode_tree()'s
#' x-axis -- this is a separate plot, not a layer added to that one. Row
#' ordering matches plot_explode_tree() exactly (via the shared
#' .order_clades() helper) so the two plots align when viewed side by side.
#'
#' @param raw_clade_dwell Tibble, one row per (sim_id, intro_node): intro_node,
#'   deme, persistence_days (the shape of dta_compute.R's
#'   _dta_clade_dwell_raw.tsv). Values come from here only -- never from
#'   clade_dwell -- so the plot reflects the actual posterior, not a summary
#'   of a summary.
#' @param clade_dwell Tibble, one row per clade: intro_clade_id, intro_node,
#'   deme, inferred_intro_date, clade_size (the shape of
#'   _dta_clade_dwell.tsv). Used only for row ordering and the min_clade
#'   filter.
#' @param min_clade Minimum clade_size to plot. Default 1.
#'
#' @return A ggplot object.
#' @export
plot_clade_persistence <- function(raw_clade_dwell, clade_dwell, min_clade = 1) {
  required_raw <- c("intro_node", "deme", "persistence_days")
  missing_raw <- setdiff(required_raw, names(raw_clade_dwell))
  if (length(missing_raw) > 0L) {
    stop(
      "raw_clade_dwell is missing required columns: ",
      paste(missing_raw, collapse = ", ")
    )
  }

  required_clade <- c("intro_clade_id", "intro_node", "deme", "inferred_intro_date", "clade_size")
  missing_clade <- setdiff(required_clade, names(clade_dwell))
  if (length(missing_clade) > 0L) {
    stop(
      "clade_dwell is missing required columns: ",
      paste(missing_clade, collapse = ", ")
    )
  }

  ordered <- .order_clades(clade_dwell, min_clade)

  plot_df <- raw_clade_dwell |>
    dplyr::inner_join(
      dplyr::select(ordered, "intro_node", "deme", "intro_clade_id", "local_clade_num"),
      by = c("intro_node", "deme")
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$persistence_days,
      y = .data$local_clade_num,
      group = .data$intro_clade_id
    )
  ) +
    ggplot2::geom_boxplot(
      orientation = "y",
      fill = "grey70",
      color = "grey40"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$deme), scales = "free_y") +
    ggplot2::scale_y_continuous(breaks = function(x) {
      seq(ceiling(x[1]), floor(x[2]))
    }) +
    ggplot2::labs(x = "Persistence (days)", y = "Introduction") +
    ggplot2::theme_classic()
}
