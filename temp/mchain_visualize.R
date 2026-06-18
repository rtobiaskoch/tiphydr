#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  N E   --------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot Ne trajectories over time, one line per chain, faceted and colored by deme.
# df : long-format summary dataframe (output of mchain_clean / mchain_ne_transform)
#      required columns: chain, deme, month, ne_median, ne_lower, ne_upper
plot_ne_trajectories = function(df, chains = "all", deme_colors = NULL) {

  # Filter to requested chains before plotting
  if (!identical(chains, "all")) {
    df = dplyr::filter(df, chain %in% chains)
  }

  # Detect the time column dynamically — whatever remains after removing the
  # fixed columns. Handles both "month" (default config) and "year" (R019.3).
  fixed_cols <- c("chain", "deme", "ne_median", "ne_lower", "ne_upper")
  time_col   <- setdiff(names(df), fixed_cols)[[1]]

  # tidyr::extract always produces character columns; coerce to integer so
  # scale_x_continuous works and ordering is numeric not alphabetical.
  df <- dplyr::mutate(df, !!time_col := as.integer(.data[[time_col]]))

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x      = .data[[time_col]],
      y      = ne_median,
      color  = deme,
      fill   = deme,
      # group ensures a separate line per chain within each deme panel
      group  = interaction(chain, deme)
    )
  ) +
    # 95% credibility ribbon per chain
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ne_lower, ymax = ne_upper),
      alpha = 0.15,
      color = NA
    ) +
    # Median trajectory line per chain
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ deme) +
    # Show every other tick to reduce label crowding
    ggplot2::scale_x_continuous(
      breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2)
    ) +
    ggplot2::labs(
      x     = time_col,
      y     = "Effective Population Size (Ne)",
      color = "Deme",
      fill  = "Deme"
    ) +
    ggplot2::scale_color_manual(values = deme_colors, na.value = "grey70") +
    ggplot2::scale_fill_manual(values = deme_colors, na.value = "grey70") +
    ggplot2::theme_classic()
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  M I G  --------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot a chord diagram of migration rates between demes.
# Aggregates across chains by taking the median of per-chain medians.
# Uses circlize::chordDiagram (base R graphics, not ggplot2).
#
# df : long-format summary dataframe (output of mchain_mig_clean)
#      required columns: chain, from, to, mig_median
plot_mig = function(df, chains = "all", deme_colors = NULL) {

  # Filter to requested chains before aggregating across chains
  if (!identical(chains, "all")) {
    df = dplyr::filter(df, chain %in% chains)
  }

  # If time_index is present (multi-rate model), preserve it through aggregation
  # so each rate period gets its own chord diagram panel.
  has_time_index = "time_index" %in% names(df)
  group_vars     = if (has_time_index) c("from", "to", "time_index") else c("from", "to")

  # Aggregate across chains: median of per-chain median migration rates
  mig_agg = df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(value = median(mig_median), .groups = "drop")

  # Use caller-supplied palette so chord sector colors match all other figures.
  # Fall back to building from the demes present in this plot if not supplied.
  if (is.null(deme_colors)) {
    demes       = sort(unique(c(mig_agg$from, mig_agg$to)))
    deme_colors = build_deme_palette(demes)
  }

  rate_periods = if (has_time_index) sort(unique(mig_agg$time_index)) else NA_character_

  # Side-by-side panels when multiple rate periods exist
  if (has_time_index) graphics::par(mfrow = c(1L, length(rate_periods)))

  purrr::walk(rate_periods, function(period) {
    df_period = if (has_time_index) {
      mig_agg |>
        dplyr::filter(time_index == period) |>
        dplyr::select(from, to, value)   # chordDiagram expects these three columns
    } else {
      dplyr::select(mig_agg, from, to, value)
    }

    circlize::circos.clear()
    # Expand canvas so labels drawn outside the ring don't clip at [-1,1].
    circlize::circos.par(canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5))
    circlize::chordDiagram(
      df_period,
      grid.col          = deme_colors,
      transparency      = 0.4,
      directional       = 1,
      direction.type    = "arrows",
      link.arr.type     = "big.arrow",
      annotationTrack   = "grid",
      preAllocateTracks = list(track.height = 0.05)
    )

    circlize::circos.trackPlotRegion(
      track.index = 1,
      panel.fun = function(x, y) {
        sector_name = circlize::get.cell.meta.data("sector.index")
        # Place text at the outer edge of the label track plus a small gap,
        # so labels sit clear of the chord ring.
        y_out = circlize::CELL_META$ylim[2] + 0.5
        circlize::circos.text(
          circlize::CELL_META$xcenter,
          y_out,
          sector_name,
          facing     = "bending.outside",
          niceFacing = TRUE,
          adj        = c(0.5, 0),
          cex        = 0.9,
          col        = unname(deme_colors[sector_name]),
          font       = 1,            # plain weight — closest to light in base R
          family     = "Helvetica Neue"
        )
      },
      bg.border = NA
    )

    if (has_time_index) graphics::title(paste("Rate period", period), line = -1)
    circlize::circos.clear()
  })

  if (has_time_index) graphics::par(mfrow = c(1L, 1L))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  M I G  H E A T M A P -----------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot a tiled heatmap of median migration rates, faceted by chain.
# Each tile represents the from -> to migration rate for one chain,
# allowing visual comparison of migration patterns across chains.
#
# df     : long-format summary dataframe (output of mchain_mig_clean)
#          required columns: chain, from, to, mig_median
# chains : character vector of chain names to include, or "all"
plot_mig_heat = function(df, chains = "all", deme_colors = NULL) {

  # Filter to requested chains before plotting
  if (!identical(chains, "all")) {
    df = dplyr::filter(df, chain %in% chains)
  }

  # When time_index is present (multi-rate model), facet_grid puts rate periods
  # on rows and chains on columns so both dimensions are visible at once.
  # When absent (single-rate), fall back to facet_wrap by chain only.
  has_time_index = "time_index" %in% names(df)

  facet_layer = if (has_time_index) {
    ggplot2::facet_grid(
      rows = ggplot2::vars(time_index),
      cols = ggplot2::vars(chain),
      labeller = ggplot2::labeller(
        time_index = function(x) paste("Rate period", x)
      )
    )
  } else {
    ggplot2::facet_wrap(~ chain)
  }

  df |>
    ggplot2::ggplot(ggplot2::aes(x = to, y = from, fill = mig_median)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(option = "plasma", name = "Median\nMig. Rate") +
    facet_layer +
    ggplot2::labs(x = "To", y = "From") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      panel.spacing = ggplot2::unit(1, "lines")
    )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  T R A C E --------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot MCMC trace lines for one or more parameters to diagnose chain convergence,
# paired with a bar plot of ESS (Effective Sample Size) per chain per parameter.
#
# Trace plot: a well-mixed chain looks like a "hairy caterpillar" — stable,
# overlapping lines. Trending or drifting traces indicate lack of convergence.
#
# ESS bar plot: coda::effectiveSize() uses spectral density estimation — the
# same estimator as Tracer. The dashed red line marks ESS = 200, Tracer's
# minimum threshold for acceptable mixing.
#
# df     : combined raw chain dataframe (output of mchain_combine)
#          required columns: Sample, chain, and any requested params
# params : character vector of column names to plot
#          defaults to overall posterior and tree height as standard diagnostics
# chains : character vector of chain names to include, or "all"
plot_trace = function(df, params = c("posterior", "Tree.height"), chains = "all") {
  
  # Guard: check structural columns are present
  if (!all(c("Sample", "chain") %in% names(df))) {
    stop("dataframe must contain 'Sample' and 'chain' columns (output of mchain_combine).")
  }
  
  # Filter to requested chains before building color palette and plotting
  if (!identical(chains, "all")) {
    df = dplyr::filter(df, chain %in% chains)
  }
  
  # Resolve param names case-insensitively against actual column names so that
  # e.g. "Tree.height" finds "tree.height" and vice versa.
  actual_cols <- names(df)
  params <- purrr::map_chr(params, function(p) {
    hit <- actual_cols[tolower(actual_cols) == tolower(p)]
    if (length(hit) > 0) hit[1] else p
  })

  # Warn if any requested params are still missing after case-insensitive resolution
  missing_params = setdiff(params, actual_cols)
  if (length(missing_params) > 0) {
    warning("params not found in df and will be skipped: ",
            paste(missing_params, collapse = ", "))
  }

  # Build chain color palette from WesAnderson "Zissou1"
  # Expand with colorRampPalette if chain count exceeds palette size
  chain_ids    = unique(df$chain)
  pal_raw      = wesanderson::wes_palette("Zissou1", n = 5)  # 5-color base palette
  pal          = grDevices::colorRampPalette(pal_raw)(length(chain_ids))
  chain_colors = stats::setNames(pal, chain_ids)
  
  # Pivot to long format once — reused by both the trace plot and ESS calculation
  df_long = df |>
    dplyr::select(Sample, chain, dplyr::any_of(params)) |>
    tidyr::pivot_longer(
      cols      = -c(Sample, chain),
      names_to  = "parameter",
      values_to = "value"
    )
  
  # --- Trace plot -----------------------------------------------------------
  p_trace = df_long |>
    ggplot2::ggplot(ggplot2::aes(x = Sample, y = value, color = chain)) +
    # One thin, semi-transparent line per chain — overlap reveals mixing
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.3) +
    # Free y scale: parameters are on very different scales
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(values = chain_colors) +
    ggplot2::labs(
      x     = "MCMC Sample",
      y     = "Value",
      color = "Chain"
    ) +
    ggplot2::theme_classic()
  
  # --- ESS calculation ------------------------------------------------------
  # coda::effectiveSize() returns a named numeric of length 1 per vector;
  # as.numeric() strips the name so dplyr::summarise doesn't complain
  ess_df = df_long |>
    dplyr::group_by(chain, parameter) |>
    dplyr::summarise(
      ess = as.numeric(coda::effectiveSize(value)),
      .groups = "drop"
    )
  
  # --- ESS bar plot ---------------------------------------------------------
  p_ess = ess_df |>
    ggplot2::ggplot(ggplot2::aes(x = chain, y = ess, fill = chain)) +
    ggplot2::geom_col() +
    # Tracer's minimum acceptable ESS threshold
    ggplot2::geom_hline(
      yintercept = 200,
      linetype   = "dashed",
      color      = "red",
      linewidth  = 0.5
    ) +
    # Free y: ESS scales differ across parameters
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = chain_colors) +
    ggplot2::labs(
      x    = "Chain",
      y    = "ESS",
      fill = "Chain"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # --- Combine side by side -------------------------------------------------
  # wrap_plots() is the namespace-safe patchwork alternative to `|`
  patchwork::wrap_plots(p_trace, p_ess, ncol = 2, widths = c(3, 1))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  G L M  -------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot violin distributions of MascotGLM covariate scaler posteriors.
# Ne GLM and migration GLM are shown as separate panels stacked vertically.
#
# Scalers are log-scale effect sizes from the GLM prior: positive = covariate
# increases the rate, negative = decreases. The dashed line at y = 0 marks
# no effect. Columns selected by prefix: NeGLM.scaler.* and migrationGLM.scaler.*
#
# df     : combined raw chain dataframe (output of mchain_combine)
#          must contain Sample, chain, and any NeGLM/migrationGLM columns
# chains : character vector of chain names to include, or "all"
plot_mascot_glm = function(df, chains = "all") {

  # Guard: check structural columns are present
  if (!all(c("Sample", "chain") %in% names(df))) {
    stop("df must contain 'Sample' and 'chain' columns (output of mchain_combine).")
  }

  # Filter to requested chains before plotting
  if (!identical(chains, "all")) {
    df = dplyr::filter(df, chain %in% chains)
  }

  # Build Zissou1 palette keyed by chain — same as plot_trace
  chain_ids    = unique(df$chain)
  pal_raw      = wesanderson::wes_palette("Zissou1", n = 5)
  pal          = grDevices::colorRampPalette(pal_raw)(length(chain_ids))
  chain_colors = stats::setNames(pal, chain_ids)

  # --- Reshape GLM scaler columns to long format ----------------------------
  # Selects all columns whose names contain "migrationGLM.scaler" or "NeGLM.scaler"
  glm_long = df |>
    dplyr::select(
      Sample, chain,
      dplyr::contains("migrationGLM.scaler"),
      dplyr::contains("NeGLM.scaler")
    ) |>
    tidyr::pivot_longer(
      cols      = -c(Sample, chain),
      names_to  = "full_name",
      values_to = "value"
    ) |>
    dplyr::mutate(
      # GLM type from column name prefix
      type = dplyr::case_when(
        stringr::str_detect(full_name, "^migrationGLM") ~ "Migration GLM",
        stringr::str_detect(full_name, "^NeGLM")        ~ "Ne GLM"
      ),
      # Covariate name: last token after the final underscore
      # e.g. "NeGLM.scaler_temperature" -> "temperature"
      cov = stringr::str_extract(full_name, "[^_]+$")
    )

  # Split into Ne and migration subsets for separate panels
  ne_long  = dplyr::filter(glm_long, type == "Ne GLM")
  mig_long = dplyr::filter(glm_long, type == "Migration GLM")

  # --- Shared violin layer helper -------------------------------------------
  # y = 0 is the null (no effect) reference line; covariates on x, chains as fill
  make_violin = function(dat, title) {
    ggplot2::ggplot(dat, ggplot2::aes(x = cov, y = value, fill = chain)) +
      ggplot2::geom_violin(alpha = 0.6, scale = "width") +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype   = "dashed",
        color      = "grey40"
      ) +
      ggplot2::scale_fill_manual(values = chain_colors, name = "Chain") +
      ggplot2::labs(
        title = title,
        x     = "Covariate",
        y     = "Scaler (effect size)"
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  p_ne  = make_violin(ne_long,  "Ne GLM Scaler Posteriors")
  p_mig = make_violin(mig_long, "Migration GLM Scaler Posteriors")

  # Stack Ne above migration; collect shared legend on the right
  patchwork::wrap_plots(p_ne, p_mig, ncol = 1) +
    patchwork::plot_layout(guides = "collect")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- P L O T  T R E E  -----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot an MCC phylogeny from a BEAST NEXUS file, colored by deme.
#
# For tips, deme is extracted from the "|"-delimited taxon label (field 3):
#   e.g. "AF202541|1999-08-28|Northeast|NY/NewYork|Unknown|unknown"
# For internal nodes, deme is the BEAST 'max' annotation (most probable
# ancestral deme from the discrete phylogeographic model).
#
# HPD intervals: 95% HPD on node heights (height_0.95_HPD) are drawn as
# horizontal range bars via ggtree::geom_range(). For tip nodes, sampling
# dates are fixed so HPD widths are negligible; uncertainty is widest at
# deep internal nodes.
#
# The x-axis is anchored to calendar years using the most recent sampling
# date (mrsd) auto-detected from taxon labels, unless overridden.
#
# tree_path : path to a BEAST MCC NEXUS tree file (.tree or .nexus)
# mrsd      : most recent sampling date as "YYYY-MM-DD"; auto-detected if NULL
plot_tree = function(beast_tree, demes = NULL, deme_colors = NULL, mrsd = NULL) {

  # --- Build deme annotation for all nodes ---------------------------------
  # Pull annotation data directly from the treedata slots to avoid
  # tidytree::as_tibble() triggering a "missing parent/label" warning when
  # it tries to reconstruct the tbl_tree edge matrix.
  #
  # @data   : data frame of BEAST annotations keyed by node number
  #           (internal nodes only; 'max' = most probable ancestral deme)
  # @phylo  : standard phylo object; tip.label holds the taxon strings
  annot_df    = beast_tree@data                # internal node annotations
  tip_labels  = beast_tree@phylo$tip.label     # taxon strings for all tips
  n_tips      = length(tip_labels)

  # Node numbers: tips are 1:n_tips, internal nodes are (n_tips+1):N
  tip_nodes   = seq_len(n_tips)

  # Auto-detect which pipe-delimited field contains deme names by finding
  # the first column whose values intersect the known demes set.
  # This replaces the hardcoded deme_loc config key.
  tip_fields  = stringr::str_split_fixed(tip_labels, "[|]", n = Inf)
  deme_col    = which(apply(tip_fields, 2, function(col) any(col %in% demes)))[1]
  if (is.na(deme_col)) {
    stop("plot_tree: no taxon label field matched any known deme. ",
         "Check deme_col_regex in config and taxon label format.")
  }
  tip_deme    = tip_fields[, deme_col]

  # Build a single data frame covering every node for %<+% attachment
  tip_meta = data.frame(
    node = tip_nodes,
    deme = tip_deme,
    stringsAsFactors = FALSE
  )

  internal_meta = data.frame(
    node = as.integer(annot_df$node),   # @data$node is character; coerce to match tip_nodes
    deme = annot_df$max,                # most probable ancestral deme (BEAST annotation)
    stringsAsFactors = FALSE
  )

  tree_meta = dplyr::bind_rows(tip_meta, internal_meta)

  # --- Auto-detect MRSD from taxon labels if not provided ------------------
  # MRSD anchors the x-axis so that tip positions equal sampling dates.
  # "YYYY-MM-DD" strings sort lexicographically so max() returns the latest.
  if (is.null(mrsd)) {
    mrsd = tip_labels |>
      purrr::map_chr(~ stringr::str_split(.x, "[|]")[[1]][[2]]) |>
      max()

    message("Auto-detected MRSD: ", mrsd)
  }

  # --- Recode compound BEAST annotations to NA ----------------------------
  # Internal nodes can carry compound states like "South+NorthCentral"
  # when BEAST is uncertain between two demes. These are not discrete
  # locations; recode them to NA so they fall through to na.value (grey).
  tree_meta = dplyr::mutate(
    tree_meta,
    deme = dplyr::if_else(grepl("+", deme, fixed = TRUE), NA_character_, deme)
  )

  # --- Palette ---------------------------------------------------------------
  # Use caller-supplied deme_colors (from the shared deme_colors target) so all
  # figures use an identical deme → color mapping. Fall back to building from
  # the demes vector if called outside the pipeline.
  if (is.null(deme_colors)) {
    if (is.null(demes)) demes = sort(unique(na.omit(tree_meta$deme)))
    deme_colors = build_deme_palette(demes)
  }

  # --- Plot -----------------------------------------------------------------
  # %<+% attaches tree_meta (keyed by 'node') to the ggtree object so that
  # the 'deme' column is available as a ggplot aesthetic across all layers.
  # Called via ggtree::`%<+%`() because the operator lives in ggtree's
  # namespace and is not available unqualified when ggtree is not library()-loaded.
  # mrsd converts the default root-distance x-axis to calendar years.
  # Suppress ggtree error on negative branch lengths — a known BEAST2 MCC
  # tree artifact where independently computed median node heights can place
  # a parent younger than its child.
  options(ignore.negative.edge = TRUE)
  p_base = ggtree::ggtree(beast_tree, mrsd = mrsd)
  ggtree::`%<+%`(p_base, tree_meta) +
    ggplot2::aes(color = deme) +

    # 95% HPD on node heights as horizontal range bars.
    # Bars are narrow for tips (fixed dates) and widest for deep internal nodes.
    ggtree::geom_range(
      range = "height_0.95_HPD",
      color = "grey50",
      alpha = 0.4,
      size  = 1.5
    ) +

    # Tip circles colored by deme
    ggtree::geom_tippoint(size = 1.5) +

    ggplot2::scale_color_manual(
      values   = deme_colors,
      na.value = "grey70",
      name     = "Deme"
    ) +

    # theme_tree2() adds a clean x-axis (year scale) without grid clutter
    ggtree::theme_tree2() +
    ggplot2::labs(x = "Year")
}
