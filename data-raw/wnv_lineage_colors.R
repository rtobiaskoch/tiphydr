# Generates data/wnv_lineage_colors.rda: default base-family colors for
# assign_lineage_colors(). Derived from lineage_colors.csv (a Dark2-style
# palette), dropping the standalone NS2A/NS4B rows from that source file —
# under assign_lineage_colors()'s shading scheme, NY10_NS2A/NY10_NS4B are
# lighter tints of NY10's own green rather than independent hues, so no
# separate base color is needed for them. "unknown"/"ambiguous" are handled
# by assign_lineage_colors()'s own dedicated parameters, not this table.

wnv_lineage_colors <- tibble::tribble(
  ~lineage, ~color,
  "NY10",   "#1B9E77",
  "WN02",   "#D95F02",
  "SW03",   "#7570B3",
  "NY99",   "#E7298A"
)

usethis::use_data(wnv_lineage_colors, overwrite = TRUE)
