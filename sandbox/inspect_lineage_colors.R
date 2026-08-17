# Visual inspection of assign_lineage_colors() output against the real
# wnv_mut_orf table — one swatch per lineage, in the order they appear in
# wnv_mut_orf (family -> partials -> cross-family hybrids), so the
# darkest-to-lightest shading within a family and the hybrid blends are
# checkable by eye.

devtools::load_all(".")
library(ggplot2)

palette <- assign_lineage_colors(wnv_mut_orf)

# Lock in wnv_mut_orf's own row order (family, then its partials, then the
# cross-family hybrids involving it) rather than ggplot's default alphabetical
# ordering, so related shades sit next to each other.
palette$lineage <- factor(palette$lineage, levels = rev(palette$lineage))

swatch_plot <- ggplot(palette, aes(x = 1, y = lineage, fill = color)) +
  geom_tile(width = 1, height = 0.9, color = "white", linewidth = 1) +
  geom_text(
    aes(label = paste(lineage, color)),
    hjust = 0,
    x = 1.55,
    size = 3.5
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.4, 4), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10)) +
  labs(title = "assign_lineage_colors(wnv_mut_orf)")

write.csv(palette, "sandbox/lineage_colors.csv", row.names = FALSE)
ggsave(
  "sandbox/lineage_colors_preview.png",
  swatch_plot,
  width = 6,
  height = 8,
  dpi = 150
)
