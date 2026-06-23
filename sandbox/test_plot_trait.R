library(dplyr)
devtools::load_all(".")

tree <- ape::read.tree("temp/test25_dated_dayadj.timetree.nwk")
node_probs <- read.delim("temp/test25_dta_node_probs.tsv")
f <- Biostrings::readDNAStringSet("temp/test25.fasta")

tip_dates <- extract_metadata(f) %>%
  select(taxa, V2) %>%
  rename(
    tip_label = taxa,
    date = V2
  ) %>%
  mutate(date = lubridate::ymd(date))

tree_df <- build_tree_df(tree, node_probs, tip_dates)
plot_tree_trait(tree_df, "inferred_state")

e <- explode_tree(tree, node_probs)
intro <- e$introductions %>%
  group_by(intro_clade_id) %>%
  mutate(
    persistence = if_else(
      last_sample_date - inferred_intro_date > 1,
      "persistent",
      "transient"
    )
  )
plot_explode_tree(e$introductions)
