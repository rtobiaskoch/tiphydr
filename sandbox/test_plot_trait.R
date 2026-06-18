tree = ape::read.tree("temp/test25_dated_dayadj.timetree.nwk")
node_probs = read.delim("temp/test25_dta_node_probs.tsv")
f = Biostrings::readDNAStringSet("temp/test25.fasta")

tip_dates = extract_metadata(f) %>%
  select(taxa, V2) %>%
  rename(
    tip_label = taxa,
    date = V2
  ) %>%
  mutate(date = lubridate::ymd(date))

tree_df <- build_tree_df(tree, node_probs, tip_dates)
plot_tree_trait(tree_df, "inferred_state")
