---
editor_options: 
  markdown: 
    wrap: 80
---

# Functions to Build
--------------------------------------------------------------------------------
## split_gene.R
take gff and generate separate fasta files for each gene for input sequences

## fasta_name_switch.R
### input:
- DNAstringset 
- match pattern
- replacement name

## explode_tree.R  ✅ DONE
Implemented in R/explode_tree.R with reusable helpers:
most_likely_deme.R, node_dates_from_timetree.R, build_tree_df.R,
detect_introductions.R. Detects introductions into all demes using the
continued-transmission definition (reintroductions counted separately) and
returns a tidy introductions table plus split subtrees per multi-tip clade.


## build_tree_df.

## plot_explode_tree.R 
Take output of explode_tree.R x$introductions and make plot using wnv-foco_phylo_intro_plots.Rmd as an example.
Final output should be the p_intro2 plot. 
group_by intro_clade_id and deme
facet by deme
add color of inferred_intro_date point to be the source inferred_intro_source


## DEFINE_LINEAGES.R
- [x] create a simpler function just with the input being the absolute amino acid position and residue as required input. User can use the resolve positions function
- [x] save existing function as define_lineages_dev.R in sandbox



