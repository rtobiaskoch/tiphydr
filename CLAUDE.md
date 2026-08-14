Repo: TIdy PHYlodynamics in R TiPhydR ("Ti-fighter")

PURPOSE:
develop suite of R functions for:
1. Faster manipulation:
  a. more efficient cleaning of sequence data that often contains associated metadata
      that needs to be filtered, replaced etc based off of the sequence file name.

2. Working with Mutations
  a. get positions of particular amino acid sequences. useful for start and stop codons
  a. use mutation for lineage calling

3. Working with Trees:
  a. manipulate, split and modify trees based off ancestral node and tip info.

OUT-OF-SCOPE
modeling, imputing and altering raw data 

To prevent reinventing the wheel ensure no current functions exist that accomplish the same thing I am proposing
Use Tidyverse styling

class preference
convert data to dataframes when possible 

Packages:
manipulation: tidyverse
sequences: Biostrings
plotting: ggplot
phylo: treeio, tidytree


