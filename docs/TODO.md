---
editor_options: 
  markdown: 
    wrap: 80
---

# Functions to Build

## fasta_read.R 
### desc: 
-reads in list of sequences

### input: 
- path to sequences

### output: 
- list of biostring objects

---

## fasta_combine.R
### desc:
-   combines list of biostring objects into a single biostring

### input:
-   list of biostring objects \### output:
-   combined biostring object


## extract_metadata.R
### desc:
-   takes metadata in fasta header (eg. ID|date|col) and creates a dataframe

### input:
-   BIostrings
-   list : list of column names for metadata
--------------------------------------------------------------------------------

## fasta_nest.R
### desc:
-   takes dataframe with metadata and biostrings data and nests biostring data
    as a column like sf geo object does in the metadata so you can easily manipulate
    using tidyverse

### input
-   df: dataframe with metdata
-   biostring: biostring data
-   id : column name that has the unique iq pattern match for the sequences name

--------------------------------------------------------------------------------

## fasta_unnest.R
### desc:
-   splits biostrings data metadata

--------------------------------------------------------------------------------

## fasta_trim_ref.R:

### desc:
-   aligns then trims to ref ensure that all sequences are the same length
-   align using mafft

### input:
-   fasta (required) : sequences to be aligned
-   ref(required) : reference sequences to trim to
-   alignment (optional) : existing alignment to use to computational efficiency
-   drop (optional) : logical for sequences in existing alignment that are not
    in fasta to be dropped. T = drop

### messages:
-   number of sequences dropped or added if alignment provided
- stats of original avg length then new trimmed length

--------------------------------------------------------------------------------

## define_lineage.R (use NY10.R as a guide)

### DESC

-   provide list of mutations with associated lineage to return dataframe with the fasta name strain and the lineage

### input: 
- alignment (required) : using an aligned file with all sequences that
are the same length 
- ref (required) : reference fasta with to confirm 
- muts (required) : dataframe with lineage then associated mutations long format
example: dataframe(lineage = c("NY10", "NY10""), pos = c(1331, 2513), aa = c("K", "M"))
- type (optional): defaulta "nuc". options c("aa",
"nuc")

### output:
- dataframe: with strain and lineage

### messages:
- if no matches return "unknown" and print list of number of unknown lineages
-print table with lineage and number ID'd


--------------------------------------------------------------------------------

#dont do yet 
## explode_tree.R
### input: 
- tree 

#dont do yet
#plot_explode_tree.R
