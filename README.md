# tiphydr — Tidy Phylodynamics in R

An R package for fast, tidy manipulation of sequence data and associated metadata. Functions are designed to compose left-to-right in a pipeline: read → combine → nest with metadata → align and trim → assign lineages.

---

## Installation

```r
# Install Bioconductor dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

# Install tiphydr from GitHub
devtools::install_github("rtobiaskoch/tiphydr")
```

MAFFT must be installed on your system for `fasta_trim_ref()`:
- macOS: `brew install mafft`
- Linux: `sudo apt install mafft`
- Or download from <https://mafft.cbrc.jp/alignment/software/>

---

## Functions

### `fasta_read(paths, verbose = TRUE)`

Read one or more FASTA files. Returns a `DNAStringSet` for a single path, or a named list of `DNAStringSet` objects for multiple paths.

```r
library(tiphydr)

# Single file → DNAStringSet
seqs <- fasta_read("data/wnv_sequences.fasta")

# Multiple files → named list
files <- fasta_read(c("data/batch1.fasta", "data/batch2.fasta"))
```

---

### `fasta_combine(fasta_list, verbose = TRUE)`

Concatenate a list of `DNAStringSet` objects into one. Warns on duplicate sequence names.

```r
# Combine the two batches from above
all_seqs <- fasta_combine(files)

# Or combine ref + sequences
combined <- fasta_combine(list(ref, seqs))
```

---

### `fasta_nest(df, biostring, id_col, delim, verbose = TRUE)`

Join a `DNAStringSet` to a metadata dataframe as a `sequence` list-column — one sequence per row, analogous to an `sf` geometry column. Matching splits sequence names on `delim` and looks for an exact match to the `id_col` value in any resulting part.

```r
metadata <- read.csv("data/metadata.csv")
seqs     <- fasta_read("data/sequences.fasta")

# FASTA names like "WNV|2021|CO|CSU123"
# metadata has column strain_id = "CSU123"
nested <- fasta_nest(metadata, seqs, id_col = "strain_id", delim = "|")

# sequence column is now a DNAStringSet list-column
nested$sequence[[1]]

# Works with dplyr verbs
library(dplyr)
nested |>
  filter(location == "Colorado") |>
  mutate(year = lubridate::year(date))
```

---

### `fasta_unnest(nested_df, verbose = TRUE)`

Inverse of `fasta_nest()`. Extracts the `sequence` list-column back into a standalone `DNAStringSet`, dropping rows with no matched sequence.

```r
result <- fasta_unnest(nested)
result$metadata   # tibble without sequence column
result$biostring  # DNAStringSet
```

---

### `fasta_trim_ref(fasta, ref, alignment = NULL, drop = FALSE, verbose = TRUE)`

Align sequences to a reference using MAFFT and trim to reference coordinates. After trimming, alignment position N equals ungapped reference position N — the coordinate invariant required by `define_lineage()`.

```r
ref  <- fasta_read("data/wnv_ref.fasta")
seqs <- fasta_read("data/wnv_sequences.fasta")

# Full alignment from scratch
trimmed <- fasta_trim_ref(seqs, ref)

# Add new sequences to an existing alignment (MAFFT --add profile mode)
# existing columns are preserved exactly
updated <- fasta_trim_ref(new_seqs, ref, alignment = trimmed)

# Remove sequences not in new_seqs from output
updated_clean <- fasta_trim_ref(new_seqs, ref, alignment = trimmed, drop = TRUE)
```

---

### `define_lineage(alignment, ref, muts, type = "nuc", verbose = TRUE)`

Assign lineage labels to sequences based on diagnostic mutations. Matching is strict — all mutations for a lineage must be present. When a sequence matches multiple lineages, the most specific (most mutations required) is assigned.

```r
# Mutation table: long format with lineage, pos (reference coordinates), residue
muts <- data.frame(
  lineage = c("NY10",  "NY10",  "NS2A_only"),
  pos     = c(1331L,   2513L,   1331L),
  residue = c("K",     "M",     "K")
)

# type = "aa": translate from first ATG in ref, check amino acid positions
result <- define_lineage(trimmed, ref, muts, type = "aa")
head(result)
#>                  strain lineage
#> 1  WNV|2021|CO|CSU001    NY10
#> 2  WNV|2020|CO|CSU002    NS2A_only
#> 3  WNV|2019|CO|CSU003    unknown

# type = "nuc" (default): check nucleotide positions directly
nuc_muts <- data.frame(
  lineage = c("LineageA", "LineageA"),
  pos     = c(10L, 25L),
  residue = c("A", "C")
)
result_nuc <- define_lineage(trimmed, ref, nuc_muts, type = "nuc")
```

---

## Full Pipeline Example

```r
library(tiphydr)
library(dplyr)

# 1. Read data
ref      <- fasta_read("data/wnv_ref.fasta")
seqs     <- fasta_read("data/wnv_sequences.fasta")
metadata <- read.csv("data/metadata.csv")
muts     <- read.csv("data/lineage_mutations.csv")

# 2. Align and trim to reference coordinates
trimmed <- fasta_trim_ref(seqs, ref)

# 3. Assign lineages
lineages <- define_lineage(trimmed, ref, muts, type = "aa")

# 4. Join lineages to metadata and nest sequences
metadata_with_lineage <- metadata |>
  left_join(lineages, by = c("strain_id" = "strain")) |>
  fasta_nest(trimmed, id_col = "strain_id", delim = "|")

# 5. Filter to a lineage and extract sequences for downstream analysis
ny10_data   <- filter(metadata_with_lineage, lineage == "NY10")
ny10_unnest <- fasta_unnest(ny10_data)

writeXStringSet(ny10_unnest$biostring, "output/ny10_sequences.fasta")
write.csv(ny10_unnest$metadata, "output/ny10_metadata.csv", row.names = FALSE)
```

---

## Test Data

The package ships with a small synthetic WNV dataset in `inst/extdata/` for testing and examples:

```r
ref_path  <- system.file("extdata", "wnv_ref.fasta",  package = "tiphydr")
seq_path  <- system.file("extdata", "wnv_seqs.fasta", package = "tiphydr")
meta_path <- system.file("extdata", "wnv_metadata.csv", package = "tiphydr")
muts_path <- system.file("extdata", "wnv_muts.csv",   package = "tiphydr")

ref      <- fasta_read(ref_path)
seqs     <- fasta_read(seq_path)
metadata <- read.csv(meta_path)
muts     <- read.csv(muts_path)
```

Sequences are 84 nt with two mutation sites (positions 10 and 25) encoding NY10 and PARTIAL_A lineages.

---

## Package Structure

```
R/
  fasta_read.R       — FASTA file I/O
  fasta_combine.R    — concatenate DNAStringSets
  fasta_nest.R       — join sequences to metadata dataframe
  fasta_unnest.R     — extract sequences from nested dataframe
  fasta_trim_ref.R   — align and trim to reference (MAFFT)
  define_lineage.R   — assign lineage labels by diagnostic mutations
```

---

## Dependencies

| Package | Role |
|---------|------|
| `Biostrings` (Bioconductor) | Sequence representation and operations |
| `ape` | DNAbin conversion for MAFFT input |
| `ips` | MAFFT wrapper |
| `purrr` | Functional iteration |
| `stringr` | String splitting for sequence name matching |
| `dplyr` | Dataframe operations |
| `tibble` | Tibble construction |
