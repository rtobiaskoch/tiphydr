library(purrr)
library(dplyr)
library(Biostrings)

split_mdata_fasta <- function(df, fasta, id_col, out_dir = "3_output/lineage_outputs") {
  # Ensure output directory exists
  dir.create(out_dir, showWarnings = FALSE)
  
  # Split the dataframe by lineage
  lineage_list <- split(df, df$lineage)
  
  # Use map to iterate over each lineage subset
  map(names(lineage_list), function(lin) {
    df_sub <- lineage_list[[lin]]
    
    
    #format for treetime
    df_date <- df_sub %>% 
      transmute(strain = csu_id,
                date = trap_date)
    
    # Get IDs to match fasta names
    ids <- df_sub[[id_col]]
    
    # Match and extract sequences
    seqs_sub <- fasta[names(fasta) %in% ids]
    
    # Define output paths
    df_path <- file.path(out_dir, paste0(lin, ".csv"))
    df_date_path <- file.path(out_dir, paste0(lin, "_dates.tsv"))
    fasta_path <- file.path(out_dir, paste0(lin, ".fasta"))
    
    # Save outputs
    write.csv(df_sub, df_path, row.names = FALSE)
    write.table(df_date, df_date_path, quote=FALSE, sep='\t', row.names = F)
    writeXStringSet(seqs_sub, fasta_path)
    
    # Return summary info
    tibble(lineage = lin,
           n_rows = nrow(df_sub),
           n_seqs = length(seqs_sub),
           csv = df_path,
           fasta = fasta_path)
  }) %>%
    bind_rows()
}
