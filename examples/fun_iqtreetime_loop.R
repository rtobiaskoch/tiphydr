iqtreetime_loop <- function(lineages, base_dir = "3_output/lineages") {
  # Ensure treetime output directories exist
  dir.create(file.path(base_dir, "treetime"), showWarnings = FALSE)
  
  # Loop through each lineage
  purrr::walk(lineages, function(lin) {
    fasta_file <- file.path(base_dir, paste0(lin, ".fasta"))
    treefile   <- paste0(fasta_file, ".treefile") # IQ-TREE output
    
    message("Processing lineage: ", lin)
    
    # ---- Run IQ-TREE ----
    cmd_iqtree <- sprintf(
      "bash -lc 'iqtree -s %s -m GTR+G -bb 1000 -alrt 1000 -bnni -nt AUTO'",
      fasta_file
    )
    system(cmd_iqtree)
    
    # ---- Run TreeTime ----
    # Modify metadata path or skyline params as needed
    cmd_treetime <- sprintf(
      "bash -lc 'treetime --aln %s --tree %s --dates 3_output/%s_metadata_treetime.tsv --coalescent skyline --stochastic-resolve --n-skyline 96 --outdir %s/treetime/%s'",
      fasta_file,
      treefile,
      lin,
      base_dir,
      lin
    )
    
    dir.create(file.path(base_dir, "treetime", lin), showWarnings = FALSE)
    system(cmd_treetime)
  })
}
