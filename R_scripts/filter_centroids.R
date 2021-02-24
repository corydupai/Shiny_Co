# Filter centroid.fasta output from panoct to only include orthogroups that contain 10% of genomes.
# Save the centroids for these groups to a .fasta file that Kallisto.sh uses to create a species 
# transcriptome index.

library(data.table)
library(tidyverse)
library(seqinr)

filter_centroid_fasta <- function(centroids, 
                            good_ids){
  if(getName(centroids) %in% good_ids){
    return(centroids)
  } else {
    NA
  }
}

filter_centroids <- 
  function(basedir, perc = 0.1){
    panoct_dir <- paste0(basedir,"panoct_output/pangenome/results/")
    centroids_in <- paste0(panoct_dir, "centroids.fasta")
    matchtable_in <- paste0(panoct_dir, "matchtable_0_1.txt")
    matchtable_full <- paste0(panoct_dir, "matchtable.txt")
    gl_in <- paste0(panoct_dir, "genomes.list")
    data4app <- paste0(basedir,"data4app/")
    dir.create(data4app)
    
    # Read in centroids file and matchtable. Count orthologs for each centroid
    # in the matchtable and then exclude any centroids with orthologs in 
    # < 10% of genomes.
    centroids.fa <- read.fasta(centroids_in)
    
    matchtable <- fread(matchtable_in) %>%
      column_to_rownames(var = "V1") %>%
      rowSums()
    max_orthos <- max(matchtable)
    matchtable <- matchtable[matchtable>=(max_orthos * perc)]
    
    matchtable_dt <- data.table(
      centroid = paste0("centroid_", names(matchtable)),
      counts = matchtable
    )
    
    
    c_out <- lapply(centroids.fa,
                    filter_centroid_fasta,
                    good_ids = matchtable_dt$centroid)
    c_out <- c_out[!is.na(c_out)]
    
    # This section filters the locus table and outputs it
    # to the appropriate dir for shinyapp files
    locus_table <- 
      fread(matchtable_full,
            na.strings = "----------"
      )[as.integer(names(matchtable))
        ][, Centroid := paste0("centroid_",V1)
          ][, !c("V1")]
    gl <- fread(gl_in, header = FALSE) %>%
      unlist() %>%
      unname()
    
    colnames(locus_table) <- 
      c(gl,"centroid")
    
    # Write files
    fwrite(locus_table,
           paste0(data4app,"locus_table.csv"))
    
    file.copy(gl_in, paste0(data4app,"genome_list.csv")) # Copy genome_list for us in shinyapp
    
    write.fasta(sequences = getSequence(c_out),
                names = getName(c_out),
                file.out = paste0(basedir,
                                  "/filtered_transcripts.fasta"))
    return("Done")
}

basedir <- "/stor/work/Wilke/tingle/WGCNA_expanded/bacteria/E.coli/"
# basedir <- "/stor/work/Wilke/tingle/WGCNA_expanded/bacteria/K.pneumoniae/"

filter_centroids(basedir, perc = 0.1)
