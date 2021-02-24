# Need to add code to download Bioconductor packages if necessary or specify that somehow
library(seqinr)
library(data.table)
library(tidyverse)
library(DOSE)
library(parallel)
library(clusterProfiler)
library(KEGGREST)
library(qdapTools)
library(OpenImageR)

#' Calculates enrichment for specific term.
#' 
#' @param df Dataframe containing terms of interest
#' @param column_name Column containing GO or KEGG terms
#' @return Dataframe with enrichment info
#' 
get_genelist <- function(df, column_name){
  
  g_sym = ";"
  
  column_name <- sym(column_name)
  
  GL <- df %>%
    select(centroid,
           dynamicColors,
           !!column_name) %>%
    mutate(term_col = as.list(str_split(!!column_name, g_sym))) %>%
    unnest(cols = c(term_col)) %>%
    filter(!is.na(term_col) &
             term_col != "") %>%
    select(centroid,
           dynamicColors,
           term_col)
}


#' Calculates enrichment for specific term.
#' 
#' @param dc Dynamic color, AKA subnetwork name
#' @param df Dataframe containing terms of interest
#' @param Domain What domain are we looking at, e.g. KEGG, Biological function, Cellular Localization
#' @return Dataframe with enrichment info
#' 
get_enriched <-
  function(dc, df, Domain) {
    en_res <-
      enricher(gene = df %>%
                 filter(dynamicColors == dc) %>%
                 select(gene = centroid) %>%
                 unlist(),
               TERM2GENE = df %>% select(
                 term = term_col,
                 gene = centroid)
      ) %>%
      as_tibble() %>%
      mutate(dynamicColors = dc,
             Domain = Domain)
  }




#' Calculates gene enrichment scores for KEGG and GO terms and writes output.
#' 
#' @param basedir The base directory.
#' 
get_enriched_terms <- function(basedir){

kegg_annotations <- paste0(basedir, "filtered_annotations.csv")
graph_rds <- paste0(basedir,"data4app/filtered_graph_0.10.RDS")
go_annotations <- paste0(basedir, "GOterms.txt")
match_table <- paste0(basedir, "panoct_output/pangenome/results/matchtable.txt")
centroid_annotations <- paste0(basedir, "cutoff_9.fasta")

# Format KEGG terms
## Merge centroid|dynamicColors with centroid|annotations|protein
# Get Graph Information & turn into dataframe 
dynamicColors_df <- readRDS(graph_rds) %>% 
  activate(nodes) %>% 
  as.data.frame() %>% 
  rename(centroid = name)

# Format kegg_annotations, merge with dynamicColors,  seperate by pathway and 
kegg_filtered <- read.csv(kegg_annotations) %>% 
  select(-X) %>% 
  group_by(centroid) %>% 
  mutate(annotations = unique(annotations)) %>% 
  # summarise(annotations = unique(annotations)) %>% 
  merge(dynamicColors_df, by = "centroid")

# Split into two dataframes - one for pathway codes and the other for enzyme codes
kegg_enzymes <- kegg_filtered %>% 
  filter(stringr::str_detect(annotations, stringr::regex("[:blank:]\\d\\."))) %>% 
  rename(enzyme = annotations)

kegg_pathway <- kegg_filtered %>% 
  filter(!stringr::str_detect(annotations, stringr::regex("[:blank:]\\d\\."))) %>% 
  rename(pathway = annotations)

# Merge centroid|dynamicColors with loci|Go Terms
# Get centroids and annotations datafarme
centroids_to_annotation <- list2df(read.fasta(centroid_annotations)) %>% 
  select(-X1) %>% 
  mutate(centroid = X2) %>% 
  separate(`X2`, c("centroid_text", "id"), "_" ) %>% 
  select(-centroid_text) %>% 
  group_by(id) %>% 
  summarise(id_number = unique(id)) %>% 
  select(id_number)


# Merge centroids & match table to filter; make longer
filter_match_table <- read.table(match_table, header = FALSE) %>% 
  rename(id_number = V1) %>% 
  merge(centroids_to_annotation, by="id_number") %>% 
  pivot_longer(-id_number) %>% 
  select(-name) %>% 
  rename(loci = value)

# This has locus tag | Go Terms
# Merge kp_GOterms with filter_match_table which has id # | locus_tag
# Then, Merge with dy_col_df which has centroid_# | dynamicColor 
GO_terms <- read.delim(go_annotations,
                       sep="\t", 
                       header=TRUE, 
                       comment.char="#", na.strings=".", 
                       stringsAsFactors=FALSE,
                       quote="", 
                       fill=TRUE) %>% 
  mutate(loci = ID) %>% 
  merge(filter_match_table, by = "loci") %>% 
  mutate(centroid = paste0("centroid_", id_number)) %>% 
  select(-id_number) %>% 
  merge(dynamicColors_df, by = "centroid") %>% 
  separate_rows(GOTERM_BP_DIRECT, sep=",") %>% 
  separate_rows(GOTERM_CC_DIRECT, sep=",") %>% 
  separate_rows(GOTERM_MF_DIRECT, sep=",")


## Get KEGG Enrichments
kegg_pathway_list <- get_genelist(kegg_pathway, "pathway")
kegg_enzymes_list <- get_genelist(kegg_enzymes, "enzyme")


net_colors <- kegg_filtered$dynamicColors %>%
  unique() %>%
  as.character() %>%
  unlist()

KEGG_terms <-
  bind_rows(
    lapply(net_colors, 
           FUN = get_enriched,
           df = kegg_pathway_list,
           Domain = "KEGG Pathways"),
    lapply(net_colors,
           FUN = get_enriched,
           df = kegg_enzymes_list,
           Domain = "KEGG Enzymes")
  ) %>%
  # This filters and formats the output
  mutate(full_counts =
           gsub(".*\\/", "", GeneRatio) %>%
           as.numeric(),
         ratio = Count / full_counts,
         Term = ID,
         Bg_Counts = gsub("\\/.*", "", BgRatio) %>%
           as.numeric(),
         Bg_Totals = gsub(".*\\/", "", BgRatio) %>%
           as.numeric(),
         BR = Bg_Counts/Bg_Totals,
         qvalue = if_else(is.na(qvalue),0,
                          qvalue)) %>%
  filter(qvalue <= 0.05)

## Get GO Enrichments
go_bp <- get_genelist(GO_terms, "GOTERM_BP_DIRECT")
go_cc <- get_genelist(GO_terms, "GOTERM_CC_DIRECT")
go_mf <- get_genelist(GO_terms, "GOTERM_MF_DIRECT")

GO_terms <-bind_rows(
  mclapply(net_colors, 
           FUN = get_enriched,
           df = go_bp,
           Domain = "Biological Process"),
  mclapply(net_colors,
           FUN = get_enriched,
           df = go_cc,
           Domain = "Cellular Component"),
  mclapply(net_colors,
           FUN = get_enriched,
           df = go_mf,
           Domain = "Molecular Function")
) %>%
  # This filters and formats the output
  mutate(full_counts =
           gsub(".*\\/", "", GeneRatio) %>%
           as.numeric(),
         ratio = Count / full_counts,
         Term = ID,
         Bg_Counts = gsub("\\/.*", "", BgRatio) %>%
           as.numeric(),
         Bg_Totals = gsub(".*\\/", "", BgRatio) %>%
           as.numeric(),
         BR = Bg_Counts/Bg_Totals,
         qvalue = if_else(is.na(qvalue),0,
                          qvalue)) %>%
  filter(qvalue <= 0.05) 

subnetwork_out<- KEGG_terms %>%
  rbind(GO_terms) %>%
  mutate(`-log(qvalue)`= -log(qvalue))

write.csv(subnetwork_out, paste0(basedir,"data4app/enriched_terms.csv"))
return("done")
}


#' Gets KEGG terms for all genes. Also outputs loci for analysis with DAVID to get GO terms
#' 
#' @param loci Locus of interest
#' @return KEGG terms associated with locus of interest
#' 
get_annotations.func <- function(loci)
{
  print(loci)
  gene.df <- list2df(keggGet(loci))[1,]
  return(gene.df)
}

#' Gets KEGG terms for all genes. Also outputs loci for analysis with DAVID to get GO terms
#' 
#' @param basedir The base directory.
#' 
GO_and_KEGG <- function(basedir){
  # Get Bacteria Directory 
  ```{r}

  bacdir <- "K.pneumoniae/"
  bac_initials <- "kp"
  bac_name <- "Klebsiella pneumoniae"
  ```
  
  # Get Bacterial Files
  cutoff_9 <- list2df(read.fasta(paste0(basedir, "cutoff_9.fasta")))
  match.table <- read.table(paste0(basedir, "panoct_output/pangenome/results/matchtable.txt"), header = FALSE) 
  
  # Step 1:
  ## Filter `match.table` based on Filtered centroids.fasta file (cutoff_9)
  # `match.table` has a list of genes associated with each ortholog used to construct the pangenome
  # `cutoff_9` filtered out orthologs from `match.table` with less than 9 genes 
  
  # Fornat centroids file 
  centroids_to_annotation <- cutoff_9 %>% 
    select(-X1) %>% 
    mutate(centroid = X2) %>% 
    separate(`X2`, c("centroid_text", "id_number"), "_" ) %>% 
    select(-centroid_text) %>% 
    summarise(id_number = unique(id_number), centroid = unique(centroid)) 
  
  # Merge centroids & match table to filter; make longer
  filter_match_table <- match.table %>% 
    rename(id_number = V1) %>% 
    merge(centroids_to_annotation, by="id_number") %>% 
    pivot_longer(-id_number) %>% 
    select(-name) %>% 
    rename(loci = value)
  
  # Clean up environment
  remove(cutoff_9, centroids_to_annotation)
  
  # Step 2: Find genes & corresponding function
  ## Get KEGG Identifiers -- (https://www.genome.jp/kegg/catalog/org_list.html)
  # In the KEGG database, each reference geneome for each organism has an associated `T()` code. 
  # The `keggList()` function pulls gene loci identifiers affiliated with the `T()`, 
  # formatted as `organism_code:gene_loci` along with the protein annotation. 
  # This dataframe of gene loci is then filtered based off of loci in the `filter_match_table`, 
  # and stored as `bac_entries`. 
  
  #"I've written the `get_annotations.func()` to call `keggGet()` on each loci, which retrieves all associated database entries. **Unfortunately, I can't parallelize this function call--I think there's some issue with how many times I can ping the database. So this takes a veeerrry long time. Ideas** The annotations are saved as a .csv file for future reference. Then we filter these entries to only include Kegg Codes and save this as a `~filtered_annotations.csv`. "

  # Get T() Number for all reference genomes & call keggList() for organism name:gene loci keys. 
  t_codes <- as.vector(read.csv(paste0(basedir, "tcodes.csv"))$x) 
  
  # Get all org_id:loci associated with each t_code 
  reference_t <- data.frame()
  for(t in t_codes)
  {
  kegg_t <- list2df(keggList(t)) %>% 
  transmute(`protein_annotation` = X1, `org_loci` = X2)
  
  reference_t <- rbind(reference_t, kegg_t)
  }
  
  # Filter entries in reference_t by gene loci in filter_match_table
  bac_entries <- reference_t %>% 
  mutate(temp = org_loci) %>% 
  separate(temp, c("org", "loci"), sep = ":") %>% 
  merge(filter_match_table, by="loci") 
  
  # Function to get annotations and format into dataframe
  # Get annotations for each gene:loci combination
  
  ############# WRITE ####################
  bac_annotations <- lapply(bac_entries$org_loci, get_annotations.func)
  # 
  bac_annotations.df <- list2df(bac_annotations)
  # 
  write.csv(bac_annotations.df,paste0(basedir, "annotations.csv"))
  # 
  # bac_annotations.df <- read.csv(paste0(basedir, bacdir, "annotations.csv"))
  # 
  bac_filtered <- bac_annotations.df %>%
     group_by(X2) %>%
     mutate(loci = X1[1], protein = X1[2], centroid = paste0("centroid_", X2)) %>%
     filter(stringr::str_detect(X1, stringr::regex("[:blank:]\\d\\d")) |
            stringr::str_detect(X1, stringr::regex("[:blank:]\\d\\."))) %>%
     filter(!stringr::str_detect(X1, stringr::regex(bac_name))) %>%
     rename(annotations = X1, id_number = X2) %>%
     merge(bac_entries, by = c("id_number")) %>%
     select(-id_number, -org_loci, -org, -loci.x) %>%
     select(centroid, loci.y, protein, protein_annotation, annotations) %>%
     arrange(centroid) %>%
    group_by(centroid, annotations) %>%
    summarise(protein = unique(protein))
  # 
  write.csv(bac_filtered, paste0(basedir, "filtered_annotations.csv"))
  

  
  
  ## Get GO Annotations
  # To access the GO annotations for each loci, I downloaded the list of loci 
  # locally and uploaded it to the DAVID website (https://david.ncifcrf.gov/). 
  # DAVID had GO terms for 5,227 of the 397,052 loci (which I suppose is reasonable since klebsiella may 
  # just not be that extensively annotated). *Does anyone know how to do this straight from R? 
  # It doesn't take long...just annoying*. 

# Clean up list of loci from filter_match_table
bac_loci <- filter_match_table %>% 
  filter(loci != "----------") %>% 
  filter(!str_detect(loci, stringr::regex("centroid_\\d")))
  
## Write list of loci to a txt file & Download locally
write(bac_loci$loci, file = paste0(basedir, "loci.txt"))

# Upload local list of loci to DAVID website & run

## Notes - 391,825 unknown, 5227 DAVID Ids (line 1 is header) 
## Upload output txt file with GO terms
# kp_GOterms <- read.delim("/stor/work/Wilke/tingle/WGCNA_expanded/bacteria/K.pneumoniae/kp_GOterms.txt",
#                          sep="\t", 
#                          header=TRUE, 
#                          comment.char="#", na.strings=".", 
#                          stringsAsFactors=FALSE,
#                          quote="", 
#                          fill=TRUE)
return("Done")

}

