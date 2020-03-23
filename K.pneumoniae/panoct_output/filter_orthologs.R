library(dplyr)
library(tidyverse)
library(ggplot2)
library(seqinr)
library(largeList)
library(qdapTools)
library(stringr)
# library(bedr)
# library(devtools)
# library(Smisc)

# Read in match table
match.table <- read.table("/Users/t/Documents/GitHub/WGCNA_expanded/K.pneumoniae/panoct_output/pangenome/results/matchtable.txt", header=FALSE)

# Plot total genes/ortholog and total genes in set of orthologs to find cutoffs
match.table.long <- match.table %>% 
  mutate(ortholog = V1) %>% 
  select(-V1) %>% 
  pivot_longer(-ortholog, names_to="genome",
               values_to ="gene") %>% 
  mutate(gene.exists = if_else(grepl("----------",gene),0,1)) %>% 
  group_by(ortholog) %>% 
  summarize(total.genes=sum(gene.exists)) %>% 
  group_by(total.genes) %>% 
  summarize(total.orthologs=n()) %>% 
  mutate(summed = cumsum(total.orthologs)) %>% 
  # ggplot(aes(x=total.genes, y=total.orthologs)) + 
  #   geom_point() + scale_y_log10() 
  ggplot(aes(x=total.genes)) +
    geom_point(aes(y=total.orthologs, color="total orthologs")) +
    geom_point(aes(y=summed, color="total genes")) + scale_y_log10(name="log10 counts")

# Get list of orthologs based on cutoffs
cutoff = 9

ortholog.list <- match.table %>% 
  mutate(ortholog = V1) %>% 
  select(-V1) %>% 
  pivot_longer(-ortholog, names_to="genome",
               values_to ="gene") %>% 
  mutate(gene.exists = if_else(grepl("----------",gene),0,1)) %>% 
  group_by(ortholog) %>% 
  summarize(total.genes = sum(gene.exists)) %>% 
  filter(total.genes>cutoff) %>% 
  select(-total.genes) %>% 
  mutate(ortholog = strtoi(ortholog, 10L))

  ortholog.list <- unlist(ortholog.list[,1])

# Read centroids.fasta file
centroids <- read.fasta(file="/Users/t/Documents/GitHub/WGCNA_expanded/K.pneumoniae/panoct_output/pangenome/results/centroids.fasta",
                        seqtype = "DNA", as.string = TRUE)

  #typeof(centroids)
  #df <- data.frame(t(sapply(centroids,c)))
centroid.df <- list2df(centroids)

filter.centroid <- centroid.df %>% 
  mutate(chargaff = X1, centroid = strtoi(str_sub(X2, 10), 10L)) %>% 
  select(-X1, -X2) %>% 
  group_by(centroid) %>% 
  filter(centroid %in% ortholog.list) 

setwd("~/Documents/GitHub/WGCNA_expanded/K.pneumoniae/panoct_output/pangenome/results")

write.fasta(as.list(filter.centroid$chargaff), 
            paste("centroid_",filter.centroid$centroid, sep="") ,file.out=paste("cutoff_", cutoff, sep=""), open='w')

  #filter.centroid <- unlist(filter.centroid)
  # df <- data.frame(t(sapply(filter.centroid,c)))
  # df.1 <- df %>% 
  #   mutate(count = "1") %>% 
  #   pivot_longer(count, names_to = "centroid", values_to = "chargaff")

  # cutoff_37 <- read.fasta(file="/Users/t/Desktop/cutoff_37", seqtype="DNA", as.string=FALSE)
  # head(cutoff_37)
  # 



#---------------------------------------------

  # Determine number of genes per ortholog
  # counter = 1
  # num.genes <- c()
  # while(counter <= dim(match.table)[1])
  # {
  #   col.num = 2
  #   valid.num = 0
  #   while(col.num <= 75)
  #   {
  #     if(as.vector(match.table[counter,col.num]) != "----------")
  #     {
  #       valid.num = valid.num+1
  #     }
  #     col.num <- col.num+1
  #   }
  #  num.genes <- append(num.genes, valid.num)
  #  counter <- counter+1
  # }
  # match.table <- cbind(match.table, num.genes)
  # names(match.table) [76] <- "total"
  # 
  # # Filter poor matches
  # cutoff = 0 
  # cutoff.vect <- c()
  # sum.vect <- c()
  # while(cutoff < 75)
  # {
  #   sum.vect <- append(sum.vect,match.table %>% 
  #     filter(total > cutoff) %>% 
  #     nrow())
  #   cutoff.vect <- append(cutoff.vect, cutoff)
  #   cutoff <- cutoff +1
  # }
  # filter.ortho <- cbind(cutoff.vect, sum.vect)
  # colnames(filter.ortho) <- c("cutoff.number", "num.orthologs")
  # plot(filter.ortho, main = "Remaining orthologs by cutoff number", pch = 1, cex = 0.4)
  # filter.ortho
  
  
  
  
  
  
  
