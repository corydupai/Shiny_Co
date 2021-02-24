#install.packages("genomes")
#Sys.getenv("R_LIBS_USER")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("biomartr")
#install.packages("data.table")
require("biomartr")
require("data.table")
source("source/getGBFF.R")

stuff_in <- fread("../K.pneumoniae/genome_accessions.txt")

getGBFF(db="genbank",stuff_in[1,1],gunzip=T)

#searched <- genomes::esearch( stuff_in[1,1], "nucleotide")
#efetch(searched , db = "nucleotide", rettype = "gb", retmode = "text", showURL = TRUE, "testing.gb)
