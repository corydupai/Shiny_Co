# Shiny Co
This repository houses an analysis pipeline and Shiny application that processes & visualizes Coexpression Network Data for some of the most highly studied bacterial species.

To process your data:
- Download genomes for species of interest and make a pangenome.
- Download RNAseq fastqs for samples of interest and map to pangenome via Kallisto.
- Run RNAseq results through DESeq2 normalization script (gets rid of genes that are uncommon and samples with low read counts or low percentage of reads mapping to pangenome).
- Run normalized results through WGCNA pipeline to create gene coexpresson networks.

To run the ShinyCo app:
- Put data in correct place (will clarify)
- Run app, select data, and have fun :)

