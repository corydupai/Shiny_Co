# WGCNA_expanded
Code to perform WGCNA analyses on bacterial datasets from SRA.

Basic pipeline will be as follows:

1. Search SRA using the Organism and Strategy options and download runinfo and accession list from SRA.
2. Add files along with list of nanopore extensions to a folder in the WGCNA_expanded repository. Also update the notes section to include sample collection date/time.
3. Run get_fastq.bash with your folder specified via the -f option.
4. Make transcriptome for your organism (need to plan how to do this).
5. Run Kallisto on your reads.
6. R script to normalize reads and do WGCNA (code not included yet).

