#!/bin/bash

## Use -f folder to specify organism of interest
while getopts ":f:" opt; do
  case $opt in
    f) folder="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&1
    ;;
  esac
done


#bad_genomes=("GCA_000219945.2", "GCA_000585875.1", "GCA_000585895.1")
# echo ${bad_genomes[*]}
# wget -O genbank_bacterial_assemblies.txt https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# wget -O refseq_bacterial_assemblies.txt https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
# awk -F '\t' '{print $20}' refseq_bacterial_assemblies.txt > refseq_bac_assemblies.txt
# awk -F '\t' '{print $20}' genbank_bacterial_assemblies.txt > genbank_bac_assemblies.txt
# rm genbank_bacterial_assemblies.txt
# rm refseq_bacterial_assemblies.txt
# mkdir -p $folder/genomes
filename=$folder/genome_accessions.txt
output=$folder/genomes
panoct=$folder/panoct_output

# # Make blank list
# echo "" > $output/gb.list
# # Loop through all SRA IDs
# while read -r line;
# do
#   name="$line"
# 	name=${name%$'\r'}
# 	name=( $name )
#   FTP_base=$(grep $name genbank_bac_assemblies.txt)
# 	FTP_file=$FTP_base/${FTP_base##*/}_genomic.gbff.gz
#   wget -O $output/$name.gbff.gz $FTP_file
#   [ -s $output/$name.gbff.gz ]
#   if [[ $? -gt 0 ]];then
#     # echo $output/$name.gbff.gz
#     rm -rf $output/$name.gbff.gz
#     name="GCF_${name#GCA_}"
#     # echo $name
#     # echo "LOOK HERE ______------______------"
#     FTP_base=$(grep $name refseq_bac_assemblies.txt)
#     # echo $FTP_base
# 	  FTP_file=$FTP_base/${FTP_base##*/}_genomic.gbff.gz
# 	  wget -O $output/$name.gbff.gz $FTP_file
#   fi
#   gunzip -f $output/$name.gbff.gz
#   counterz=$(grep -c locus_tag $output/$name.gbff) # Only count gbff files with locus annotations
#   if [[ "$counterz" -eq 0 ]];then
# 		echo "bad string ${name}"
# 	else
# 		echo "${WORK}/tingle/WGCNA_expanded/${output}/${name}.gbff" >> $output/gb.list
# 	fi
# # This feeds the appropriate file into the while loop
# done < "$filename"

# sed -i '1d' $output/gb.list
mkdir -p $panoct/pangenome

parse_genbank_files.pl -l $output/gb.list -o $panoct --both --no_dos2unix
# check_att_files.pl --genomes_list $panoct/check_att_files.input --att_dir $panoct --fasta_dir $panoct  --resolve

run_panoct.pl -g $panoct/check_att_files.input --att_dir $WORK/tingle/WGCNA_expanded/$panoct --working_dir $panoct/pangenome --no_grid --use_nuc --fasta_dir $panoct
