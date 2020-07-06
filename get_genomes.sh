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


bad_genomes=("GCA_000219945.2", "GCA_000585875.1", "GCA_000585895.1")
echo ${bad_genomes[*]}
# wget -O bacterial_assemblies.txt https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# awk -F '\t' '{print $20}' bacterial_assemblies.txt > bac_assemblies.txt
# rm bacterial_assemblies.txt
mkdir -p $folder/genomes
filename=$folder/genome_accessions.txt
output=$folder/genomes
panoct=$folder/panoct_output
# download_ncbi_annotation.pl -a $folder/genome_accessions.txt -o $output

# echo "" > $output/gb.list
# Loop through all SRA IDs
# while read -r line;
# do
#         
#         name="$line"
# 	name=${name%$'\r'}
# 	#name=${name%$'\s'*}
# 	name=( $name )
#         FTP_base=$(grep $name bac_assemblies.txt)
# 	FTP_file=$FTP_base/${FTP_base##*/}_genomic.gbff.gz
# echo $name
# wget -O $output/$name.gbff.gz $FTP_file
# gunzip -f $output/$name.gbff.gz
# 	if [[ "${bad_genomes[@]}" =~ "${name}" ]];then
# 		echo "bad string ${name}"
# 	else
# 		#echo "${name}"
# 		echo "${WORK}/tingle/WGCNA_expanded/${output}/${name}.gbff" >> $output/gb.list
# 	fi
# # This feeds the appropriate file into the while loop
# done < "$filename"

# sed -i '1d' $output/gb.list
# mkdir -p $panoct/pangenome

parse_genbank_files.pl -l $output/gb.list -o $panoct --both --no_dos2unix

run_panoct.pl -g $panoct/check_att_files.input  --att_dir $WORK/tingle/WGCNA_expanded/$panoct --working_dir $panoct/pangenome --no_grid --use_nuc --fasta_dir $panoct
