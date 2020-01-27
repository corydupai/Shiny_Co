#!/bin/bash

# Use -f folder to specify organism of interest
while getopts ":f:" opt; do
  case $opt in
    f) folder="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&1
    ;;
  esac
done


filename=$folder/SraAccList.txt
mkdir $folder/fastq

# Loop through all SRA IDs
while read -r line; 
do
	
	name="$line"
	echo $name
	# See if SRA ID is for a nanopore read
	file_count=$(grep $name $folder/nanopore.txt)
	size=${#file_count} 
	
	# If ID isn't nanopore, use fasterq-dump to get fastq files
	if [[ $size -eq 0 ]]; then
		echo $fnew
		fasterq-dump $name -O $folder/fastq -e 8
	# If ID is nanopore, use fastq-dump to get fastq files
	else
		fastq-dump --table SEQUENCE $name -O $folder/fastq
	fi
# This feads whatever
done < "$filename"