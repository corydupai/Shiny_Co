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


filename=$folder/Sra_real.txt
mkdir $folder/fastq

# Loop through all SRA IDs
while read -r line; 
do
	
	name="$line"
	# See if SRA ID is for a nanopore read
	file_count=$(grep $name $folder/nanopore.txt)
	size=${#file_count} 
	
	exists=$(ls $folder/fastq | grep $name)
	exists_size=${#exists}
	
	testing=$(echo "$exists" | grep "_1")
        testing1=${#testing}

	t=$(echo "$exists" | grep "_2")
	t2=${#t}

	if [[ $exists_size -eq 0 ]]; then

                # If ID isn't nanopore, use fasterq-dump to get fastq files
                if [[ $size -eq 0 ]]; then
                        echo "Fasterq $name"
                        #fasterq-dump $name -O $folder/fastq -e 8 -f &
			#prefetch $name
			#fastq-dump --split-files $name -O $folder/fastq --skip-technical
                # If ID is nanopore, use fastq-dump to get fastq files
                else
                        echo "WEIRDOOOOO $name"
                        #fastq-dump --table SEQUENCE $name -O $folder/fastq
                fi
#	elif [[ $testing1 -ne 0 ]] && [[ $t2 -eq 0 ]]; then
#		testing=$(echo "$exists" | grep "_1")
#		testing1=${#testing}
#		echo $testing
#		echo $name
	fi

# This feeds the appropriate file into the while loop
done < "$filename"
