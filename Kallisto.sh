#!/bin/bash

while getopts ":f:" opt; do
  case $opt in
    f) folder="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&1
    ;;
  esac
done

DIR=Kallisto_results
counter=0
thread=1
# Make this directory and any parent directories that are missing
mkdir -p $folder/$DIR/stdout
filename=$folder/Sra_real.txt

kallisto index -i $folder/$DIR/index.idx $folder/filtered_transcripts.fasta

#This loops through a text file with SRA accessions line by line
while read -r line;
do
	name="$line"
	#echo $name
	size=$(ls $folder/fastq | grep -c $name)

	#echo $size
	# Call script that runs kallisto
	bash k_layer.sh $folder $name $size &

	# Run 44 samples in parallel. We can probably up this number.
	counter=$((counter + 1))
	remainder=$((counter % 20))
	if [ "$remainder" -eq "0" ]
	then
		wait
	fi
done < "$filename"
