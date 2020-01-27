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
filename=$folder/SraAccList.txt

while read -r line; 
do
	name="$line"
	echo $name
	file_count=$(grep $fname $folder/nanopore.txt)
	size=${#file_count} 
	
	# Call script that runs kallisto
	bash k_layer.sh $folder $name $size &
	
	# Run 44 samples in parallel. We can probably up this number.
	counter=$((counter + 1))
	remainder=$((counter % 44))
	if [ "$remainder" -eq "0" ]
	then
		wait
	fi
done < "$filename"