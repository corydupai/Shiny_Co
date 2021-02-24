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

counter=1
# echo $folder
cd $folder/fastq
for file in *.fastq
do
	# remainder=$((counter % 1))
	echo $file
	gzip -f -1 $file
	# if [ "$remainder" -eq "0" ]; then
	# 	wait
	# fi
done