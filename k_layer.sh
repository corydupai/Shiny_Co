#!/bin/bash
# Get input arguments
folder=$1
name=$2
size=$3
DIR=$folder/Kallisto_results

# Run Kallisto with different settings for paired end or single reads
	if [[ $size -eq 0 ]]; then
		echo $folder/fastq/${name}.fastq.gz
		kallisto quant -i  $DIR/index.idx --single -o $DIR/$name \
		-l 75 -t $thread -s 1 $folder/fastq/${name}.fastq.gz \
		> $DIR/stdout/${name}.txt
	else
		echo $folder/fastq/${name}_1.fastq.gz
		kallisto quant -i  $DIR/index.idx -o $DIR/$name -t $thread \
		$folder/fastq/${name}_1.fastq.gz $folder/fastq/${name}_2.fastq.gz \
		> $DIR/stdout/${name}.txt
	fi