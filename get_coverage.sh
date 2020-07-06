#!/bin/bash
# Purpose: Determine sample coverage (pseuodaligned reads v total reads)

# Specify stdout folder from Kallisto_results
rm "/stor/work/Wilke/tingle/WGCNA_expanded/coverage.tsv"
COUNT=0
RUN=""
COVERAGE=""
FILES=/stor/work/Wilke/tingle/WGCNA_expanded/K.pneumoniae/Kallisto_results/stdout/*
for f in $FILES
do
  RUN=$(basename $f .txt)
  COVERAGE=$(grep "processed" $f)
  echo "$RUN  $COVERAGE" >> /stor/work/Wilke/tingle/WGCNA_expanded/coverage.tsv
  #echo $RUN
  #echo grep processed $f 
  #let "COUNT++"
done #> /stor/work/Wilke/tingle/WGCNA_expanded/coverage.csv

#echo $COUNT
