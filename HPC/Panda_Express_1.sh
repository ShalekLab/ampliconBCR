#!/bin/bash

#first arg should be the full path to the folder containing the PE fastqs
fastqdir=$1

#if the path is provided, run pandaseq on every fastq file (and its pair) in that folder
if [[ -d $fastqdir ]]; then
    for i1 in $fastqdir/*_R1_001.fastq.gz; do
	i2=${i1/R1_001/R2_001}
	paired=${i1/R1_001.fastq.gz/paired.fastq}
	pandaseq -F -f "$i1" -r "$i2" -w "$paired" -g log.txt -l 200
	echo "$paired is done"
    done
#throwing an error
else	
    echo "User didn't supply full path to directory"

fi

#make an output directory and move the new files there
mkdir paired_fastqs/
mv *paired.fastq paired_fastqs/

