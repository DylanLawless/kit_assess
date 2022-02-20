#!/bin/bash
# Dylan Lawless
# 20220219

cd ../data/raw

# This command does:
# print filename then read file content
# count number of lines
# divide number of fastq lines by 4 to give number of reads
# output read_counts to file

printf "counting reads...\n"

for file in *.fastq.gz
do 
	printf "$file " & gzcat $file | \
		wc -l  | \
		awk '{x=$1/4; print x}' 
done > ../processed/read_counts/read_counts.txt

printf "read counts output to data/processed/read_counts/read_counts.txt\n"
