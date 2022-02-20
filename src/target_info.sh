#!/bin/bash
# Dylan Lawless
# 20220219

set -e

# Use biomart to get info about each target.
# Convert target.txt to the required format (convert "tab" to ":").
# The output can be used on biomart at
# http://grch37.ensembl.org/biomart/martview/

input="../raw"
output="../processed/target_info"

# note the tab delim may not be recognised corretly on different OS
cat $input/AH_S1_target.txt | sed 's/	/:/g' > $output/AH_S1_target_info.txt
cat $input/CH_S2_target.txt | sed 's/	/:/g' > $output/CH_S2_target_info.txt
