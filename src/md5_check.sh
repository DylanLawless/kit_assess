#!/bin/bash
# Dylan Lawless
# 20220219

set -e

output="../processed/metadata/raw.md5sum"

# record of raw data received
printf "User: Lawless\n" > $output

date >> $output

find ../raw -type f -exec md5sum {} \; >> $output


