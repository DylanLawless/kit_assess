#!/bin/bash
#SBATCH --job-name=SamSortParallel1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=30:00:00
#SBATCH --mem=60gb
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
##SBATCH --reservation gr-fe

# Load all modules required
module load gcc/7.3.0
module load intel/18.0.2
module load picard/2.18.14
module load parallel

# Directories
BAMDIR=/scratch/lawless/spss/3.bam/1
SORTDIR=/scratch/lawless/spss/4.sort/1

echo STARTING AT `date`

ls "$BAMDIR"/*.bam |
 xargs -n 1 basename |
  awk -F "." '{print $1 | "sort -u"}' |
   parallel --joblog log1 -j 5 picard SortSam \
    I="$BAMDIR"/{}.bam \
    O="$SORTDIR"/{}.sort.bam \
    TMP_DIR=./tmp \
    SO=coordinate CREATE_INDEX=TRUE &&\
    echo $SORTDIR/{}.sort.bam >> completed_list_sort_b1.txt

echo FINISHED AT `date`
