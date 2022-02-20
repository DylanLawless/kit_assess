#!/bin/bash -l
#SBATCH --job-name=sort
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=04:00:00
#SBATCH -o ../log/data/processed/sort.%J.out
#SBATCH -e ../log/data/processed/sort.%J.err

set -e
echo STARTING AT `date`

# Load all modules required
module load intel/19.0.5
module load samtools/1.10
module load parallel
module load picard/2.20.8

# Software
# picard as module

# Directories
rawdata=/work/gr-fe/lawless/target_kit_quality_AH_CH/data/raw
BASEDIR=/work/gr-fe/lawless/target_kit_quality_AH_CH/data/processed
filelist=$BASEDIR/file_list/fastq_list
finlist=$BASEDIR/file_list/finished_list
BAMDIR=$BASEDIR/2.bam
SORTDIR=$BASEDIR/3.sort


# Sort with Picard Samtools
for line in `cat $filelist`
do
    picard SortSam \
    I=$BAMDIR/$line.bam \
    O=$SORTDIR/$line.sort.bam \
    SO=coordinate CREATE_INDEX=TRUE &&\
echo $SORTDIR/$line.sort.bam >> $finlist
done
