#!/bin/bash -l
#SBATCH --job-name=trim
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 2
#SBATCH --nodes 1
#SBATCH --mem 10GB
#SBATCH --time 02:00:00
#SBATCH -o ../log/data/processed/trim.%J.out
#SBATCH -e ../log/data/processed/trim.%J.err

set -e
echo started at `date`

module load intel/19.1.1
module load parallel
module load fastqc/0.11.7

# software
trim_galore=/work/gr-fe/lawless/tool/TrimGalore-0.5.0/trim_galore
cutadapt=/home/lawless/.local/bin/cutadapt

# data directories
BASEDIR=/work/gr-fe/lawless/target_kit_quality_AH_CH/data
FASTQ=$BASEDIR/raw/*_R1.fastq.gz
TRIM=$BASEDIR/processed/1.trim

printf "\ntrimgalore\n" >> $BASEDIR/processed/log/completion_list

parallel $trim_galore \
--path_to_cutadapt $cutadapt \
-q 20 \
--illumina --paired --gzip \
--fastqc \
-o $TRIM {} {=s/_R1\./_R2\./=} ::: $FASTQ

echo sample_{} {=s/_R1\./_R2\./=} has completed >> $BASEDIR/processed/log/completion_list
