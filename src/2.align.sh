#!/bin/bash -l
#SBATCH --job-name=align
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=03:00:00
#SBATCH -o ../log/data/processed/align.%J.out
#SBATCH -e ../log/data/processed/align.%J.err

set -e
echo STARTING AT `date`

# Load all modules required
module load intel/19.0.5
module load samtools/1.10
module load parallel

# Software
BWA=/work/gr-fe/lawless/tool/bwa-0.7.17/bwa
# samtools as module

# Directories
BASEDIR=/work/gr-fe/lawless/target_kit_quality_AH_CH/data
TRIM=$BASEDIR/processed/1.trim
BAM=$BASEDIR/processed/2.bam
filelist=$BASEDIR/processed/file_list/fastq_list

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta

# -t $SLURM_NTASKS_PER_NODE \
for line in `cat $filelist`
do
	$BWA mem \
	-t 2 \
	-M $hg1kv37 \
	$TRIM/$line\_R1_val_1.fq.gz \
	$TRIM/$line\_R2_val_2.fq.gz \
	-v 1 -R '@RG\tID:$line\tSM:$line\tPL:ILLUMINA\tLB:$line_exome' \
	-M | samtools view -S -b > $BAM/$line.bam 
done

echo FINISHED at `date`
