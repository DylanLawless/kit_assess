#!/bin/sh
#SBATCH --job-name=Bwa_Mem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=18:00:00
#SBATCH --mem=50gb
#SBATCH --output=BwaMem.%J.out
#SBATCH --error=BwaMem.%J.err
 
# Load all modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2
module load samtools/1.8

LIST1=/work/gr-fe/lawless/script/pri/sample_list/align/done1.txt
LIST2=/work/gr-fe/lawless/script/pri/sample_list/align/done2.txt
LIST3=/work/gr-fe/lawless/script/pri/sample_list/align/done3.txt
LIST4=/work/gr-fe/lawless/script/pri/sample_list/align/done4.txt
LIST5=/work/gr-fe/lawless/script/pri/sample_list/align/done5.txt
LIST6=/work/gr-fe/lawless/script/pri/sample_list/align/done6.txt
LIST7=/work/gr-fe/lawless/script/pri/sample_list/align/done7.txt
LIST8=/work/gr-fe/lawless/script/pri/sample_list/align/done8.txt
LIST9=/work/gr-fe/lawless/script/pri/sample_list/align/done9.txt
LIST10=/work/gr-fe/lawless/script/pri/sample_list/align/done10.txt
LIST11=/work/gr-fe/lawless/script/pri/sample_list/align/done11.txt
LIST12=/work/gr-fe/lawless/script/pri/sample_list/align/done12.txt

# Directories
TRIMDIR=/scratch/lawless/pri/parallel
BAMDIR=/scratch/lawless/pri/parallel/3.bam
SORTDIR=/scratch/lawless/pri/parallel/3.sort

# Tools
bwa=/work/gr-fe/lawless/tool/bwa-0.7.17/bwa

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta


# align with bwa
for line in `cat $LIST11`
do
$bwa mem \
-t $SLURM_NTASKS_PER_NODE \
-M $hg1kv37 \
$TRIMDIR/2.trimredo/$line.R1_val_1.fq.gz \
$TRIMDIR/2.trimredo/$line.R2_val_2.fq.gz \
-v 1 -R '@RG\tID:$line\tSM:$line\tPL:ILLUMINA\tLB:$line_exome' \
-M | samtools view -Sb - > $BAMDIR/$line.bam 
done
