#!/bin/bash -l
#SBATCH --job-name=mapping
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=03:00:00
#SBATCH -o ../log/data/processed/mapping.%J.out
#SBATCH -e ../log/data/processed/mapping.%J.err

set -e
echo STARTING AT `date`

# Load all modules required
module load intel/19.0.5
module load samtools/1.10
module load parallel

# software
BWA=/work/gr-fe/lawless/tool/bwa-0.7.17/bwa
# samtools as module
qualimap=/work/gr-fe/lawless/tool/qualimap_v2.2.1/qualimap

# Directories
BASEDIR=/work/gr-fe/lawless/target_kit_quality_AH_CH/data
TRIM=$BASEDIR/processed/1.trim
BAM=$BASEDIR/processed/2.bam
SORT=$BASEDIR/processed/3.sort
MAPPING=$BASEDIR/processed/mapping
filelist=$BASEDIR/processed/file_list/fastq_list

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta

# get mapping summary
# <https://genomics.sschmeier.com/ngs-mapping/>
for line in `cat $filelist`
do
	samtools flagstat $SORT/$line.sort.bam > $MAPPING/$line.sort.mapping
done

# read depth for at all positions of the reference genome, e.g. how many reads are overlapping the genomic position.
for line in `cat $filelist`
do
	samtools depth $SORT/$line.sort.bam | gzip > $MAPPING/$line.sort.depth.txt.gz
done

# need a target list bed file named to match $filelist
# qualimap requires the bed file with 6 columns with NO header included:
# chrom	chromStart	chromEnd	name	score	strand

# cat the file, skip header line, add 3 columns of non-data
cat $BASEDIR/raw/AH_S1_target.txt |\
	tail -n +2 |\
	sed 's/$/	\.	\.	\./' \
	> $BASEDIR/processed/file_list/AH_S1_L001_target.bed

cat $BASEDIR/raw/CH_S2_target.txt |\
	tail -n +2 |\
	sed 's/$/	\.	\.	\./' \
	> $BASEDIR/processed/file_list/CH_S2_L001_target.bed


# QualiMap examines sequencing alignment data in SAM/BAM files according to the features of the mapped reads and provides an overall view of the data that helps to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.
for line in `cat $filelist`
do
	 $qualimap bamqc \
		 --feature-file $BASEDIR/processed/file_list/$line\_target.bed \
		 --paint-chromosome-limits \
		 -bam $SORT/$line.sort.bam
done

# move the output to mapping dir
for line in `cat $filelist`
do
	mv $SORT/$line.sort_stats $MAPPING/
done

echo FINISHED at `date`



