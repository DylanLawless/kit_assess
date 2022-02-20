#!/bin/bash -l
#SBATCH --job-name=pipeline
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=04:00:00
#SBATCH -o ../log/data/processed/pipeline.%J.out
#SBATCH -e ../log/data/processed/pipeline.%J.err

# Load all modules required
module load intel/19.0.5
module load samtools/1.10
module load parallel
module load picard/2.20.8

# Directories
rawdata=/work/gr-fe/lawless/target_kit_quality_AH_CH/data/raw
BASEDIR=/work/gr-fe/lawless/target_kit_quality_AH_CH/data/processed
filelist=$BASEDIR/file_list/fastq_list
finlist=$BASEDIR/file_list/finished_list
BAMDIR=$BASEDIR/2.bam
SORTDIR=$BASEDIR/3.sort
DEDUP=$BASEDIR/4.dedup
REALTAR=$BASEDIR/5.realtar
INDELREALN=$BASEDIR/6.indelrealn
BASERECAL=$BASEDIR/7.baserecal
PRINTBAM=$BASEDIR/9.printbam
GVCF=$BASEDIR/10.gvcf

# Tools
bwa=/work/gr-fe/lawless/tool/bwa-0.7.17/bwa

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
K1GINDEL=/work/gr-fe/lawless/ref/b37/1000G_phase1.indels.b37.vcf
MILLS1000G=/work/gr-fe/lawless/ref/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
# SSV5=/work/gr-fe/lawless/ref/SureSelectAllExonV5/S04380110_Regions_b37.bed

echo started at `date`

# # align with bwa
# for line in `cat $filelist`
# do
# 	$BWA mem \
# 	-t 2 \
# 	-M $hg1kv37 \
# 	$TRIM/$line\_R1_val_1.fq.gz \
# 	$TRIM/$line\_R2_val_2.fq.gz \
# 	-v 1 -R '@RG\tID:$line\tSM:$line\tPL:ILLUMINA\tLB:$line_exome' \
# 	-M | samtools view -S -b > $BAMDIR/$line.bam 
# done

# Sort with Picard Samtools
for line in `cat $filelist`
do
    picard SortSam \
    I=$BAMDIR/$line.bam \
    O=$SORTDIR/$line.sort.bam \
    SO=coordinate CREATE_INDEX=TRUE &&\
echo $SORTDIR/$line.sort.bam >> $finlist
done

# for line in `cat $filelist`
# do
#     picard MarkDuplicates \
#     I=$SORTDIR/$line.sort.bam \
#     O=$DEDUP/$line.sort.deup.bam \
#     M=$DEDUP/$line.sort.deup.metrics CREATE_INDEX=TRUE && \
#     echo $DEDUP/$line.sort.dedup.bam >> $finlist
# done

#### 5 Create indel realigner targets
# java -jar $GATK \
# -T RealignerTargetCreator \
#     -R $hg1kv37 \
#     -known $K1GINDEL \
#     -known $MILLS1000G \
#     -I $DEDUP/$line.sort.deup.bam \
#     -o $REALTAR/$line.sort.deup.intervals &&\
#     echo $REALTAR/$line.sort.dedup.intervals >> $REALTARLIST &&\

#### 6 Indel realignment
#### 7 Recalibrate Base Quality Scores (GATK) Get Recalibration Model
#### 8 optional check for base recalibration
#### 9 Print final reads after applying BQSR
#### 10 HC
# for line in `cat $LIST1`
# do
# java -jar $GATK \
# -T IndelRealigner \
#     -R $hg1kv37 \
#     -known $K1GINDEL \
#     -known $MILLS1000G \
#     -I $DEDUP/$line.sort.deup.bam \
#     -targetIntervals $REALTAR/$line.sort.deup.intervals \
#     -o $INDELREALN/$line.sort.deup.indelrealn.bam && \
#     echo $INDELREALN/$line.sort.deup.indelrealn.bam >> $REALNLIST &&\
# java -jar $GATK \
# -T BaseRecalibrator \
#     -R $hg1kv37 \
#     -knownSites $K1GINDEL \
#     -knownSites $MILLS1000G \
#     -knownSites $DBSNP \
#     -o $BASERECAL/$line.sort.deup.indelrealn.recal.grp \
#     -I $INDELREALN/$line.sort.deup.indelrealn.bam \
#     -nct 16 && \
#     echo $BASERECAL/$line.sort.deup.indelrealn.recal.grp >> $BASERECALLIST &&\
# java -jar $GATK \
# -T PrintReads \
#     -R $hg1kv37 \
#     -I $INDELREALN/$line.sort.deup.indelrealn.bam \
#     -BQSR $BASERECAL/$line.sort.deup.indelrealn.recal.grp \
#     -o $PRINTBAM/$line.sort.deup.indelrealn.recal.bam \
#     --disable_indel_quals && \
#     echo $PRINTBAM/$line.sort.deup.indelrealn.recal.grp >> $PRINTBAMLIST
# java -jar $GATK \
#     -T HaplotypeCaller \
#     --emitRefConfidence GVCF \
#     -R $hg1kv37 \
#     -D $DBSNP \
#     -stand_call_conf 30 \
#     -I $PRINTBAM/$line.sort.deup.indelrealn.recal.bam \
#     -o $GVCF/$line.sort.deup.indelrealn.recal.HC.g.vcf \
#     -L $SSV5 -ip 30 && \
#     echo $GVCF/$line.HC.gvcf >> $GVCFLIST
# done
