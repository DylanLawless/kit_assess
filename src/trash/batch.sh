#!/bin/sh
#SBATCH --job-name=realrecal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --mem=50gb
#SBATCH --out=pipeline.%J.out
#SBATCH --error=pipeline1.%J.err
 
# Load all modules required
module load gcc/6.4.0
module load gcc/7.3.0
module load intel/18.0.2
module load samtools/1.8
module load picard/2.18.14

rawdata=/work/gr-fe/lawless/target_kit_quality_AH_CH/data/raw

LIST1=/work/gr-fe/lawless/script/pri/sample_list/align/done1.txt
FINLIST=/scratch/lawless/pri/parallel/3.sort/complete.txt
DEDUPLIST=/scratch/lawless/pri/parallel/4.dedup/complete.txt
REALTARLIST=/scratch/lawless/pri/parallel/5.realtar/complete.txt
REALNLIST=/scratch/lawless/pri/parallel/6.indelrealn/complete.txt
BASERECALLIST=/scratch/lawless/pri/parallel/7.baserecal/complete.txt
PRINTBAM=/scratch/lawless/pri/parallel/9.printbam/complete.txt
GVCF=/scratch/lawless/pri/parallel/10.gvcf/complete.txt

# Directories
TRIMDIR=/scratch/lawless/pri/parallel
BAMDIR=/scratch/lawless/pri/parallel/3.bam
SORTDIR=/scratch/lawless/pri/parallel/3.sort
DEDUP=/scratch/lawless/pri/parallel/4.dedup
REALTAR=/scratch/lawless/pri/parallel/5.realtar
INDELREALN=/scratch/lawless/pri/parallel/6.indelrealn
BASERECAL=/scratch/lawless/pri/parallel/7.baserecal
PRINTBAM=/scratch/lawless/pri/parallel/9.printbam
GVCF=/scratch/lawless/pri/parallel/10.gvcf

# Tools
bwa=/work/gr-fe/lawless/tool/bwa-0.7.17/bwa

# References
hg1kv37=/work/gr-fe/lawless/ref/b37/human_g1k_v37.fasta
GATK=/work/gr-fe/lawless/tool/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
K1GINDEL=/work/gr-fe/lawless/ref/b37/1000G_phase1.indels.b37.vcf
MILLS1000G=/work/gr-fe/lawless/ref/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
DBSNP=/work/gr-fe/lawless/ref/b37/dbSnp146.b37.vcf.gz
SSV5=/work/gr-fe/lawless/ref/SureSelectAllExonV5/S04380110_Regions_b37.bed

echo started at `date`

# align with bwa
for line in `cat $fastq_list`
do
  $bwa mem \
  -t $SLURM_NTASKS_PER_NODE \
  -M $hg1kv37 \
  $TRIMDIR/2.trimredo/$line.R1_val_1.fq.gz \
  $TRIMDIR/2.trimredo/$line.R2_val_2.fq.gz \
  -v 1 -R '@RG\tID:$line\tSM:$line\tPL:ILLUMINA\tLB:$line_exome' \
  -M | samtools view -Sb - > $BAMDIR/$line.bam 
done

# Sort with Picard Samtools
# for line in `cat $LIST3`
# do
#     picard SortSam \
#     I=$BAMDIR/$line.bam \
#     O=$SORTDIR/$line.sort.bam \
#     SO=coordinate CREATE_INDEX=TRUE &&\
# echo $SORTDIR/$line.sort.bam >> $FINLIST
# done

# for line in `cat $LIST3`
# do
#     picard MarkDuplicates \
#     I=$SORTDIR/$line.sort.bam \
#     O=$DEDUP/$line.sort.deup.bam \
#     M=$DEDUP/$line.sort.deup.metrics CREATE_INDEX=TRUE && \
#     echo $DEDUP/$line.sort.dedup.bam >> $FINLIST
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
for line in `cat $LIST1`
do
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
java -jar $GATK \
-T PrintReads \
    -R $hg1kv37 \
    -I $INDELREALN/$line.sort.deup.indelrealn.bam \
    -BQSR $BASERECAL/$line.sort.deup.indelrealn.recal.grp \
    -o $PRINTBAM/$line.sort.deup.indelrealn.recal.bam \
    --disable_indel_quals && \
    echo $PRINTBAM/$line.sort.deup.indelrealn.recal.grp >> $PRINTBAMLIST
java -jar $GATK \
    -T HaplotypeCaller \
    --emitRefConfidence GVCF \
    -R $hg1kv37 \
    -D $DBSNP \
    -stand_call_conf 30 \
    -I $PRINTBAM/$line.sort.deup.indelrealn.recal.bam \
    -o $GVCF/$line.sort.deup.indelrealn.recal.HC.g.vcf \
    -L $SSV5 -ip 30 && \
    echo $GVCF/$line.HC.gvcf >> $GVCFLIST
done
