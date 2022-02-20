# Comparison of two target sequencing approaches

## About
Author: Dylan Lawless
Date: 20220220

## Task
The complete task should contain the following files
1. readme.txt
2. AH_S1_L001_R1.fastq.gz
3. AH_S1_L001_R2.fastq.gz
4. AH_S1_target.txt
5. CH_S2_L001_R1.fastq.gz
6. CH_S2_L001_R2.fastq.gz
7. CH_S2_target.txt

The two datasets (S1 and S2) of paired-end short reads are from the SAME human DNA NGS library for clinical diagnosis of solid tumor. 
They were obtained through two different target sequencing approaches with corresponding target region file provided (hg19). 
Please evaluate their performances as much as you can and compile your results into a task report to submit back.

1. the task has to be finished in one week by yourself
2. if the two approaches are wrapped into two commercial NGS products, 
which one would you chose in your lab for clinical diagnose? And why?

## Configuration
* Local: macOS v11.6
* Remote: Red Hat Enterprise Linux Server 7.6 (Maipo)
* R: R v4.1.0 Camp Pontanezen
* fastqc: v0.11.7
* compiler intel: (19.0.5 and 19.1.1)
* samtools: v1.10
* bwa: v0.7.17
* picard: v2.20.8
* R libraries: versions unlisted

## Protocol order
1. src/md5sum.ch
2. src/read_count.sh
3. fastQC: manual run all fastq and save to ./processed/fastqc
3. fastqcr: Assess fastqc rerpots further
4. src/target_info.sh
5. src/1.trim.sh
6. src/2.align.sh
7. src/3.sort.sh
8. src/mapping.sh

## md5sum
* description: A record of data integrety and date received.
* script: src/md5sum.ch
* output: data/processed/metadata/raw.md5sum

## Read count
* description: Check data logic.
Count number of reads by reading fastq lines and dividing by 4. 
Non-whole numbers indicate a truncation.
File sizes that differ unexpectedly may indicate truncation of library error.
* script: src/read_count.sh
* output: data/read_counts/read_counts.txt

## FastQC
* description: FastQC was run all fastq with html reports saved to ./processed/fastqc. Cite FastQC: Andrews S (2010).
* software: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.
* output: data/fastqc/*.html 

## FastqcR
* description: The reports generated from fastQC were also assessed with fastqcr.
* software: <https://rpkgs.datanovia.com/fastqcr/index.html>
* script: src/fastqcr.R
* output: data/fastqcr/qc_report.html (and p1-6.png)

## Data exploration
* description: To learn what target regions are captured in this dataset, the coordinates 
(hg19) were annotated using Ensembl Biomart. The data foratting was changed from
"tab" delimiter to ":".
The data was annotated using the web interface at:
<http://grch37.ensembl.org/biomart/martview/3b67c8bb1f8be31e245e76a1b368fc8b>

Method:
* List of coordinates with the format: 1:100:200 (Chr:start:stop).
* Dataset used: Ensembl Genes 105, Human genes (GRCh37.p13).
* Filters: Region -> Multiple regions -> input list of coordinates
* Output Attributes:
* Gene stable ID, Transcript stable ID, Gene name, Gene start (bp), Gene end (bp), Chromosome/scaffold name, Gene % GC content, HGNC symbol, HGNC ID.
* Save output as TSV.
* Output: ./processed/target_info/file.target_info_biomart.csv
* Output: ./processed/target_info/file.target_info_biomartii_simple.csv (no trascript info).
Note that non-coding region are missing as this query is only for gene information.

* software: <http://grch37.ensembl.org/biomart/martview/3b67c8bb1f8be31e245e76a1b368fc8b>
* script: src/target_info.sh (prep only, manual completion)
* output: data/processed/target_info (prep only, manual completion)

## Trim adaptor
* description: According to fastQC, this dataset used Illumina adaptors.
These were trimmed for alignment, almost most downsteam tools do not require pre-trimming.
TrimGalore was used for this process.
<https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>
Trim Galore also requires installaion of cutadapt <https://github.com/marcelm/cutadapt/>.

* script: src/1.trim
* log:  ../log/data/processed/ (not pulled)
* output: data/processed/1.trim

Make Plots: data/processed/1.trim/CH_S2_L001_R1.fastq.gz_trimming_report.txt

## Align
* description: data was aliged to human_g1k_v37.fasta using BWA MEM. Note that other references may be better suited to assess the method if more information were available.
* script: src/2.align
* log:  ../log/data/processed/ (not pulled)
* output: data/processed/2.align

## Sort
* description: Aligned reads are sorted to allow for downstream processing.
* script: src/3.sort
* log:  ../log/data/processed/ (not pulled)
* output: data/processed/2.sort

## Alignment mapping
samtools flagstat: get mapping summary <https://genomics.sschmeier.com/ngs-mapping/>¬
samtools depth:  read depth for at all positions of the reference genome, e.g. how many reads are overlapping the genomic position.
qualimap:  QualiMap examines sequencing alignment data in SAM/BAM files according to the features of the mapped reads and provides an overall view of the data that helps to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.

* script: src/mapping.sh
* log:  ../log/data/processed/ (not pulled)
* output: data/processed/mapping (mapping, depth, stats).

## Coverage
coverage_plot.R

## Notes
Deduplication was not performed as we are only looking at the read data for now.

## References
Andrews S (2010)."FastQC: a quality control tool for high throughput sequence data". Available online at:
<http://www.bioinformatics.babraham.ac.uk/projects/fastqc>
