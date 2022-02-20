# Dylan Lawless
# 20220219

# <https://rpkgs.datanovia.com/fastqcr/index.html>
# install.packages("fastqcr")
library(fastqcr)

# First, run fastQC and save the report (zip, html, etc.) The zip is then read by fastqcr
# You can install fastqc from R and run it here also, if you prefer. 

# add the fastqc report directory
qc.dir <- "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqc"
list.files(qc.dir)
qc <- qc_aggregate(qc.dir)
qc

library(dplyr)
qc %>%
	select(sample, module, status) %>%    
	filter(status %in% c("WARN", "FAIL")) %>%
	arrange(sample)

summary(qc)

qc_stats(qc)

# Build a report
qc_report(qc.dir, result.file = "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqcr/qc_report",
			 experiment = "Comparison of AH_S1 and CH_S2",
			 template = NULL,
			 interpret = TRUE)

# Plot individual files.
# I don't think you can get plots for multi samples together
qc.file1 <- "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqc/AH_S1_L001_R1_fastqc.zip"
qc.file2 <- "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqc/AH_S1_L001_R2_fastqc.zip"
qc.file3 <- "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqc/CH_S2_L001_R1_fastqc.zip"
qc.file4 <- "~/Desktop/target_kit_quality_AH_CH/data/processed/fastqc/CH_S2_L001_R2_fastqc.zip"

# interpret
qc <- qc_read(qc.file1)

p1 <- qc_plot(qc, "Per base sequence quality")
p2 <- qc_plot(qc, "Per sequence quality scores")
p3 <- qc_plot(qc, "Per base sequence content")
p4 <- qc_plot(qc, "Per sequence GC content")
p5 <- qc_plot(qc, "Sequence duplication levels")
p6 <- qc_plot(qc, "Adapter content")

qc_plot(qc, "Overrepresented sequences")
qc_plot(qc, "Per base N content")
qc_plot(qc, "Sequence length distribution")
qc_plot(qc, "Kmer content")

require(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2, left = "", bottom = "")

df.lst <- list(qc.file1, qc.file2, qc.file3, qc.file4)

plotdata <- function(x) { qc_plot(qc_read(x), "Per base sequence quality") }
png("../data/processed/fastqcr/p1.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()

plotdata <- function(x) { qc_plot(qc_read(x), "Per sequence quality scores") }
png("../data/processed/fastqcr/p2.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()

plotdata <- function(x) { qc_plot(qc_read(x), "Per base sequence content") }
png("../data/processed/fastqcr/p3.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()

plotdata <- function(x) { qc_plot(qc_read(x), "Per sequence GC content") }
png("../data/processed/fastqcr/p4.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()

plotdata <- function(x) { qc_plot(qc_read(x), "Sequence duplication levels") }
png("../data/processed/fastqcr/p5.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()

plotdata <- function(x) { qc_plot(qc_read(x), "Adapter content") }
png("../data/processed/fastqcr/p6.png", width = 1200, height = 600)
do.call(grid.arrange,lapply(df.lst, plotdata))
dev.off()






