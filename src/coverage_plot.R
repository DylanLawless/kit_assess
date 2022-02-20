# code produced from using Dave Tang's blog
# <https://davetang.org/muse/2013/09/07/creating-a-coverage-plot-in-r/>

#the smallest CAGE BAM file from ENCODE
setwd("./")
system("wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/wgEncodeRikenCageHchCellPapAlnRep1.bam")
# First let’s read in a BAM file using Bioconductor’s Rsamtools and convert it into a data frame (since I’m more familiar with data frames):

#install if necessary
##source("http://bioconductor.org/biocLite.R")
##biocLite("Rsamtools")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rsamtools")


#load library
library(Rsamtools)

#read in entire BAM file
# bam <- scanBam("wgEncodeRikenCageHchCellPapAlnRep1.bam")
bam <- scanBam("../data/processed/3.sort/AH_S1_L001.sort.bam")

#names of the BAM fields
names(bam[[1]])
# [1] "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar"
# [9] "mrnm"   "mpos"   "isize"  "seq"    "qual"

# Col	Field	Description		
# 1	QNAME	Query (pair) NAME	
# 2	FLAG	bitwise FLAG	
# 3	RNAME	Reference sequence NAME	
# 4	POS	1-based leftmost POSition/coordinate of clipped sequence	
# 5	MAPQ	MAPping Quality (Phred-scaled)	
# 6	CIAGR	extended CIGAR string	
# 7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)	
# 8	MPOS	1-based Mate POSition	
# 9	ISIZE	Inferred insert SIZE	
# 10	SEQ	query SEQuence on the same strand as the reference	
# 11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)	
# 12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE	

#distribution of BAM flags
table(bam[[1]]$flag)

#      0       4      16 
#1472261  775200 1652949

#function for collapsing the list of lists into a single list
#as per the Rsamtools vignette
.unlist <- function (x){
  ## do.call(c, …) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)
#[1] 3900410      13

# Now that we’ve stored our BAM file into a data.frame we can make our coverage plot.



#use chr22 as an example
#how many entries on the negative strand of chr22?
table(bam_df$rname == '4' & bam_df$flag == 16)

# FALSE    TRUE 
#3875997   24413

#function for checking negative strand
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}

#test neg function with subset of chr22
test <- subset(bam_df, rname == '4')
dim(test)
#[1] 56426    13
table(apply(as.data.frame(test$flag), 1, check_neg))
#number same as above
#FALSE  TRUE 
#32013 24413

#function for checking positive strand
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}

#check pos function
table(apply(as.data.frame(test$flag), 1, check_pos))
#looks OK
#FALSE  TRUE 
#24413 32013

#store the mapped positions on the plus and minus strands
chr22_neg <- bam_df[bam_df$rname == '4' &
                      apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
]
length(chr22_neg)
#[1] 24413
chr22_pos <- bam_df[bam_df$rname == '4' &
                      apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
]
length(chr22_pos)
#[1] 32013

#calculate the densities
chr22_neg_density <- density(chr22_neg)
chr22_pos_density <- density(chr22_pos)

# we have no neg strand, but i want it to work for testing
# set positive to == negative
chr22_neg_density <- chr22_pos_density

#display the negative strand with negative values
chr22_neg_density$y <- chr22_neg_density$y * -1

plot(chr22_pos_density,
     ylim = range(c(chr22_neg_density$y, chr22_pos_density$y)),
     main = "Coverage plot of mapped CAGE reads",
     xlab = "Chromosome 22",
     col = 'blue',
     lwd=2.5)
lines(chr22_neg_density, lwd=2.5, col = 'red')

chr22_neg_density$y



# ----
# If you want to focus on a specific region, follow the steps above until you create the two vectors chr22_neg and chr22_pos. Remember that these vectors store the starting position of all the reads mapped onto chromosome 22. So to focus on a specific region we just need to shrink or subset this vector:
#check our two vectors
length(chr22_neg)
#[1] 24413
length(chr22_pos)
#[1] 32013

#what is actually inside these two vectors?
head(chr22_pos)
[1] 16096446 16124065 16147449 16165745 16339913 16364755

#get a summary of the mapped positions
summary(chr22_pos)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#16100000 31480000 35790000 35040000 39080000 51240000

#say I am interested in region 31480000 to 39080000
lower <- 31480000
upper <- 39080000

chr22_pos_interest <- chr22_pos[chr22_pos > lower & chr22_pos < upper]
#check how many entries we have
length(chr22_pos_interest)
#[1] 16220

#do the same for the negative strand
chr22_neg_interest <- chr22_neg[chr22_neg > lower & chr22_neg < upper]
length(chr22_neg_interest)
#[1] 7744

#now continue with the code above
#but with our two new vectors of interest
chr22_neg_density <- density(chr22_neg_interest)
chr22_pos_density <- density(chr22_pos_interest)

#display the negative strand with negative values
chr22_neg_density$y <- chr22_neg_density$y * -1

plot(chr22_pos_density,
     ylim = range(c(chr22_neg_density$y, chr22_pos_density$y)),
     main = "Coverage plot of mapped CAGE reads",
     xlab = "Chromosome 22",
     col = 'blue',
     lwd=2.5,
     type='h'
)
lines(chr22_neg_density, lwd=2.5, col = 'red', type='h')
