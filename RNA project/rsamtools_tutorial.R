if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")
BiocManager::install("GenomicFeatures")
BiocManager::install("RNAseqData.HNRNPC.bam.chr14")

library(Rsamtools)

# The names of the aln list are basically the names used in the BAM specification. Here is a quick list of some important ones
# qname: The name of the read.
# rname: The name of the chromosome / sequence / contig it was aligned to.
# strand: The strand of the alignment.
# pos: The coordinate of the left-most part of the alignment.
# qwidth: The length of the read.
# mapq: The mapping quality of the alignment.
# seq: The actual sequence of the alignment.
# qual: The quality string of the alignment.
# cigar: The CIGAR string (below).
# flag: The flag (below).
# 
# The CIGAR string
# The "CIGAR" is how the BAM format represents spliced alignments. For example, the format stored 
# the left most coordinate of the alignment. To get to the right coordinate, you have to parse the 
# CIGAR string. In this example "36M" means that it has been aligned with no insertions or deletions. 
# If you need to work with spliced alignments or alignments containing insertions or deletions, you 
# should use the GenomicAlignments package. 
# 
# An alignment may have a number of "flags" set or unset. These flags provide information about the 
# alignment. The flag integer is a representation of multiple flags simultanously. An example of a 
# flag is indicating (for a paired end alignment) whether both pairs have been properly aligned. 
# For more information, see the BAM specification.

which <- IRangesList(seq1 = IRanges(1000, 2000), seq2 = IRanges(c(100, 1000), c(1000, 2000)), compress=FALSE)
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)

bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
bam <- scanBam(bamFile, param=param)
class(bam)
names(bam)

#Each element is itself a list, containing the elements specified by the what and
#tag arguments to ScanBamParam.
class(bam[[1]])
#[1] "list"
names(bam[[1]])
#[1] "rname"  "strand" "pos"    "qwidth" "seq"

#The elements are either basic R or IRanges data types
sapply(bam[[1]], class)
#rname         strand            pos         qwidth            seq 
#"factor"       "factor"      "integer"      "integer" "DNAStringSet"


#A paradigm for collapsing the list-of-lists into a single list is
lst <- lapply(names(bam[[1]]), function(elt) {
  do.call(c, unname(lapply(bam, "[[", elt))) })
names(lst) <- names(bam[[1]])
#This might be further transformed, e.g., to a DataFrame (a data structure defined in the IRanges package), with
df <- do.call("DataFrame", lst)
head(df)
#DataFrame with 6 rows and 5 columns
#rname    strand       pos    qwidth                     seq
#<integer> <integer> <integer> <integer>          <DNAStringSet>
#  1         1         1       970        35 TATTAGGAAA...ACTATGAAGA
#2         1         1       971        35 ATTAGGAAAT...CTATGAAGAG
#3         1         1       972        35 TTAGGAAATG...TATGAAGAGA
#4         1         1       973        35 TAGGAAATGC...ATGAAGAGAC
#5         1         1       974        35 AGGAAATGCT...TGAAGAGACT
#6         1         2       975        35 GGAAATGCTT...GAAGAGACTA

# na19240url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/NA19240/alignment/NA19240.chrom6.SLX.maq.SRP000032.2009_07.bam"
# which <- GRanges(seqnames = "6", ranges = IRanges(1e+05, 110000))
# param <- ScanBamParam(which = which)
# na19240bam <- scanBam(na19240url, param = param)

###############
# http://www.bioconductor.org/packages/2.13/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
# na19240bam = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/NA19240/alignment/NA19240.chrom6.SLX.maq.SRP000032.2009_07.bam"
# which <- RangesList("6"=IRanges(100000L, 110000L))
# param <- ScanBamParam(which=which, what=scanBamWhat())
# na19240bam <- scanBam(na19240url, param=param)

#The BAM file in the previous example includes an index, represented by a separate file with extension .bai:
list.files(dirname(bamFile), pattern = "ex1.bam(.bai)?")
#[1] "ex1.bam"     "ex1.bam.bai"

# library(GenomicFeatures)
# library(featuretoolsR)
# bamRanges <- local({fl <- system.file("extdata", "CaffeineTxdb.sqlite", package = "Rsamtools") 
# transcripts(loadFeatures(fl))
# #I don't know why but loadFeatures doesn't work
# })
# slxMaq09 <- local({fl <- system.file("extdata", "slxMaq09_urls.txt", package = "Rsamtools")
# readLines(fl)})
# 
# bamExperiment <- list(description="Caffeine metabolism views on 1000 genomes samples", created=date())
# 
# bv <- BamViews(slxMaq09, bamRanges=bamRanges, bamExperiment=bamExperiment)
# metadata(bamSamples(bv)) <- list(description="Solexa/MAQ samples, August 2009", created="Thu Mar 25 14:08:42 2010")

#####################################
#This is from the second url (https://kasperdanielhansen.github.io/genbioconductor/html/Rsamtools.html)
seqinfo(bamFile)
aln <- scanBam(bamFile)
length(aln)
aln <- aln[[1]]
names(aln)
lapply(aln, function(xx) xx[1])

#BAM files can be extremely big and it is there often necessary to read in parts of the file. 
#You can do this in different ways
# Read a set number of records (alignments).
# Only read alignments satisfying certain criteria.
# Only read alignments in certain genomic regions.
# Let us start with the first of this. By specifying yieldSize when you use BamFile(), 
#every invocation of scanBam() will only read yieldSize number of alignments. 
#You can then invoke scanBam() again to get the next set of alignments; 
#this requires you to open() the file first (otherwise you will keep read the same alignments).

yieldSize(bamFile) <- 1
open(bamFile)
scanBam(bamFile)[[1]]$seq
close(bamFile)
yieldSize(bamFile) <- NA

quickBamFlagSummary(bamFile)
countBam(bamFile)


BamFile(bamFile)
seqinfo(bamFile)

##############################


##https://bioconductor.org/packages/devel/bioc/manuals/Rsamtools/man/Rsamtools.pdf
## BamFile options.
##
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
bf <- BamFile(fl)
bf
## When 'asMates=TRUE' scanBam() reads the data in as
## pairs. See 'asMates' above for details of the pairing
## algorithm.
asMates(bf) <- TRUE
## When 'yieldSize' is set, scanBam() will iterate
## through the file in chunks.
yieldSize(bf) <- 500
## Some applications append a filename (e.g., NCBI Sequence Read
## Archive (SRA) toolkit) or allele identifier to the sequence qname.
## This may result in a unique qname for each record which presents a
## problem when mating paired-end reads (identical qnames is one
## criteria for paired-end mating). 'qnamePrefixEnd' and
## 'qnameSuffixStart' can be used to trim an unwanted prefix or suffix.
qnamePrefixEnd(bf) <- "/"
qnameSuffixStart(bf) <- "."
##
## Reading Bam files.
##
fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)
(bf <- BamFile(fl))
head(seqlengths(bf)) # sequences and lengths in BAM file

if (require(RNAseqData.HNRNPC.bam.chr14)) {
  bfl <- BamFileList(RNAseqData.HNRNPC.bam.chr14_BAMFILES)
  bfl
  bfl[1:2] # subset
  bfl[[1]] # select first element -- BamFile
  ## merged across BAM files
  seqinfo(bfl)
  head(seqlengths(bfl))
}
length(scanBam(fl)[[1]][[1]]) # all records
bf <- open(BamFile(fl)) # implicit index
bf
identical(scanBam(bf), scanBam(fl))
close(bf)
## Use 'yieldSize' to iterate through a file in chunks.
bf <- open(BamFile(fl, yieldSize=1000))
while (nrec <- length(scanBam(bf)[[1]][[1]]))
  cat("records:", nrec, "\n")
close(bf)
## Repeatedly visit multiple ranges in the BamFile.
rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
bf <- open(BamFile(fl))
sapply(seq_len(length(rng)), function(i, bamFile, rng) {
  param <- ScanBamParam(which=rng[i], what="seq")
  bam <- scanBam(bamFile, param=param)[[1]]
  alphabetFrequency(bam[["seq"]], baseOnly=TRUE, collapse=TRUE)
}, bf, rng)
close(bf)



fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)
##
## scanBam
##
res0 <- scanBam(fl)[[1]] # always list-of-lists
names(res0)
length(res0[["qname"]])
lapply(res0, head, 3)
table(width(res0[["seq"]])) # query widths
table(res0[["qwidth"]], useNA="always") # query widths derived from cigar
table(res0[["cigar"]], useNA="always")
table(res0[["strand"]], useNA="always")
table(res0[["flag"]], useNA="always")
which <- IRangesList(seq1=IRanges(1000, 2000),
                     seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which, what=scanBamWhat())
res1 <- scanBam(fl, param=p1)
names(res1)
names(res1[[2]])
p2 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
res2 <- scanBam(fl, param=p2)
p3 <- ScanBamParam(
  what="flag", # information to query from BAM file
  flag=scanBamFlag(isMinusStrand=FALSE))
length(scanBam(fl, param=p3)[[1]]$flag)
##
## idxstatsBam
##
idxstatsBam(fl)
##
## filterBam
##
param <- ScanBamParam(
  flag=scanBamFlag(isUnmappedQuery=FALSE),
  what="seq")
dest <- filterBam(fl, tempfile(), param=param)
countBam(dest) ## 3271 records
## filter to a single file
filter <- FilterRules(list(MinWidth = function(x) width(x$seq) > 35))
dest <- filterBam(fl, tempfile(), param=param, filter=filter)
countBam(dest) ## 398 records
res3 <- scanBam(dest, param=ScanBamParam(what="seq"))[[1]]
table(width(res3$seq))
## filter 1 file to 2 destinations
filters <- list(
  FilterRules(list(long=function(x) width(x$seq) > 35)),
  FilterRules(list(short=function(x) width(x$seq) <= 35))
)
destinations <- replicate(2, tempfile())
dest <- filterBam(fl, destinations, param=param, filter=filters)
lapply(dest, countBam)
##
## sortBam
##
sorted <- sortBam(fl, tempfile())
##
## scanBamParam re-orders 'which'; recover original order
##
gwhich <- as(which, "GRanges")[c(2, 1, 3)] # example data
cnt <- countBam(fl, param=ScanBamParam(which=gwhich))
reorderIdx <- unlist(split(seq_along(gwhich), seqnames(gwhich)))
cnt
cnt[reorderIdx,]



fls <- system.file("extdata", "ex1.bam", package="Rsamtools",
                   mustWork=TRUE)
rngs <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(9, 9)),
                ranges = c(IRanges(seq(10000, 90000, 10000), width=500),
                           IRanges(seq(100000, 900000, 100000), width=5000)),
                Count = seq_len(18L))
v <- BamViews(fls, bamRanges=rngs)
v
v[1:5,]
bamRanges(v[c(1:5, 11:15),])
bamDirname(v) <- getwd()
v






#https://rdrr.io/bioc/Rsamtools/man/scanBam.html
fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)

##
## scanBam
##

res0 <- scanBam(fl)[[1]] # always list-of-lists
names(res0)
length(res0[["qname"]])
lapply(res0, head, 3)
table(width(res0[["seq"]])) # query widths
table(res0[["qwidth"]], useNA="always") # query widths derived from cigar
table(res0[["cigar"]], useNA="always")
table(res0[["strand"]], useNA="always")
table(res0[["flag"]], useNA="always")

which <- IRangesList(seq1=IRanges(1000, 2000),
                     seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which, what=scanBamWhat())
res1 <- scanBam(fl, param=p1)
names(res1)
names(res1[[2]])

p2 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
res2 <- scanBam(fl, param=p2)

p3 <- ScanBamParam(
  what="flag",           # information to query from BAM file 
  flag=scanBamFlag(isMinusStrand=FALSE))
length(scanBam(fl, param=p3)[[1]]$flag)

##
## idxstatsBam
##

idxstatsBam(fl)

##
## filterBam
##

param <- ScanBamParam(
  flag=scanBamFlag(isUnmappedQuery=FALSE),
  what="seq")
dest <- filterBam(fl, tempfile(), param=param)
countBam(dest)  ## 3271 records

## filter to a single file
filter <- FilterRules(list(MinWidth = function(x) width(x$seq) > 35))
dest <- filterBam(fl, tempfile(), param=param, filter=filter)
countBam(dest)  ## 398 records
res3 <- scanBam(dest, param=ScanBamParam(what="seq"))[[1]]
table(width(res3$seq))

## filter 1 file to 2 destinations
filters <- list(
  FilterRules(list(long=function(x) width(x$seq) > 35)),
  FilterRules(list(short=function(x) width(x$seq) <= 35))
)
destinations <- replicate(2, tempfile())
dest <- filterBam(fl, destinations, param=param, filter=filters)
lapply(dest, countBam)

##
## sortBam
##

sorted <- sortBam(fl, tempfile())

##
## scanBamParam re-orders 'which'; recover original order
##

gwhich <- as(which, "GRanges")[c(2, 1, 3)]    # example data
cnt <- countBam(fl, param=ScanBamParam(which=gwhich))
reorderIdx <- unlist(split(seq_along(gwhich), seqnames(gwhich)))
cnt
cnt[reorderIdx,]





##############################
######## URL's ###############
##############################

#http://bioconductor.org/help/course-materials/2010/BioC2010/Exercises-Rsamtools.pdf
#https://kasperdanielhansen.github.io/genbioconductor/html/Rsamtools.html
#http://biobits.org/samtools_primer.html (picture ngs_overview.png)
#http://www.bioconductor.org/packages/2.13/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
#https://bioconductor.org/packages/devel/bioc/manuals/Rsamtools/man/Rsamtools.pdf
