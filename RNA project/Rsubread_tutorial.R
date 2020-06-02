### R code from vignette source 'Rsubread.Rnw'

###################################################
### code chunk number 1: Rsubread.Rnw:70-73
###################################################
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)


###################################################
### code chunk number 2: Rsubread.Rnw:90-92
###################################################
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align.stat <- align(index="reference_index",readfile1=reads,output_file="alignResults.BAM",phredOffset=64)


###################################################
### code chunk number 3: Rsubread.Rnw:99-103
###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")
align.stat2 <- align(index="reference_index",readfile1=reads1,readfile2=reads2,
output_file="alignResultsPE.BAM",phredOffset=64)


###################################################
### code chunk number 4: Rsubread.Rnw:125-135
###################################################
ann <- data.frame(
GeneID=c("gene1","gene1","gene2","gene2"),
Chr="chr_dummy",
Start=c(100,1000,3000,5000),
End=c(500,1800,4000,5500),
Strand=c("+","+","-","-"),
stringsAsFactors=FALSE)
ann
fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE


###################################################
### code chunk number 5: Rsubread.Rnw:142-144
###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",annot.ext=ann,isPairedEnd=TRUE)
fc_PE


###################################################
### code chunk number 6: Rsubread.Rnw:164-166
###################################################
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]


###################################################
### code chunk number 7: Rsubread.Rnw:191-192
###################################################
propmapped("alignResults.BAM")


##################################################
# tutorial

#Summarize single-end reads using built-in RefSeq annotation for mouse genome 'mm10' ('mm10'
#is the default inbuilt genome annotation):
  
featureCounts(files="mapping_results_SE.sam")

#Summarize single-end reads using a user-provided GTF annotation file:
  
featureCounts(files="mapping_results_SE.sam",annot.ext="annotation.gtf",
                isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id")

#Summarize single-end reads using 5 threads:
  
featureCounts(files="mapping_results_SE.sam",nthreads=5)

#Summarize BAM format single-end read data:
  
featureCounts(files="mapping_results_SE.bam")

#Summarize multiple libraries at the same time:
  
featureCounts(files=c("mapping_results1.bam","mapping_results2.bam"))

#Summarize paired-end reads and counting fragments (instead of reads):
  
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE)

#Count fragments satisfying the fragment length criteria, eg. [50bp, 600bp]:
  
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,checkFragLength=TRUE,
                minFragLength=50,maxFragLength=600)

#Count fragments which have both ends successfully aligned without considering the fragment
#length constraint:
  
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,requireBothEndsMapped=TRUE)

#Exclude chimeric fragments (exons from different transcripts are mixed together in a fusion product known as chimeric RNA) from fragment counting:
  
featureCounts(files="mapping_results_PE.bam",isPairedEnd=TRUE,countChimericFragments=FALSE)