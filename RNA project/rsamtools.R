library(Rsamtools)

bampath <- "~/Nur/Traineeship/GSM1381219.bam"
indexBam("~/Nur/Traineeship/GSM1381219.bam")
scanBam(bampath)
countBam(bampath)
quickBamFlagSummary(bamFile)

