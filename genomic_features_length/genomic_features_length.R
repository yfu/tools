pdf("genomic_features_length.pdf")
feature_fn <- "/data/fuy2/piPipes/common/mm9/UCSC.refSeq.Exon.bed.gz"
feature <- read.table(feature_fn)

len <- feature$V3 - feature$V2
len <- len[len < 1000]
hist(len, breaks=1000, xlim=c(0, 1000), main="length dist of Exons")

for (i in c("3UTR", "5UTR", "Intron")) {
    feature_fn <- paste("/data/fuy2/piPipes/common/mm9/UCSC.refSeq.", i, ".bed.gz", sep="")
    ## UCSC.refSeq.3UTR.bed.gz
    feature <- read.table(feature_fn)
    len_raw <- feature$V3 - feature$V2
    lim=10000
    len <- len_raw[len_raw < lim]
    hist(len, breaks=1000, xlim=c(0, lim), main=paste("length dist of ", i, sep=""))

    lim=1000
    len <- len_raw[len_raw < lim]
    hist(len, breaks=1000, xlim=c(0, lim), main=paste("length dist of ", i, sep=""))

}

dev.off()
