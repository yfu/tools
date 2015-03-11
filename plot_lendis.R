#!/usr/bin/env Rscript
## Given a file (with the 1st column as the read length and 2nd column as the number of reads), this script gives out the plot of length distribution

library(ggplot2)

args <- commandArgs(TRUE)
lendis <- args[1]
output <- args[2]

## lendis <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/Hi5.nodavirus.ox.lendis"
## output <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/test.pdf"

a <- read.table(lendis)

if( ncol(a) == 2) {
    ## This type of input
    ## 18	351
    ## 19	344
    ## 20	231
    ## 21	210
    write("The table has two columns", stderr())
    my.total <- sum(a$V2)
    my.median <- round( median(rep(a$V1, times=a$V2)), 2)
    my.mean <- round( mean(rep(a$V1, times=a$V2)), 2)
    
    pdf(output)
    ggplot(a, aes(x=V1, y=V2)) + geom_bar(stat="identity") + xlab("Read length") + ylab("Number of reads") + ggtitle(paste("Read length distribution", "\nTotal: ", my.total, " Median: ", my.median, " Mean: ", my.mean, sep=""))
    dev.off()
}

if (ncol(a) == 3) {
    ## This type of input
    ## 18	351	+
    ## 19	344	+
    ## 20	231	+
    ## 21	210	+
    write("The table has three columns. I will use the strand information", stderr())
    idx.plus <- a$V3 == "+"
    idx.minus <- a$V3 == "-"
    plus <- a[idx.plus, ]
    minus <- a[idx.minus, ]
    my.total.plus  <- sum(plus$V2)
    my.total.minus <- sum(minus$V2)
    my.total <- my.total.plus + my.total.minus
    my.median.plus  <- median(rep(plus$V1,  times= plus$V2))
    my.median.minus <- median(rep(minus$V1, times=minus$V2))
    my.median <- median(rep(a$V1, times=a$V2))
    my.mean.plus <- round( mean( rep(plus$V1, times=plus$V2) ), digits=2)
    my.mean.minus <- round( mean( rep(minus$V1, times=minus$V2) ), digits=2)
    my.mean <- mean( rep(a$V1, times=a$V2) )
    
    a[idx.minus, ]$V2 <- -a[idx.minus, ]$V2

    pdf(output)
    my.title <- paste("Read length distribution", "\nTotal: ", my.total, "(", my.total.plus, "+", my.total.minus, ")", " Median: ", my.median, " Mean: ", my.mean, sep="", "\n", "Median(+): ", my.median.plus, " Mean(+): ", my.mean.plus, " Median(-): ", my.median.minus, " Mean(-): ", my.mean.minus)

    minus$V2 <- -minus$V2
    p <- ggplot() + geom_bar(data=minus, aes(x=V1, y=V2, group=V3, fill=V3), stat="identity", position=) + xlab("Read length") + ylab("Number of reads") + ggtitle(my.title) + geom_bar(data=plus, aes(x=V1, y=V2, group=V3, fill=V3), stat="identity", position=)
    print(p)
    dev.off()
}
