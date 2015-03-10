#!/usr/bin/env Rscript
## Given a file (with the 1st column as the read length and 2nd column as the number of reads), this script gives out the plot of length distribution

library(ggplot2)

args <- commandArgs(TRUE)
lendis <- args[1]
output <- args[2]

a <- read.table(lendis)
total <- sum(a$V2)

my.total <- sum(a$V2)
my.median <- round( median(rep(a$V1, times=a$V2)), 2)
my.mean <- round( mean(rep(a$V1, times=a$V2)), 2)

pdf(output)
ggplot(a, aes(x=V1, y=V2)) + geom_bar(stat="identity") + xlab("Read length") + ylab("Number of reads") + ggtitle(paste("Read length distribution", "\nTotal: ", my.total, " Median: ", my.median, " Mean: ", my.mean, sep=""))
dev.off()

