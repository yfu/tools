#!/usr/bin/env Rscript
## The first column is the coverage by lariat supporting reads (i.e. those reads crossing BP with at leats 5nt on the left and the right)
## the second column is the # sequencing error

library(ggplot2)

args <- commandArgs(TRUE)
e <- args[1]
output <- args[2]

a <- read.table(e)
colnames(a) <- c("pos.on.lar", "n.covered", "n.err")
a$err.rate <- a[, 3] / a[, 2]
pdf(output)
plot(a$pos.on.lar, a$err.rate, xlab="position on lariat", ylab="error rate", main=paste("Error rates calculated from lariat supporting reads", "\n", e, sep=""), pch=16)
dev.off()
