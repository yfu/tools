#!/usr/bin/env Rscript

## This script assume an input file that has three columns, with the 1st
## column specifying the coordinate on the gene, 2nd column specifying
## signals on the plus strand and 3rd column specifying the signals on
## the minus strand

## Usage: plot_gene_signal.R my.4col.depth my.title my.plot.pdf

## Notice that the input table should cover the whole region
plot_gene_signal <- function(a, tt, output) {
    ## write("The input file has 3 columns.", stderr())
    ## adapted from https://github.com/bowhan/piPipes/blob/c98b12324553ca487fb634b584ca213300044903/bin/piPipes.R
    pdf (output)
    par (bty="n")
    # plot (a$V2,a$V3, xlim=c(0,nrow(a)), ylim=c(min(a$V4), max(a$V3)) , type='n', xlab=paste("Gene body", nrow(a), sep=" "), ylab="Signal", tck=0.01, main=tt)
    plot (a$V2,a$V3, xlim=c(0,nrow(a)), ylim=range(pretty(c(a$V4, a$V3))) , type='n', xlab=paste("Gene body", nrow(a), sep=" "), ylab="Signal", tck=0.01, main=tt)    
    points (a$V2, a$V3, col="blue", type="s")
    points (a$V2, a$V4, col="red", type="s")
    abline (h=0, lty=2)
    dev.off()
}

args <- commandArgs(TRUE)
input <- args[1]
my.title <- args[2]
output <- args[3]
## input = "/data/fuy2/gfp_pirna/results/2015-03-07/signals/test.GSV6.combined.depth"
## my.title <- "test.GSV6"
## output <- "test.GSV6.combined.depth.pdf"
df <- read.table(input)
df$V2 <- as.integer(df$V2)
df$V3 <- as.numeric(df$V3)
if (ncol(df) <= 3) {
    write("The input file only has 3 columns.", stderr())
    df$V4 <-0
} else {
    df$V4 <- as.numeric(df$V4)
}
plot_gene_signal(df, my.title, output)
