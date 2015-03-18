#!/usr/bin/env Rscript

## This script plots each RNA species as a rectangle. The height of the rectangle represents
## the normalized abundance. The width represents the length of the RNA species.
## Author: Yu Fu

## Usage: plot_rna_rectangles.R input.bed2 200 output.pdf
## This second paramter specifies how long each window is for the detailed plots.

plot_rna_rectangles <- function (fn, each.length=100) {
    library(ggplot2)
    ## fn <- "Hi5.nodavirus.unox.bed2"
    a <- read.table(fn)
    a$ymin <- 0
    a$ymax <- 0
    idx <- a$V6 == "+"
    ## Set up how to randomly generate positions of rectangles
    r_max = max(a$V4/a$V5) / 2
    
    sum.plus <- sum(idx)
    a[idx, ]$ymin <- runif(sum.plus, 0, r_max)
    a[idx, ]$ymax <- a[idx, ]$ymin + a[idx, ]$V4 / a[idx, ]$V5
    sum.minus <- sum(!idx)
    a[!idx, ]$ymax <- -runif(sum.minus, 0, r_max)
    a[!idx, ]$ymin <- a[!idx, ]$ymax - a[!idx, ]$V4 / a[!idx, ]$V5
    
    a$V6 <- relevel(a$V6, "+")  # The default levels could be different. This will make sure that + is blue and - is red
    ## a$ymin[idx] <- -a$ymin[idx]
    ## a$V4[idx] <- -a$V4[idx]

    ret = list()
    p <- ggplot(data=a, aes(xmin=V2, xmax=V3, ymin=ymin, ymax=ymin+V4/V5, fill=V6)) + geom_rect() + scale_fill_manual(values=c("blue", "red")) + guides(fill=guide_legend(title="strand")) + ggtitle("Rectangles") + xlab("coordinate")
    ret[[length(ret) + 1]] <- p
    
    ## Divide the data into 500-nt chunks
    l <- each.length
    brk <- floor(max(a$V2) / l)
    aa <- cut(a$V2, breaks=brk)

    ## pdf("test.pdf")
    for(i in levels(aa)) {
        idx <- aa == i
        df <- a[idx, ]
        p <- ggplot(data=df, aes(xmin=V2, xmax=V3, ymin=ymin, ymax=ymin+V4/V5, fill=V6)) + geom_rect(alpha=0.3) + scale_fill_manual(values=c("blue", "red")) + guides(fill=guide_legend(title="strand")) + ggtitle("Rectangles") + xlab("coordinate")
        ret[[length(ret)+1]] <- p
        ## print(p)
    }
    ret
}

args <- commandArgs(TRUE)
fn <- args[1]
each.length <- as.integer(args[2])
output <- args[3]

## fn <- "Hi5.nodavirus.unox.bed2"
ps <- plot_rna_rectangles(fn, each.length=each.length)
pdf(output)
for (p in ps) {
    print(p)
}
dev.off()
