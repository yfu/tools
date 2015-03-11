#!/usr/bin/env Rscript

## Read as input the table from pps_simple.py and output a PDF file
## Keep in mind that in the input, i=9 means 10nt overlap

library(ggplot2)

args <- commandArgs(TRUE)
fn <- args[1]
output <- args[2]

## fn <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/Hi5.nodavirus.ox.pps"
## output <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/test.pdf"

a <- read.table(fn)
if( ncol(a) == 2) {
    ## This type of input
    ## 18	351
    ## 19	344
    ## 20	231
    ## 21	210
    pdf(output)
    idx <- a$V1==9
    tenth <- a[idx, 2]
    others <- a[!idx, 2]
    mean.others <- mean(others)
    sd.others <- sd(others)
    zscore <- (tenth - mean.others) / sd.others
    hi.pos <- a[max(a$V2) == a$V2, 1] + 1
    my.title <- paste("Ping-Pong z-score: ", zscore, "\nPosition with max signal: ", hi.pos, sep="")
    ## Add 1 to the 2nd column here to make the plots more visually pleasing...
    p <- ggplot(a, aes(x=V1, y=V2+1)) + geom_point() + geom_line() + xlab("# nucleotide overlap (10 means 10-nt overlap)") + ylab("raw signal") + ggtitle(my.title)
    print(p)
    dev.off()
} else {
    write("Wrong format: I want a two-column table from pps_simple.py", stderr())
}

