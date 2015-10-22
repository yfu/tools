#!/usr/bin/env Rscript

## Read as input the table from pps_simple.py and output a PDF file
## Keep in mind that in the input, i=9 means 10nt overlap

library(ggplot2)

args <- commandArgs(TRUE)
fn <- args[1]
output <- paste(fn, ".pdf", sep="")
# output <- args[2]
## fn <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/Hi5.nodavirus.unox.pps"
## output <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/test.pdf"

a <- read.table(fn)

pps_plot <- function(a) {
    if( ncol(a) == 2) {
        idx <- a$V1 >= 0
        pos <- a[idx, ]
        ## pos$V1 = pos$V1 + 1
        neg <- a[!idx, ]
        ## This type of input
        ## 18	351
        ## 19	344
        ## 20	231
        ## 21	210

        a$color = "0"
        a$color[ a$V1 %%  5 == 0 ] = "5"
        a$color[ a$V1 %% 10 == 0 ] = "10"
        a$color[ a$V1 %% 20 == 0 ] = "20"
        pdf(output)
        idx <- a$V1==9
        tenth <- a[idx, 2]
        others <- a[!idx, 2]
        mean.others <- mean(others)
        sd.others <- sd(others)
        zscore <- (tenth - mean.others) / sd.others
        ## hi.pos <- a[max(a$V2) == a$V2, 1]
        max.pos.positive <- head( pos[ order(-pos$V2), "V1"], n=10)
        max.pos.positive <- paste(max.pos.positive, collapse=", ")
        max.pos.negative <- head( neg[ order(-neg$V2), "V1"], n=10)
        max.pos.negative <- paste(max.pos.negative, collapse=", ")
        my.title <- paste("Ping-Pong z-score: ", zscore, "\nPositions with max signals: \n", max.pos.positive, "\n", max.pos.negative, sep="")
        ## Add 1 to the 2nd column here to make the plots more visually pleasing...
        ## p <- ggplot()
        p <- ggplot(a, aes(x=V1, y=V2)) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("raw signal") + ggtitle(my.title) + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1)
        p1 <- p + geom_bar(stat="identity")        
        p2 <- p + geom_bar(stat="identity", aes(fill=color, color=color)) + scale_fill_manual(values=c("black", "red", "green", "blue")) + scale_color_manual(values=c("black", "red", "green", "blue"))
        ## p <- p + geom_point(aes(x=V1, y=V2), data=neg) + geom_line(aes(x=V1, y=V2))
        p3 <- p2 + scale_y_log10()
        print(p1)
        print(p2)
        print(p3)
        dev.off()
    } else {
        write("Wrong format: a two-column table from pps_simple.py is expected", stderr())
    }
}

pps_plot(a)
