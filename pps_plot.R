#!/usr/bin/env Rscript

## Read as input the table from pps_simple.py and output a PDF file
## Keep in mind that in the input, i=9 means 10nt overlap
## Usage: pps_plot.R test.bed2.pp_hist 0 29 (use signals at [0, 29] except for the 9th position as the background)
## Usage: pps_plot.R test.bed2.pp_hist 1000000 (use everything except for the 9th position as the background and tell it that there is 1 million reads in total)
## Usage: pps_plot.R test.bed2.pp_hist (use everything except for the 9th position as the background

library(ggplot2)

args <- commandArgs(TRUE)
args <- c("/data/fuy2/gfp_pirna/results/2016-02-24-confirming-dm3-dm6-and-plotting-wiggle/updown2k/phasing/test/test.Watson.phasing_hist", "10000000", "0", "29")
fn <- args[1]
output <- paste(fn, ".pdf", sep="")
bg.start <- -Inf
bg.end <- Inf
## Total number of reads
total <- 1

if(length(args) == 2) {
    total <- as.numeric(args[2])
}
if(length(args)==4) {
    total <- as.numeric(args[2])
    bg.start <- as.numeric(args[3])
    bg.end <- as.numeric(args[4])
    sprintf("Using [%d, %d] except for the 9th position as the background!", bg.start, bg.end)
} else {
    ## bg.start <- -Inf
    ## bg.end <- Inf
    sprintf("Using all signals except for the 9th position as the background", bg.start, bg.end)
}

sprintf("Using %d reads to normalize", total)

add.zscore.col <- function(my.df) {
    ## Add a column of z-scores to the data frame
    ## my.df$zscore <- 0
    for (i in 1:nrow(my.df)) {
        bg <- my.df[, 1:2]
        idx <- bg[, 1] >= bg.start & bg[, 1] <= bg.end
        ## dont forget the remove itself
        idx[i] <- FALSE
        bg <- bg[idx, ]
        ## If the position in interest happens to be between bg.start and bg.end
        ## we need to remove that from the bg
        cur.score <- my.df[i, 2]
        if (cur.score >= bg.start & cur.score <= bg.end) {
            bg <- bg[-i, ]
        }
        my.mean <- mean(bg[, 2])
        my.sd <- sd(bg[, 2])
        my.df$zscore[i] <- (cur.score - my.mean) / my.sd
    }
    my.df
}


# output <- args[2]
## fn <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/Hi5.nodavirus.unox.pps"
## output <- "/data/fuy2/cl/results/2015-03-09/bowtie_mapping/test.pdf"
a <- read.table(fn)
## total <- 15093295

pps_plot <- function(a, total=1) {
    if( ncol(a) == 2) {
        a <- add.zscore.col(a)
        print("The third column is the zscore")
        print(a)
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
        idx2 <- a$V1 != 9 & a$V1 >= bg.start & a$V1 <= bg.end
        ## print( a[idx2, ] )
        others <- a[idx2, 2]
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
        q <- ggplot(a, aes(x=V1, y=zscore)) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("z-score given the background set by the user") + ggtitle(my.title) + geom_bar(stat="identity") + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1)
        print(p1)
        print(p2)
        print(p3)
        print(q)
        # print(q + scale_y_log10())
        ## Normalized value
        if(total!=1) {
            pn1 <- ggplot(a, aes(x=V1, y=V2 / (total * total)*1e3 )) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("pairs per thousand pairs") + ggtitle(my.title) + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1) + geom_bar(stat="identity")
            pn2 <- ggplot(a, aes(x=V1, y=V2 / (total * total)*1e6 )) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("pairs per million pairs") + ggtitle(my.title) + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1) + geom_bar(stat="identity")
            pn3 <- ggplot(a, aes(x=V1, y=V2 / (total * total)*1e9 )) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("pairs per billion pairs") + ggtitle(my.title) + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1) + geom_bar(stat="identity")
            pn4 <- ggplot(a, aes(x=V1, y=V2 / (total * total)*1e12 )) + xlab("# nucleotide overlap (x=9: 10-nt overlap)") + ylab("pairs per trillion pairs") + ggtitle(my.title) + geom_vline(xintercept = 0, color="gray") + theme(aspect.ratio=1) + geom_bar(stat="identity")
            print(pn1)
            print(pn2)
            print(pn3)
            print(pn4)
        }
        dev.off()
    } else {
        write("Wrong format: a two-column table from pps_simple.py is expected", stderr())
    }
}

pps_plot(a, total=total)
