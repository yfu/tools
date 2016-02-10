
args <- commandArgs(TRUE)
## Input is just a two-column file: the first column being the read names and second columns as the lengths

## fn <- args[1]
## output <- args[2]

## fn <- "m151028_174034_42183_c100894422550000001823208304231690_s1_p0.lendis"
## a <- read.table(fn)

plot_fasta_size <- function(a) {
    library(ggplot2)
    my.median <- median(a$V2)
    my.mean <- mean(a$V2)
    n.reads <- nrow(a)
    my.title <- paste("# reads: ", n.reads, " median: ", my.median, " mean: ", my.mean)
    p <- ggplot(a, aes(x=V2)) + geom_histogram(binwidth=100) + ggtitle(my.title)
    list(p, p + xlim(c(0, 5000)))
}

for(i in list.files(".", "^m.+.lendis")) {
    a <- read.table(i)    
    my.plots <- plot_fasta_size(a)
    my.plots <- plot_fasta_size(a)
    pdf(paste(i, ".pdf", sep=""))
    for( p in my.plots ) {
        print(p)
    }
    dev.off()
}



## pdf("test.pdf")
## my.plots <- plot_fasta_size(a)
## for(p in my.plots) {
##     print(p)
## }
## dev.off()

