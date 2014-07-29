sample.name <- "Phil.SRA.wt.ox.ovary.inserts"
a <- read.table("pp.dataframe")

my.df <- data.frame(zscore=0)
a <- cbind(a, my.df)

for (i in 1:nrow(a)) {
    my.rows <- c(1:(i-1), (i+1):nrow(a))
    my.num <- a[my.rows, 2]
    my.mean <- mean(my.num)
    my.sd <- sd(my.num)
    z <- (a[i, 2] - my.mean) / my.sd
    if(is.na(z)) { z = 0 }
    a[i, 3] = z
}

library(ggplot2)

pdf.name <- paste(sample.name, ".pp_signal.pdf", sep="")
pdf(pdf.name)
ggplot(a, aes(x=V1, y=zscore)) + geom_line() + xlab("offset (10 means there is a 10-nt overlap)") + ylab("z-score") + ggtitle(sample.name)
dev.off()
