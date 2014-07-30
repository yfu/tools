library(ggplot2)

sample.name <- "Hi5.ox"
a <- read.table("pp.dataframe")

my.df <- data.frame(zscore=0)
a <- cbind(a, my.df)

pdf.name <- paste(sample.name, ".pp_raw_signal.pdf", sep="")
pdf(pdf.name)
p <- ggplot(a, aes(x=V1, y=V2)) + geom_line() + geom_point() + xlab("offset (10 means there is a 10-nt overlap)") + ylab("raw signal") + ggtitle(sample.name)
print(p)
dev.off()

for (i in 1:nrow(a)) {
    my.rows <- 1:nrow(a)
    my.rows <- my.rows[-i]
    my.num <- a[my.rows, 2]
    my.mean <- mean(my.num)
    my.sd <- sd(my.num)
    z <- (a[i, 2] - my.mean) / my.sd
    # message(my.sd)
    if(is.na(z)) { z = 0 }
    a[i, 3] = z
}

library(ggplot2)

pdf.name <- paste(sample.name, ".pp_zscore.pdf", sep="")
pdf(pdf.name)
ggplot(a, aes(x=V1, y=zscore)) + geom_line() + geom_point() + xlab("offset (10 means there is a 10-nt overlap)") + ylab("z-score") + ggtitle(sample.name)
dev.off()
