get.mean <- function(d) {
    my.mean <- mean(d[c(2,3,4)])
    my.mean
}

get.sd <- function(d) {
    my.sd <- sd(d[c(2,3,4)])
    my.sd
}

get.zscore <- function(df) {
    ## Get the zscore (for a two-coloumn table)
    bg.start = 0
    bg.end = 29
    idx <- df[, 1]==9
    tenth <- df[idx, 2]
    idx2 <- df[, 1] != 9 & df[, 1] >= bg.start & df[, 1] <= bg.end
    ## print( a[idx2, ] )
    others <- df[idx2, 2]
    mean.others <- mean(others)
    sd.others <- sd(others)
    zscore <- (tenth - mean.others) / sd.others
}

pdf("pps_with_errorbar.GSV6.uniq.pdf")
for(i in c("SRS.42A18AT1.br123.GSV6", "SRS.42A18BT1.br123.GSV6", "SRS.42A18CT1.br123.GSV6", "SRS.GS1AT1.br123.GSV6", "SRS.GS1BT1.br123.GSV6", "SRS.GS1CT1.br123.GSV6")) {
    ## lib = "SRS.42A18AT1.br123.GSV6.uniq"

    lib = i;
    fn <- paste(lib, ".uniq.pp_hist.norm.all", sep="")
    a <- read.table(fn)
    
    colnames(a) <- c("pos", paste("br", c(1,2,3), sep=""))

zscores.3 <- c(0, 0, 0)
for(j in c(2,3,4)) {
    df <- a[, c(1, j)]
    my.zscore <- get.zscore(df)
    zscores.3[j-1] <- my.zscore
}
mean.zscores <- mean(zscores.3)
sd.zscores <- sd(zscores.3)
    
    my.means <- apply(a, 1, get.mean)
    my.sds <- apply(a, 1, get.sd)
    
    my.ymin <- my.means - my.sds
    my.ymax <- my.means + my.sds
    
    df <- data.frame(pos=a$pos, y=my.means, ymin=my.ymin, ymax=my.ymax)

zscore <- get.zscore(df[, c(1,2)])

    ## idx <- df$pos==9
    ## tenth <- df[idx, 2]
    ## bg.start = 0
    ## bg.end = 29
    ## idx2 <- df$pos != 9 & df$pos >= bg.start & df$pos <= bg.end
    
    ## others <- df[idx2, "y"]
    ## mean.others <- mean(others)
    ## sd.others <- sd(others)
    ## zscore <- (tenth - mean.others) / sd.others

my.title <- paste(lib, "\nZ-score (from the averaged signals, 0-29 as background): ", zscore, "\n", "mean Z-score: ", mean.zscores, "\n", "sd Z-score: ", sd.zscores, "\n Z-score for individual library: ", paste(round(zscores.3, digits=3), collapse=", "), sep="")

    library(ggplot2)
    
    p <- ggplot(df, aes(x=pos, y=my.means, ymin=ymin, ymax=ymax)) + geom_bar(stat="identity") + geom_errorbar() + ggtitle(my.title) + theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())

print(p)
}
dev.off()
