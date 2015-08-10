## libs <- c("SRS.42A18AT1.br1.tr1.ox.ovary","SRS.42A18AT1.br2.tr1.ox.ovary","SRS.42A18AT1.br3.tr1.ox.ovary","SRS.GS1AT1.br1.tr1.ox.ovary","SRS.GS1AT1.br2.tr1.ox.ovary","SRS.GS1AT1.br3.tr1.ox.ovary","SRS.42A18BT1.br1.tr1_2.ox.ovary","SRS.42A18BT1.br2.tr1.ox.ovary","SRS.42A18BT1.br3.tr1.ox.ovary","SRS.GS1BT1.br1.tr1.ox.ovary","SRS.GS1BT1.br2.tr1.ox.ovary","SRS.GS1BT1.br3.tr1.ox.ovary","SRS.42A18CT1.br1.tr1.ox.ovary","SRS.42A18CT1.br2.tr1.ox.ovary","SRS.42A18CT1.br3.tr1.ox.ovary","SRS.GS1CT1.br1.tr1.ox.ovary","SRS.GS1CT1.br2.tr1.ox.ovary","SRS.GS1CT1.br3.tr1_2.ox.ovary")

libs <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary", "SRS.42A18BT1.br1.tr1_2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary", "SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary", "SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary", "SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary", "SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary", "SRS.GS1DT1.br1.tr1_2.ox.ovary", "SRS.GS1DT1.br2.tr1.ox.ovary", "SRS.GS1DT1.br3.tr1.ox.ovary", "SRS.GS1ET1.br1.tr1.ox.ovary", "SRS.GS1ET1.br2.tr1.ox.ovary", "SRS.GS1ET1.br3.tr1.ox.ovary", "SRS.GS1FT1.br1.tr1.ox.ovary", "SRS.GS1FT1.br2.tr1.ox.ovary", "SRS.GS1FT1.br3.tr1_2.ox.ovary", "SRS.GS1GT1.br1.tr1.ox.ovary", "SRS.GS1GT1.br2.tr1.ox.ovary", "SRS.GS1GT1.br3.tr1.ox.ovary")

d <- "/home/fuy2/data/gfp_pirna/results/2015-07-15-check-all-replicates"
## CHANGE here to switch between using all mappers or unique mappers
au = "uniq"
# libs <- libs[seq(1, 6)]
df <- data.frame()
for (i in libs) {
    one.sample <- data.frame()
    for (j in c("transposon", "cluster", "construct")) {
        fn <- paste(d, "/", i, ".", au, ".", j, ".abundance.normalized_by_allxmirna", sep="")
        tmp <- read.table(fn)
        if (j == "transposon") {
            colnames (tmp) = c("name", "group", paste(i, ".S", sep=""), paste(i, ".AS", sep=""))
            tmp$group <- paste("tg.", tmp$group, sep="")
            one.sample <- tmp
        } else {
            colnames(tmp) <- c("name", paste(i, ".S", sep=""), paste(i, ".AS", sep=""))
            tmp$group <- j
            tmp <- tmp[, c(1,4,2,3)]
            one.sample <- rbind(one.sample, tmp)
        }
    }
    if (length(df) == 0) {
        df <- one.sample
    } else {
        # all=TRUE becasue some libraries have nasty 0 count
        df <- merge(df, one.sample, by = c("name", "group"), all=TRUE)
    }    
}
df$name <- as.character(df$name)
df <- df[ order(df$name), ]

### Some constructs may not be visible because they have no unique mappers
## For those constructs with NA count, set them to 0
df[is.na(df)] = 0
## idx <- df$group=="construct"
## df.tmp <- df[idx, seq(3, ncol(df))]
## df.tmp[df.tmp < 1] = 1
## df[idx, seq(3, ncol(df))] = df.tmp

library(ggplot2)
library(ggthemes)
library(reshape)
library(gridExtra)
library(scales)
library(grid)

## lim1 = floor   (min( df[, seq(3, ncol(df))])) - 0.5
## lim2 = ceiling (max( df[, seq(3, ncol(df))])) + 0.5

lim <- max( df[, seq(3, ncol(df))])
lim <- 10^ceiling(log10(lim))
idx <- grepl("tg", df$group)
df.tg <- df[idx, ]
idx <- grepl("construct", df$group)
df.construct <- df[idx, ]
idx <- grepl("cluster", df$group)
df.cluster <- df[idx, ]

one_plot <- function(df, cn1, cn2, name1, name2, lim) {
    ## cn1 and cn2: colname1 and colname2
    sample <- df
    gg=ggplot( sample, aes_string(x = cn1, y = cn2, color="group") ) +
        theme (
            plot.margin=unit(c(1,1,0,0),"lines"),
            axis.text=element_text (size=4),
            axis.title=element_text(size=6),
            legend.margin=unit(0,"lines"),
            panel.margin=unit(0, "lines"),
            axis.ticks.margin=unit(0,"lines"),
            legend.key.size=unit(0.5,"lines")
            ) +
                scale_x_log10 ( limits = c(1,lim), breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x) ) ) +
                                   scale_y_log10 ( limits = c(1,lim), breaks = trans_breaks("log10", function(x) 10^x),
                                                  labels = trans_format("log10", math_format(10^.x) ) ) +
                                                          annotation_logticks () +
                scale_size_manual( values=c(0,8) ) +
                    geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') +
                        geom_abline (intercept=log10(2), slope=1, colour="darkgrey", linetype='twodash') +
                            geom_abline (intercept=-log10(2), slope=1, colour="darkgrey", linetype='twodash') +

                    theme_few () +
                        scale_fill_continuous(guide = "legend") +
                        geom_point(size=5, alpha=0.75, na.rm=T) +
                                    scale_colour_manual(values=c("lightblue","black","darkgreen","red","blue","orange")) +
                                        guides(colour = guide_legend (title=expression (paste ("type")), title.position = "top")) +
                                            xlab ( substitute ( paste(italic(name1)), list(name1=name1, name2=name2))) +
                                                ylab ( substitute ( paste(italic(name2)), list(name1=name1, name2=name2))) +
                                                    coord_fixed(ratio=1, xlim = c(1, lim), ylim = c(1, lim))
    gg
}


## g1.x <- c("SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary")
## g1.y <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary")
## g2.x <- c("SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary")
## g2.y <- c("SRS.42A18BT1.br1.tr1_2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary")
## g3.x <- c("SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary")
## g3.y <- c("SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary")
## g1.x.lab <- "SRS.GS1AT1"
## g1.y.lab <- "SRS.42A18AT1"
## g2.x.lab <- "SRS.GS1BT1"
## g2.y.lab <- "SRS.42A18BT1"
## g3.x.lab <- "SRS.GS1CT1"
## g3.y.lab <- "SRS.42A18CT1"

g1 <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary")
g2 <- c("SRS.42A18BT1.br1.tr1_2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary")
g3 <- c("SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary")
g4 <- c("SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary")
g5 <- c("SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary")
g6 <- c("SRS.GS1CT1.br1.tr1.ox.ovary",  "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary")
g7 <- c("SRS.GS1DT1.br1.tr1_2.ox.ovary", "SRS.GS1DT1.br2.tr1.ox.ovary", "SRS.GS1DT1.br3.tr1.ox.ovary")
g8 <- c("SRS.GS1ET1.br1.tr1.ox.ovary", "SRS.GS1ET1.br2.tr1.ox.ovary", "SRS.GS1ET1.br3.tr1.ox.ovary")
g9 <- c("SRS.GS1FT1.br1.tr1.ox.ovary", "SRS.GS1FT1.br2.tr1.ox.ovary", "SRS.GS1FT1.br3.tr1_2.ox.ovary")
g10 <- c("SRS.GS1GT1.br1.tr1.ox.ovary", "SRS.GS1GT1.br2.tr1.ox.ovary", "SRS.GS1GT1.br3.tr1.ox.ovary")

n1 <- "SRS.42A18AT1"
n2 <- "SRS.42A18BT1"
n3 <- "SRS.42A18CT1"
n4 <- "SRS.GS1AT1"
n5 <- "SRS.GS1BT1"
n6 <- "SRS.GS1CT1"
n7 <- "SRS.GS1DT1"
n8 <- "SRS.GS1ET1"
n9 <- "SRS.GS1FT1"
n10<- "SRS.GS1GT1"

df.mean <- df[, c(1, 2)]
## for (g in list(g1.x, g1.y, g2.x, g2.y, g3.x, g3.y)) {
for (g in list(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10)) {
    s <- apply(df[, paste(g, ".S", sep="")], 1, mean)
    df.mean <- cbind(df.mean, s)
    n <- gsub("(SRS.[^.]+).+", "\\1", g)[1]
    colnames(df.mean)[ ncol(df.mean) ] <- paste(n, ".S", sep="")
    as <- apply(df[, paste(g, ".AS", sep="")], 1, mean)
    df.mean <- cbind(df.mean, as)
    n <- gsub("(SRS.[^.]+).+", "\\1", g)[1]
    colnames(df.mean)[ ncol(df.mean) ] <- paste(n, ".AS", sep="")
    df.mean <- cbind(df.mean, s+as)
    colnames(df.mean)[ ncol(df.mean) ] <- paste(n, ".SAS", sep="")    
}

df.p.fdr <- data.frame(name=df$name)
f <- function(x) {
    ## t.test(x[c(1,2,3)], x[c(4,5,6)]) [["p.value"]]
    if (all(x==x[1])) {
        1
    } else {
        t.test(x[c(1,2,3)], x[c(4,5,6)]) [["p.value"]] 
    }
}

f2 <- function(x) {
    ## given 6 numbers, calc the std of the first 3 and last 3
    c(sd( x[c(1,2,3)] ), sd( x[c(4,5,6)] ) )
}

f3 <- function(x) {
    ## given 6 numbers, calc the mean of the first 3 and last 3
    c(mean( x[c(1,2,3)] ), mean( x[c(4,5,6)] ) )
}

## Given a pair of samples, calculate the p-vals and FDRs
# df.output <- cbind(name=as.character(df$name), group=df$group)
## for(xy in list(list(g1.x, g1.y), list(g2.x, g2.y), list(g3.x, g3.y))) {
output.p.fdr <- function(d, x.group, y.group) {
    x <- x.group
    n.x <- gsub("(SRS.[^.]+).+", "\\1", x[[1]])[1]
    y <- y.group
    n.y <- gsub("(SRS.[^.]+).+", "\\1", y[[1]])[1]
    s.x <- d[, paste(x, ".S", sep="")]
    s.y <- d[, paste(y, ".S", sep="")]
    as.x <- d[, paste(x, ".AS", sep="")]
    as.y <- d[, paste(y, ".AS", sep="")]
    a.xy <- cbind(s.x, s.y)
    as.xy <- cbind(as.x, as.y)
    ## Total = sense + antisense. The colnames are not accurate here...
    sas.xy <- a.xy + as.xy
    sas.std <- t(apply(sas.xy, 1, f2))
    colnames(sas.std) <- paste("std.", c(n.x, n.y), sep="")
    sas.mean <- t(apply(sas.xy, 1, f3))
    colnames(sas.mean) <- paste("mean.", c(n.x, n.y), sep="")
    ## Add a very small random number to every element so that t.test does not complain...
    p <- apply(sas.xy, 1, f)
    p.adj <- p.adjust(p, method="fdr")
    p.adj <- data.frame(p.adj)
    colnames(p.adj) <- paste(n.x, ".vs.", n.y, ".fdr", sep="")
    p <- data.frame(p)
    colnames(p) <- paste(n.x, ".vs.", n.y, ".p", sep="")
    ## df.p.fdr <- cbind(df.p.fdr, p)        
    ## df.p.fdr <- cbind(df.p.fdr, p.adj)
    cbind(s.x, as.x, s.y, as.y, sas.mean, sas.std, p, p.adj)
}

a <- cbind(combn(c("g1", "g2", "g3"), m=2), combn(c("g4", "g5", "g6", "g7", "g8", "g9", "g10"), m=2))
for (i in seq(1, ncol(a))) {
    x.group <- eval(parse(text=a[1, i]))
    y.group <- eval(parse(text=a[2, i]))
    n.x <- gsub("(SRS.[^.]+).+", "\\1", x.group[[1]])[1]
    n.y <- gsub("(SRS.[^.]+).+", "\\1", y.group[[1]])[1]    
    tmp <- output.p.fdr(df, x.group, y.group)
    message(paste(n.x, ".vs.", n.y, ".table.txt", sep=""))
    write.table(tmp, file=paste(n.x, ".vs.", n.y, ".txt", sep=""), row.names=FALSE, quote=FALSE)
}

lim <- max( df.mean[, seq(3, ncol(df.mean))])
lim <- 10^ceiling(log10(lim))
idx <- grepl("tg", df.mean$group)
df.mean.tg <- df.mean[idx, ]
idx <- grepl("construct", df.mean$group)
df.mean.construct <- df.mean[idx, ]
idx <- grepl("cluster", df.mean$group)
df.mean.cluster <- df.mean[idx, ]
      
# pdf(paste(au, ".scatterplots.transposon.cluster.construct.pdf", sep=""))
## for (i in list( list(g1.x.lab, g1.y.lab), list(g2.x.lab, g2.y.lab), list(g3.x.lab, g3.y.lab))) {
output.plot <- function(x, y) {
    ## x <- i[[1]]
    ## y <- i[[2]]
    ret <- list()
    for (j in c(".S", ".AS", ".SAS")) {
        cn1 <- paste(x, j, sep="")
        cn2 <- paste(y, j, sep="")
        ## Just use the global variables: df.mean.*
        p1 <- one_plot(df.mean.tg, cn1, cn2, cn1, cn2, lim=lim)
        p2 <- one_plot(df.mean.construct, cn1, cn2, cn1, cn2, lim=lim)
        p3 <- one_plot(df.mean.cluster, cn1, cn2, cn1, cn2, lim=lim)
        ## list(p1, p2, p3)
        ret <- c(ret, list(p1, p2, p3))
    }
    ret
}

pdf(paste(au, ".scatterplots.transposon.cluster.construct.pdf", sep=""))
b <- cbind(combn(c("n1", "n2", "n3"), m=2), combn(c("n4", "n5", "n6", "n7", "n8", "n9", "n10"), m=2))
## b <- b[, c(1,2)]
for (i in seq(1, ncol(b))) {
    x <- eval(parse(text=b[1, i]))
    y <- eval(parse(text=b[2, i]))
    ps <- output.plot(x, y)
    for (p in ps) {
        print(p)
    }
}
dev.off()
