## echo "SRS.42A18AT1.br1.tr1.ox.ovary.fq.gz SRS.42A18AT1.br2.tr1.ox.ovary.fq.gz SRS.42A18AT1.br3.tr1.ox.ovary.fq.gz SRS.GS1AT1.br1.tr1.ox.ovary.fq.gz SRS.GS1AT1.br2.tr1.ox.ovary.fq.gz SRS.GS1AT1.br3.tr1.ox.ovary.fq.gz SRS.42A18BT1.br1.tr1_2.ox.ovary.fq.gz SRS.42A18BT1.br2.tr1.ox.ovary.fq.gz SRS.42A18BT1.br3.tr1.ox.ovary.fq.gz SRS.GS1BT1.br1.tr1.ox.ovary.fq.gz SRS.GS1BT1.br2.tr1.ox.ovary.fq.gz SRS.GS1BT1.br3.tr1.ox.ovary.fq.gz SRS.42A18CT1.br1.tr1.ox.ovary.fq.gz SRS.42A18CT1.br2.tr1.ox.ovary.fq.gz SRS.42A18CT1.br3.tr1.ox.ovary.fq.gz SRS.GS1CT1.br1.tr1.ox.ovary.fq.gz SRS.GS1CT1.br2.tr1.ox.ovary.fq.gz SRS.GS1CT1.br3.tr1_2.ox.ovary.fq.gz SRS.GS1CT1.br3.tr2.ox.ovary.fq.gz" | sed 's/ /\n/g' | sed 's/.fq.gz//' | sed 's/^/"/' | sed 's/$/"/' | tr '\n' ','
libs <- c("SRS.42A18AT1.br1.tr1.ox.ovary","SRS.42A18AT1.br2.tr1.ox.ovary","SRS.42A18AT1.br3.tr1.ox.ovary","SRS.GS1AT1.br1.tr1.ox.ovary","SRS.GS1AT1.br2.tr1.ox.ovary","SRS.GS1AT1.br3.tr1.ox.ovary","SRS.42A18BT1.br1.tr1_2.ox.ovary","SRS.42A18BT1.br2.tr1.ox.ovary","SRS.42A18BT1.br3.tr1.ox.ovary","SRS.GS1BT1.br1.tr1.ox.ovary","SRS.GS1BT1.br2.tr1.ox.ovary","SRS.GS1BT1.br3.tr1.ox.ovary","SRS.42A18CT1.br1.tr1.ox.ovary","SRS.42A18CT1.br2.tr1.ox.ovary","SRS.42A18CT1.br3.tr1.ox.ovary","SRS.GS1CT1.br1.tr1.ox.ovary","SRS.GS1CT1.br2.tr1.ox.ovary","SRS.GS1CT1.br3.tr1_2.ox.ovary")
## CHANGE here to switch between using all mappers or unique mappers
au = "uniq"
# libs <- libs[seq(1, 6)]
df <- data.frame()
for (i in libs) {
    one.sample <- data.frame()
    for (j in c("transposon", "cluster", "construct")) {
        fn <- paste(i, ".", au, ".", j, ".abundance.normalized_by_allxmirna", sep="")
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

g1.x <- c("SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary")
g1.y <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary")
g2.x <- c("SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary")
g2.y <- c("SRS.42A18BT1.br1.tr1_2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary")
g3.x <- c("SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary")
g3.y <- c("SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary")
g1.x.lab <- "SRS.GS1AT1"
g1.y.lab <- "SRS.42A18AT1"
g2.x.lab <- "SRS.GS1BT1"
g2.y.lab <- "SRS.42A18BT1"
g3.x.lab <- "SRS.GS1CT1"
g3.y.lab <- "SRS.42A18CT1"

df.mean <- df[, c(1, 2)]
for (g in list(g1.x, g1.y, g2.x, g2.y, g3.x, g3.y)) {
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

for(xy in list(list(g1.x, g1.y), list(g2.x, g2.y), list(g3.x, g3.y))) {
    x <- xy[[1]]
    n.x <- gsub("(SRS.[^.]+).+", "\\1", x[[1]])[1]
    y <- xy[[2]]
    n.y <- gsub("(SRS.[^.]+).+", "\\1", y[[1]])[1]
    s.x <- df[, paste(x, ".S", sep="")]
    s.y <- df[, paste(y, ".S", sep="")]
    as.x <- df[, paste(x, ".AS", sep="")]
    as.y <- df[, paste(y, ".AS", sep="")]
    a.xy <- cbind(s.x, s.y)
    as.xy <- cbind(as.x, as.y)
    ## Total = sense + antisense. The colnames are not accurate here...
    sas.xy <- a.xy + as.xy
    ## Add a very small random number to every element so that t.test does not complain...
    p <- apply(sas.xy, 1, f)
    p.adj <- p.adjust(p, method="fdr")
    p.adj <- data.frame(p.adj)
    colnames(p.adj) <- paste(n.x, ".vs.", n.y, ".fdr", sep="")
    p <- data.frame(p)
    colnames(p) <- paste(n.x, ".vs.", n.y, ".p", sep="")
    df.p.fdr <- cbind(df.p.fdr, p)        
    df.p.fdr <- cbind(df.p.fdr, p.adj)
}

write.table(df.p.fdr, paste(au, ".scatterplots.transposon.cluster.construct.p.fdr.table.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

lim <- max( df.mean[, seq(3, ncol(df.mean))])
lim <- 10^ceiling(log10(lim))
idx <- grepl("tg", df.mean$group)
df.mean.tg <- df.mean[idx, ]
idx <- grepl("construct", df.mean$group)
df.mean.construct <- df.mean[idx, ]
idx <- grepl("cluster", df.mean$group)
df.mean.cluster <- df.mean[idx, ]
      

pdf(paste(au, ".scatterplots.transposon.cluster.construct.pdf", sep=""))
for (i in list( list(g1.x.lab, g1.y.lab), list(g2.x.lab, g2.y.lab), list(g3.x.lab, g3.y.lab))) {
    x <- i[[1]]
    y <- i[[2]]
    for (j in c(".S", ".AS", ".SAS")) {
        cn1 <- paste(x, j, sep="")
        cn2 <- paste(y, j, sep="")
        p1 <- one_plot(df.mean.tg, cn1, cn2, cn1, cn2, lim=lim)
        p2 <- one_plot(df.mean.construct, cn1, cn2, cn1, cn2, lim=lim)
        p3 <- one_plot(df.mean.cluster, cn1, cn2, cn1, cn2, lim=lim)
        print(p1)
        print(p2)
        print(p3)
    }
}
dev.off()

write.table(df, file=paste(au, ".scatterplots.transposon.cluster.construct.table.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
