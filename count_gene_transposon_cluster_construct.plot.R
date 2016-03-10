## libs <- c(paste("RSRZ.CC1.br", 1:3, ".tr1.ovary", sep=""),
##           "RSRZ.GS1CT1.br{1,2,3}.tr1.ovary",
##           "RSRZ.CC2.br{1,2,3}.tr1.ovary", 
##           "RSRZ.GS1CT2.br{1,2,3}.tr1.ovary",
##           "RSRZ.CC3.br{1,2,3}.tr1.ovary",
##           "RSRZ.GS1CT3.br{1,2,3}.tr1.ovary",
##           "RSRZ.CC4.br{1,2,3}.tr1.ovary",
##           "RSRZ.GS1CT4.br{1,2,3}.tr1.ovary",
##           "RSRZ.AC1.br0.tr1.ovary",
##           "RSRZ.GS1AT1.br0.tr1.ovary",
##           "RSRZ.42A18AT1.br0.tr1.ovary",
##           "RSRZ.ZIPAT1.br0.tr1.ovary")
## libs1 <- paste("RSRZ.CC1.br", 1:3, ".tr1.ovary", sep="")
## libs2 <- paste("RSRZ.GS1CT1.br", 1:3, ".tr1.ovary", sep="")

args = commandArgs(trailingOnly=TRUE)

libs1 <- unlist(strsplit( args[1], ","))
libs2 <- unlist(strsplit( args[2], ","))
print(libs1)
output <- args[3]

## pdf("scatterplots.pdf")
pdf(output)
libs1.df <- data.frame()
libs2.df <- data.frame()
l <- length(libs1)
for (i in 1:l) {
    fn <- paste(libs1[i], ".count", sep="")
    a <- read.table(fn)
    a <- a[, c(1, 2, 3, 5)]
    fn <- paste(libs2[i], ".count", sep="")
    b <- read.table(fn)
    b <- b[, c(1, 2, 3, 5)]
    if (length(libs1.df) == 0) {
        libs1.df <- a
        libs2.df <- b
    } else {
        libs1.df <- merge(libs1.df, a, by=c("V1", "V2", "V3"))
        libs2.df <- merge(libs2.df, b, by=c("V1", "V2", "V3"))
    }
    colnames(libs1.df)[ncol(libs1.df)] <- paste(libs1[i], ".rep", i, sep="")
    colnames(libs2.df)[ncol(libs2.df)] <- paste(libs2[i], ".rep", i, sep="")
}

idx <- libs1.df$V3 == "Watson"
libs1.df$V3[idx] <- "S"
idx <- libs1.df$V3 == "Crick"
libs1.df$V3[idx] <- "AS"

idx <- libs2.df$V3 == "Watson"
libs2.df$V3[idx] <- "S"
idx <- libs2.df$V3 == "Crick"
libs2.df$V3[idx] <- "AS"

library(ggplot2)

both <- merge(libs1.df, libs2.df, by=c("V1", "V2", "V3"), all=T)

get.pval <- function(row) {
    t.test(as.numeric(row[4:6]), as.numeric(row[7:9]))$p.value  
}

idx <- is.na(both)
both[idx] <- 0

my.p.adj <- apply(both, 1, get.pval) * nrow(both)
my.p.adj[my.p.adj > 1] <- 1
both$p.adj <- my.p.adj
write.table(both, paste(output, ".dataframe.txt", sep=""), row.names=F, sep="\t", quote=F)

both$libs1.avg <- apply(both[4:6], 1, mean)
both$libs2.avg <- apply(both[7:9], 1, mean)

my.max <- max(both[, c("libs1.avg", "libs2.avg")])
my.title <- paste( "x: ", paste(libs1, collapse="\n"), "\nvs\n", "y: ", paste(libs2, collapse="\n"))

## p <- ggplot(sense, aes(x=libs1.avg, y=libs2.avg, color=V2)) + geom_point() + xlim(c(0, my.max)) + ylim(c(0, my.max)) + coord_fixed() + scale_x_log10() + scale_y_log10() + ggtitle(my.title)
for (g in c("gene", "transposon_family", "threeprimeutr", "cluster", "construct")) {
    for (strand in c("S", "AS")) {
        tmp <- both[both$V2==g & both$V3==strand, ]
        my.title.2 <- paste(my.title, "\n", strand)
        p <- ggplot(tmp, aes(x=libs1.avg, y=libs2.avg, color=V2)) + geom_point() + xlim(c(0, my.max)) + ylim(c(0, my.max)) + coord_fixed() + scale_x_log10() + scale_y_log10() + ggtitle(paste(my.title.2))
        print(p)
    }
}
dev.off()
