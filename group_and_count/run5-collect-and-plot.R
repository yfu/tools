df = data.frame()
for (i in (c("00", "02", "04", "07", "10", "12", "14", "17", "20", "42"))) {
    for (j in c("")) {
            counts <- c()
            feature.names <- c()
            fn <- list.files(pattern=paste(i, "dpp", j, "[.].*[.]count", sep=""))
            for (f in fn) {
                ## print(f)
                tmp <- read.table(f)
                count <- tmp[1,1]
                counts <- c(counts, count)
                f <- gsub("[0-9]+dpp.", "", f)
                f <- gsub(".count", "", f)
                f <- gsub("unox.", "", f)
                feature.names <- c(feature.names, f)
            }
            df <- rbind(df, counts)
            colnames(df) <- feature.names   
            rownames(df)[nrow(df)] <- paste(i, "dpp", j, sep="")
        }
}

idx <- which(colnames(df) == "norm_factor")
norm_factors <- df[, idx]

df.norm <- sweep(df, 1, norm_factors, "/")

df.norm <- df.norm * 1e6
write.table(df, file="read_counts.txt", sep="\t", quote=FALSE)
write.table(df.norm, file="read_counts_norm.txt", sep="\t", quote=FALSE)

library(reshape2)
d <- df.norm
d$dpp <- rownames(df.norm)
d <- melt(d, value.name="ppm", variable.name="category")
d$dpp <- factor(d$dpp)
d$category <- factor(d$category)

library(ggplot2)
library(gplots)

df.norm.ox <- df.norm[  grepl("^[0-9]+dpp$", rownames(df.norm)),  ]
df.norm.unox <- df.norm[  grepl("^[0-9]+dpp_unox$", rownames(df.norm)),  ]

pdf("heatmap.pdf")

heatmap.2(log10(as.matrix(df.norm.ox)), Rowv=F, Colv=F, dendrogram="none", main="Everybody and his brother: log, ox", trace="none")
heatmap.2(log10(as.matrix(df.norm.unox)), Rowv=F, Colv=F, dendrogram="none", main="Everybody and his brother: log, unox", trace="none")

idx <- grepl("rmsk", colnames(df.norm.ox))
heatmap.2( log10(as.matrix(df.norm.ox[, idx])), Rowv=F, Colv=F, dendrogram="none", main="rmsk: log, ox", trace="none", srtRow=30, srtCol=30, cexCol=0.8)

idx <- grepl("rmsk_", colnames(df.norm.ox))
heatmap.2( log10(as.matrix(df.norm.ox[, idx])), Rowv=F, Colv=F, dendrogram="none", main="rmsk: log, ox", trace="none", srtRow=30, srtCol=30, cexCol=0.8)
dev.off()


dd <- d[grepl("cluster", d$category), ]

ggplot(dd, aes(x=dpp, y=ppm, color=category, group=category)) + geom_point() + geom_line()

tmp <- data.frame(prepachytene.s=df.norm.ox$cluster.prepachytene.s, prepachytene_exon.s=df.norm.ox$cluster.prepachytene_exon.s)
rownames(tmp) <- rownames(df.norm.ox)
tmp$prepachytene_intron.s <- with(tmp, prepachytene.s - prepachytene_exon.s)
tmp$prepachytene.s <- NULL
tmp$dpp <- rownames(tmp)

tmp2<- melt(tmp, id.vars="dpp")

pdf("prepachyten_exon_intron.pdf", width=10)
ggplot(tmp2, aes(x=dpp, y=value, group=dpp, fill=variable)) + geom_bar(stat="identity", position="stack") + ggtitle("Prepachytene piRNAs mostly map to the exons of prepachytene genes")
ggplot(tmp2, aes(x=dpp, y=value, group=dpp, fill=variable)) + geom_bar(stat="identity", position="fill") + ggtitle("Prepachytene piRNAs mostly map to the exons of prepachytene genes")
dev.off()

tmp <- data.frame(protein_coding.s=df.norm.ox$protein_coding.s, protein_coding.exon.s=df.norm.ox$protein_coding.exon.s)
tmp$protein_coding.intron <- with(tmp, protein_coding.s - protein_coding.exon.s)
tmp$dpp <- rownames(df.norm.ox)
tmp$protein_coding.s <- NULL
tmp2 <- melt(tmp, id.vars="dpp")

pdf("protein_coding_exon_intron.pdf", width=10)
ggplot(tmp2, aes(x=dpp, y=value, group=dpp, fill=variable)) + geom_bar(stat="identity", position="stack") + ggtitle("Prepachytene piRNAs map to exons and introns of protein-coding genes")
ggplot(tmp2, aes(x=dpp, y=value, group=dpp, fill=variable)) + geom_bar(stat="identity", position="fill") + ggtitle("Prepachytene piRNAs map to exons and introns of protein-coding genes")
dev.off()

tmp <- data.frame(protein_coding = with(df.norm.ox, protein_coding.s+df.norm.ox$protein_coding.as),
                  protein_coding_upstream5k = with(df.norm.ox, protein_coding_upstream5k.as + protein_coding_upstream5k.s),
                  protein_coding_downstream5k = with(df.norm.ox, protein_coding_downstream5k.as + protein_coding_downstream5k.s),
                  rmsk = with(df.norm.ox, rmsk.s + rmsk.as),
                  cluster = with(df.norm.ox, cluster.s + cluster.as),
                  lincRNA = with(df.norm.ox, lincRNA.s + lincRNA.as),
                  pseudogene = with(df.norm.ox, pseudogene.s + pseudogene.as),
                  other = with(df.norm.ox, other)
           )
rownames(tmp) <- rownames(df.norm.ox)

pdf("pie.pdf")
for (i in 1:nrow(tmp)) {
    pie(as.numeric(tmp[i, ]), labels=colnames(tmp), col=rainbow(ncol(tmp)), main=paste("pie", rownames(tmp)[i]), edges=400)
}
dev.off()

tmp$dpp <- rownames(tmp)
tmp2 <- melt(tmp, id.vars="dpp", variable.name="category", value.name = "ppm")

ggplot(tmp2, aes(x=dpp, y=ppm))
