df <- data.frame()
norm_table <- data.frame()

libs <- c("SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary", "SRS.GS1CT1.br3.tr1.ox.ovary", "SRS.GS1CT1.br3.tr2.ox.ovary")
for (D in libs) {
        table <- paste(D, "/", D, ".table", sep="")
        norm_table_tmp <- read.table(table, sep="\t")
        colnames(norm_table_tmp) <- c("norm_type", D)
        if(length(norm_table)==0) {
            norm_table <- norm_table_tmp
        } else {
            norm_table <- merge(norm_table, norm_table_tmp, by="norm_type")
        }
        gcr_d <- paste(D, "/", "gene+cluster+repBase", sep="")
        cluster_count <- read.table( paste(gcr_d, "/", D, ".gene+cluster+repBase.bed2.cluster.count", sep="") )
        ## cluster_count1$type="cluster"
        gene_count <- read.table(paste(gcr_d, "/", D, ".gene+cluster+repBase.bed2.gene.count", sep=""))
        ## gene_count$type="gene"
        repBase_count <- read.table(paste(gcr_d, "/", D, ".gene+cluster+repBase.bed2.repBase.count", sep=""))
        ## repBase_count$type="repBase"
        construct_count <- read.table(paste(gcr_d, "/", D, ".gene+cluster+repBase.bed2.construct.count", sep=""))
        ## construct_count$type="construct"
        if( length(df) ==0 ) {
            df <- rbind(cluster_count, gene_count, repBase_count, construct_count)
            colnames(df) <- c("feature", D)
        } else {
            tmp <- rbind(cluster_count, gene_count, repBase_count, construct_count)
            colnames(tmp) <- c("feature", D)
            df <- merge(df, tmp, by="feature")
        }
    }

df.backup <- df

idx <- grepl("^NM_|^NR_", df$feature)
df$type="x"
df$type[idx] <- "gene"
idx <- grep("^FBgn", df$feature)
df$type[idx] <- "repBase"
idx <- grep("^cluster|^42AB|^flam", df$feature)
df$type[idx] <- "cluster"
idx <- grep("GSV6|nos-gal4-vp16|UASp-EGFP", df$feature)
df$type[idx] <- "construct"

pdf("scatter_plots.pdf")
## for (x in libs) {
for (n in norm_table$norm_type) {
    line <- which(norm_table$norm_type==n)
    ## for (x in "SRS.GS1CT1.br1.tr1.ox.ovary") {
    for(x in libs) { 
        nf1 <- norm_table[line, x]
        dx <- df[, c(x, "type")]
        dx[, 1] <- dx[, 1] / nf1 * 1e6
        for (y in libs) {
        ## for (y in "SRS.GS1CT1.br2.tr1.ox.ovary") {
            nf2 <- norm_table[line, y]
            dy <- df[, c(y, "type")]
            dy[, 1] <- dy[, 1] / nf2 * 1e6
            my <- cbind(dx, dy)
            p <- ggplot(data=my, aes(x=SRS.GS1CT1.br1.tr1.ox.ovary, y=SRS.GS1CT1.br2.tr1.ox.ovary, color=type)) + geom_point() + geom_abline(intercept = 0, slope = 1) + scale_x_log10() + scale_y_log10() + ggtitle(paste("Normalized to 1 million of ", n, sep=""))
            print(p)
        }
    }   
}
dev.off()
