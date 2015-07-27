# Do not use SRS.42A18AT1.br4.tr1.ox.ovary according to Chengjian
groups <- list(c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary"),
               c("SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary"),
               c("SRS.42A18BT1.br1.tr1.ox.ovary", "SRS.42A18BT1.br1.tr2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary"),
               c("SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary"),
               c("SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary"),
               c("SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1.ox.ovary", "SRS.GS1CT1.br3.tr2.ox.ovary"))
err <- list()
gn <- 0

# groups <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary")
# Read in all data
df <- data.frame()
norm_table <- data.frame()
for (l in unlist(groups)) {
    table <- paste(l, "/", l, ".table", sep="")
    norm_table_tmp <- read.table(table, sep="\t")
    colnames(norm_table_tmp) <- c("norm_type", l)
    if(length(norm_table)==0) {
        norm_table <- norm_table_tmp
    } else {
        norm_table <- merge(norm_table, norm_table_tmp, by="norm_type")
        }
    gcr_d <- paste(l, "/", "gene+cluster+repBase", sep="")
    cluster_count <- read.table( paste(gcr_d, "/", l, ".gene+cluster+repBase.bed2.cluster.count", sep="") )
    ## cluster_count1$type="cluster"
    gene_count <- read.table(paste(gcr_d, "/", l, ".gene+cluster+repBase.bed2.gene.count", sep=""))
    ## gene_count$type="gene"
    repBase_count <- read.table(paste(gcr_d, "/", l, ".gene+cluster+repBase.bed2.repBase.count", sep=""))
    ## repBase_count$type="repBase"
    construct_count <- read.table(paste(gcr_d, "/", l, ".gene+cluster+repBase.bed2.construct.count", sep=""))
    ## construct_count$type="construct"
    if( length(df) ==0 ) {
        df <- rbind(cluster_count, gene_count, repBase_count, construct_count)
        colnames(df) <- c("feature", l)
    } else {
        tmp <- rbind(cluster_count, gene_count, repBase_count, construct_count)
        colnames(tmp) <- c("feature", l)
        df <- merge(df, tmp, by="feature")
    }
}

idx <- grepl("^NM_|^NR_", df$feature)
df$type="x"
df$type[idx] <- "gene"
idx <- grep("^FBgn", df$feature)
df$type[idx] <- "repBase"
idx <- grep("^cluster|^42AB|^flam", df$feature)
df$type[idx] <- "cluster"
idx <- grep("GSV6|nos-gal4-vp16|UASp-EGFP", df$feature)
df$type[idx] <- "construct"
df.backup <- df

# Merge the technical replicates
df$`SRS.42A18BT1.br1.tr1_2.ox.ovary` <- df$`SRS.42A18BT1.br1.tr1.ox.ovary` + df $`SRS.42A18BT1.br1.tr2.ox.ovary`
df$`SRS.GS1CT1.br3.tr1_2.ox.ovary` <- df$`SRS.GS1CT1.br3.tr1.ox.ovary` + df$`SRS.GS1CT1.br3.tr2.ox.ovary`

df$`SRS.42A18BT1.br1.tr1.ox.ovary` <- NULL
df$`SRS.42A18BT1.br1.tr2.ox.ovary` <- NULL
df$`SRS.GS1CT1.br3.tr1.ox.ovary` <- NULL
df$`SRS.GS1CT1.br3.tr2.ox.ovary` <- NULL

norm_table.backup <- norm_table

norm_table$`SRS.42A18BT1.br1.tr1_2.ox.ovary` <- norm_table$`SRS.42A18BT1.br1.tr1.ox.ovary` + norm_table$`SRS.42A18BT1.br1.tr2.ox.ovary`
norm_table$`SRS.GS1CT1.br3.tr1_2.ox.ovary` <- norm_table$`SRS.GS1CT1.br3.tr1.ox.ovary` + norm_table$`SRS.GS1CT1.br3.tr2.ox.ovary`
norm_table$`SRS.42A18BT1.br1.tr1.ox.ovary` <- NULL
norm_table$`SRS.42A18BT1.br1.tr2.ox.ovary` <- NULL
norm_table$`SRS.GS1CT1.br3.tr1.ox.ovary` <- NULL
norm_table$`SRS.GS1CT1.br3.tr2.ox.ovary` <- NULL


## x1 <- c("SRS.GS1AT1.br1.tr1.ox.ovary", "SRS.GS1AT1.br2.tr1.ox.ovary", "SRS.GS1AT1.br3.tr1.ox.ovary")
## y1 <- c("SRS.42A18AT1.br1.tr1.ox.ovary", "SRS.42A18AT1.br2.tr1.ox.ovary", "SRS.42A18AT1.br3.tr1.ox.ovary")
## x2 <- c("SRS.GS1BT1.br1.tr1.ox.ovary", "SRS.GS1BT1.br2.tr1.ox.ovary", "SRS.GS1BT1.br3.tr1.ox.ovary")
## y2 <- c("SRS.42A18BT1.br1.tr1.ox.ovary", "SRS.42A18BT1.br1.tr2.ox.ovary", "SRS.42A18BT1.br2.tr1.ox.ovary", "SRS.42A18BT1.br3.tr1.ox.ovary")
## y3 <- c("SRS.42A18CT1.br1.tr1.ox.ovary", "SRS.42A18CT1.br2.tr1.ox.ovary", "SRS.42A18CT1.br3.tr1.ox.ovary")
## x3 <- c("SRS.GS1CT1.br1.tr1.ox.ovary", "SRS.GS1CT1.br2.tr1.ox.ovary", "SRS.GS1CT1.br3.tr1_2.ox.ovary")

x1 <- "SRS.GS1AT1"
y1 <- "SRS.42A18AT1"
x2 <- "SRS.GS1BT1"
y2 <- "SRS.42A18BT1"
x3 <- "SRS.GS1CT1"
y3 <- "SRS.42A18CT1"

pdf("genotype_xy.pdf")
ll <- list(list(x1, y1), list(x2, y2) ,list(x3, y3))
for ( norm_type in c("Total", "flam unique mappers")) {
    for (l in ll) {
        x.pat <- l[[1]]
        y.pat <- l[[2]]
        ## Get the x and y values (before normalization)
        ## x.pat <- "SRS.42A18AT1"
        idx <- grepl(paste(x.pat, ".*", sep=""), colnames(df))
        message("x axis before norm:")
        message(paste(colnames(df)[idx], sep="", collapse="\n"))
        x.val <- apply(df[, idx], 1, function(x) {sum(x)/3} )
        ## y.pat <- "SRS.42A18BT1"
        idx <- grepl(paste(y.pat, ".*", sep=""), colnames(df))
        message("y axis before norm:")
        message(paste(colnames(df)[idx], sep="", collapse="\n"))
        y.val <- apply(df[, idx], 1, function(x) {sum(x)/3} )
        xy <- data.frame(x.val, y.val)
        colnames(xy) <- c(x.pat, y.pat)
        xy$feature <- df$feature
        ## Get norm factor
        rownames(norm_table) <- norm_table$norm_type
        ## set it in the outer loop
        ## norm_type <- 'Total'
        
        idx <- grepl(paste(x.pat, ".*", sep=""), colnames(norm_table))
        x.nf <- 1e6 / apply(norm_table[norm_type, idx], 1, sum)
        idx <- grepl(paste(y.pat, ".*", sep=""), colnames(norm_table))
        y.nf <- 1e6 / apply(norm_table[norm_type, idx], 1, sum)
        ## Data for plotting
        ## pseudocount
        pc = 1
        dp <- data.frame(x.val * x.nf+1, y.val * y.nf+1)
        colnames(dp) <- c(x.pat, y.pat)
        library(ggplot2)
        dp$type <- df$type
        dp$feature <- df$feature
        my.max <- max(x.val * x.nf, y.val * y.nf)
        
        ## Notice pseudocount is 1
        p <- ggplot(dp, aes_string(x=x.pat, y=y.pat, color="type", group="type")) + geom_point() + scale_x_log10(limits=c(1, my.max)) + scale_y_log10(limits=c(1, my.max)) + ggtitle(paste("norm: ", norm_type, ", norm: ", norm_type)) + geom_abline(slope=1, intercept=0)
        print(p)
        for (t in c("cluster", "gene", "construct", "repBase")) {
            dp.tmp <- dp[dp$type == t, ]
            p <- ggplot(dp.tmp, aes_string(x=x.pat, y=y.pat)) + geom_point() + scale_x_log10(limits=c(1, my.max)) + scale_y_log10(limits=c(1, my.max)) + ggtitle(paste("norm: ", norm_type, ",feature: ", t, " ,norm: ", norm_type)) + geom_abline(slope=1, intercept=0)    
            print(p)
        }
    }
}
dev.off()
