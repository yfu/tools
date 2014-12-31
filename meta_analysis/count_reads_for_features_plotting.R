library(ggplot2)

# For gene plots, I need to have upstream and downstream regions in them, so they are treated separately
pdf("meta.all.pdf", width=10)

for (i in c("00dpp", "02dpp", "04dpp", "07dpp", "10dpp", "14dpp", "17dpp", "20dpp", "42dpp")) {
    for (j in c("Genes")) {
        basic <- read.table(paste(i, ".basic_stats", sep=""), sep="\t")
        nf <- basic[5, 2]
        a <- read.table(paste(i, ".", j, ".aggregation.piRNA.100bin.signal", sep=""))
        a$V3 <- a$V3 / nf * 1e6
        idx <- a$V1=="Antisense"
        a[idx, 3] <- -a[idx, 3]
        a$type <- "gene"

        df <- a
        u.fn <- paste(i, ".", j, "_Upstream10k",".aggregation.piRNA.100bin.signal", sep="")
        d.fn <- paste(i, ".", j, "_Downstream10k",".aggregation.piRNA.100bin.signal", sep="")
        u.exists = file.exists(u.fn)
        d.exists = file.exists(d.fn)
        if (u.exists == TRUE && d.exists == TRUE) {

            u <- tryCatch(read.table(paste(i, ".", j, "_Upstream10k",".aggregation.piRNA.100bin.signal", sep="")),
                          error = function(e) {
                              warning(paste(u.fn, " is an empty file. I will not plot upstream10k."))
                              0
                          }
                          )
            if (length(u)==1) {
                df <- df
            } else {
                u$V3 <- u$V3 / nf * 1e6
                u$V2 <- u$V2 - 100
                idx <- u$V1=="Antisense"
                u[idx, 3] <- -u[idx, 3]
                u$type <- "upstream10k"
                df <- rbind(df, u)
            }
            
            d <- tryCatch(read.table(paste(i, ".", j, "_Downstream10k",".aggregation.piRNA.100bin.signal", sep="")),
                          error = function(e) {
                              warning(paste(u.fn, " is an empty file. I will not plot upstream10k."))
                              0
                          }
                          )
            if (length(d)==1) {
                df <- df
            } else {
                d$V3 <- d$V3 / nf * 1e6
                d$V2 <- d$V2 + 100
                idx <- d$V1=="Antisense"
                d[idx, 3] <- -d[idx, 3]
                d$type <- "downstream10k"
                df <- rbind(df, d)
            }
        }
        p <- ggplot(df, aes(x=V2, y=V3, group=type, color=V1)) + geom_point() + ggtitle(paste(i, j)) + xlab("[-100:0) 5'flanking, [0-99) gene body, [99, 199) 3'flanking") + ylab("# 5'ends of small RNA reads per million genome mapping reads (-rRNA; +miRNA_hairpin)")
        print(p)
    }
}

## dev.off()

## pdf("meta_Exons.pdf", width=10)

    ## for (j in c("Genes")) {
for (j in c("Genes_Exons", "Genes_Exons_5UTR", "Genes_Exons_3UTR", "Genes_Intron")) {
    for (i in c("00dpp", "02dpp", "04dpp", "07dpp", "10dpp", "14dpp", "17dpp", "20dpp", "42dpp")) {
        basic <- read.table(paste(i, ".basic_stats", sep=""), sep="\t")
        nf <- basic[5, 2]
        
        a <- read.table(paste(i, ".", j, ".aggregation.piRNA.100bin.signal", sep=""))
        a$V3 <- a$V3 / nf * 1e6
        idx <- a$V1=="Antisense"
        a[idx, 3] <- -a[idx, 3]
        a$type <- "gene"
        
        df <- a
        
        p <- ggplot(df, aes(x=V2, y=V3, group=type, color=V1)) + geom_point() + ggtitle(paste(i, j)) + xlab("[-100:0) 5'flanking, [0-99) gene body, [99, 199) 3'flanking") + ylab("# 5'ends of small RNA reads per million genome mapping reads (-rRNA; +miRNA_hairpin)") 
        
        print(p)
        ## print(p + geom_line())
    }
}

dev.off()
