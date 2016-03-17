source("~/repo/tools/stylish_scatterplot.R")

pdf("transpoons.SAS.pdf", useDingbats=FALSE)
for (i in list.files(".", "*vs_transposon.count$")) {
    ## a <- read.table("Zamore.RSQ.fly.ago3_mut_aub_mut.rep1.SRR1924207.sorted.unique.vs_transposon.count")
    a <- read.table(i)
    S <- a[a$V2=="sense", ]
    AS <- a[a$V2=="antisense", ]
    b <- merge(S, AS, by=c("V1", "V3"))
    

    ## my.max <- max(a$V4)
    my.max <- 1e4
    p <- stylish_scatterplot(b, "V5.x", "V5.y", "V3", "sense", "antisense", i, my.max)
    print(p)
}
dev.off()
