a <- read.table("Mattick.RNaseR.HeLa.BP.rep1.fq.gz.map_to_put_lar.error_rate")
colnames(a) <- c("pos.on.lar", "n.covered", "n.err")
a$err.rate <- a[, 3] / a[, 2]
pdf("error_rate.pdf")
plot(a$pos.on.lar, a$err.rate, xlab="position on lariat", ylab="error rate", main="Error rates calculated from lariat supporting reads")
dev.off()
