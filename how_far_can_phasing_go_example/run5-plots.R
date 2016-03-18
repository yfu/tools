

a <- read.table("all.z1")

watson <- a[a$V3 == "Watson", ]
crick <- a[a$V3 == "Crick", ]

library(ggplot2)

pdf("how_far_can_phasing_go.pdf", useDingbats=FALSE)
my.title <- "Watson strand \nnegative == upstream of insertion\npositive == downstream of insertion"
ggplot(watson, aes(x=V2, y=V4, group=V1, color=V1)) + geom_point() + geom_line() + ggtitle(my.title) + ylab("z1 score") + xlab("offset")
my.title <- "Crick strand \nnegative == upstream of insertion\npositive == downstream of insertion"
ggplot(crick, aes(x=V2, y=V4, group=V1, color=V1)) + geom_point() + geom_line() + ggtitle(my.title) + ylab("z1 score") + xlab("offset")
dev.off()
