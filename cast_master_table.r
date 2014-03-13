library(plyr)
library(reshape)



#cast_master_table <- function (input,outfilename) {

args <- commandArgs (TRUE);
input<-args[1]
outfilename<- args[2]
a=read.table(input,F)
colnames(a)=c("gt","feature","count")
b=cast(a,feature~gt)

write.table(b,outfilename,append = FALSE, quote = FALSE, sep = "\t",
			eol = "\n", na = "NA", dec = ".", row.names = FALSE,
			col.names = TRUE, qmethod = c("escape", "double"),
			fileEncoding = "")
			

			

#}

