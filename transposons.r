library(DESeq2)
library(gplots)
library(ggplot2)

setwd("/home/fuy2/data/RS_raw/transposon_htseq_differential_analysis/RSQ_masterTable_All")

# suf = "FLY_TRANSPOSON_ALL.htseqcountMult.AS"; dir = "transposon_as"; processing.transposon <- TRUE
# suf = "FLY_TRANSPOSON_ALL.htseqcountMult.S"; dir = "transposon_s"; processing.transposon <- TRUE

# suf = "FLY_PIRNA_CLUSTER.htseqcountMult.AS"; dir = "pirnacluster_as"; processing.transposon <-FALSE

# suf = "FLY_PIRNA_CLUSTER.htseqcountMult.S"; dir = "pirnacluster_s"; processing.transposon <-FALSE

# suf = "FLY_GENE_TRANSPOSON_ALL.htseqcountMult.AS"; dir = "gene_transposon_as"; processing.transposon <-FALSE
suf <- "FLY_TRANSPOSON_ALL_GENE.S"; dir <- "transposon_s"; processing.transposon <- TRUE
# suf = "FLY_GENE_TRANSPOSON_ALL.htseqcountMult.S"; dir = "gene_transposon_s"; processing.transposon <-FALSE

file <-paste("RSQ_masterTable.", suf, ".raw.txt", sep="")

htv <- read.table(file,T)

htv$feature=sub("FBgn0000004_17","FBgn0000004_17.6",htv$feature)
countsTable <- htv[-1] #remove the feature column
featureName=unlist(htv[1]) #get the row names from the first column

featureName=sub('FBgn[0-9n]+\\_','',featureName,perl=TRUE)
featureName=sub("FBgnnnnnnnn_","",featureName,perl=TRUE) 
rownames(countsTable)=featureName

# trn=read.table("trn.list",F)
# colnames(trn)=c("feature")
# trn$feature=sub('FBgn[0-9n]+\\_','',trn$feature,perl=TRUE)
# trn$feature=sub("FBgnnnnnnnn_","",trn$feature,perl=TRUE)

roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# Remove those unused barcodes (or experiments, so to speak)
unused <- c("HEIL_id1", "HEIL_id2", # These two are for GFP project
            "HEIL_id5", "HEIL_id6",
            "JMAC_id1", "JMAC_id2", 
            "KNBD_id1", "KNBD_id2", 
            "OQSU_id5", "OQSU_id6", 
            "PRTV_id5", "PRTV_id6",
            "X10131619_id5", "X10131619_id6",
            "X12pmpmmp_id1", "X12pmpmmp_id4",
            "X12pmpmmp_id2", "X12pmpmmp_id3" # These are for GFP project
            )

# Drop the the unused barcodes
# Not sure if this is a bug of R, but when you do df[,-somevector], if somevector is integer(0), then the whole data frame is dropped.
if ( length(which(colnames(countsTable) %in% unused) )  != 0 ) {
    countsTable <- countsTable[, -which(colnames(countsTable) %in% unused)]
}    

row_sub = apply(countsTable, 1, function(row) all(row !=0 ))
countsTable=countsTable[row_sub,]

numofRows=nrow(countsTable)

GT=factor(c(# "GFP: w[1]; +/+; nos>EGFP/+ (F1) BR1", "GFP: w[1]; GS1/+; nos>EGFP/+ (F1) BR1",
    # Now it is 10131619's turn
    "Testis: ago3[t2]/TM6B & ago3[t3]/TM6B BR1", "Testis: ago3[t2]/ago3[t3] BR1",
    "Testis: aub[HN2]/CyO & aub[QC42]/CyO BR1", "Testis: aub[HN2]/aub[QC42] BR1",
    # Now it is 12pmpmmp's turn
    "Testis: piwi[2]/CyO BR1", "Testis: piwi[2]/piwi[2] BR1",
    # HEIL starts here
    "Testis: XO; +/+; +/ry[506] BR1", "Testis: Oregon R BR1",
    # JMAC begins here
    "Testis: XO; +/+; +/ry[506] BR2", "Testis: Oregon R BR2",
    "Testis: piwi[2]/CyO BR2", "Testis: piwi[2]/piwi[2] BR2",
    # KNBD starts here
    "Testis: XO; +/+; +/ry[506] BR3", "Testis: Oregon R BR3",
    "Testis: piwi[2]/CyO BR3", "Testis: piwi[2]/piwi[2] BR3",
    # Here comes OQSU
    "Testis: ago3[t2]/TM6B & ago3[t3]/TM6B BR2", "Testis: ago3[t2]/ago3[t3] BR2",
    "Testis: aub[HN2]/CyO & aub[QC42]/CyO BR2", "Testis: aub[HN2]/aub[QC42] BR2",
    # Now it is PRTV's turn
    "Testis: ago3[t2]/TM6B & ago3[t3]/TM6B BR3", "Testis: ago3[t2]/ago3[t3] BR3",
    "Testis: aub[HN2]/CyO & aub[QC42]/CyO BR3", "Testis: aub[HN2]/aub[QC42]"))


replicates = factor(c(
    # 10131619
    "ago3Het", "ago3Mut",
    "aubHet", "aubMut",
    # 12pmpmmp
    "piwiHet", "piwiMut",
    # HEIL
    "XO", "WT",
    # JMAC
    "XO", "WT",
    "piwiHet", "piwiMut",
    # KNBD
    "XO", "WT",
    "piwiHet", "piwiMut",
    # OQSU
    "ago3Het", "ago3Mut",
    "aubHet", "aubMut",
    # PRTV
    "ago3Het", "ago3Mut",
    "aub3Het", "aubMut"
))

type=rep("PE", length(GT))
colData=data.frame(GT, type, replicates)

# Is this CORRECT?
cc <- apply(countsTable, c(1,2), round)
# dds <- DESeqDataSetFromMatrix(countData=countsTable,colData= colData, design=~ GT)
dds <- DESeqDataSetFromMatrix(countData=cc, colData= colData, design=~ replicates)

# Use some sample as the reference sample
dds1=dds
# WT is 
colData(dds1)$replicates <- relevel (colData(dds1)$replicates,"XO")

# colData(dds1)$GT <- relevel (colData(dds1)$GT, "Testis: piwi[2]/CyO BR2", "Testis: piwi[2]/piwi[2] BR2")

# dds1 <- DESeq (dds1)

dds1 <- estimateSizeFactors(dds1)

dds1 <- estimateDispersions(dds1)

# ddsw1 <- nbinomWaldTest(ddsw1,cooksCutoff=FALSE)
dds1 <- nbinomWaldTest(dds1)

countsHM=counts(dds1,normalized=TRUE)

colnames(countsHM)=colnames(cc)

outfilename=paste("my.", suf,".normalized.counts.txt",sep="")
write.table(countsHM,outfilename,row.names=TRUE)

#pdfname=paste("RSQ.allgtwoZZ.hcluster.heatmap.normalizedcountstow1.", "dummy",".all.pdf",sep="")
#pdf(pdfname,width=10,height=20)

# heatmap.2(log2(countsHM+1))


                                        #col = hmcol,
#                          Rowv = TRUE, Colv = TRUE, scale="column",
#                          dendrogram="both", trace="none", margin=c(13,10))
#dev.off()

trn.group <- read.table("Zamore.group")
colnames(trn.group)=c("feature", "group")
trn.group$feature=sub('FBgn[0-9n]+\\_','',trn.group$feature,perl=TRUE)
trn.group$feature=sub("FBgnnnnnnnn_","",trn.group$feature,perl=TRUE)
rownames(trn.group) = trn.group[[1]]

group <- data.frame(trn.group[,2])
colnames(group) = c("group")

group$group <- as.character(group$group)

rownames(group) = trn.group[[1]]

if (processing.transposon == TRUE){
    gene.rows <- sapply(row.names(countsHM), function(x) grepl("CG", x))
    counts.with.group <- countsHM[!gene.rows,] 
    counts.with.group <- merge(counts.with.group, group, by="row.names", all.x=TRUE)
    counts.with.group$group[is.na(counts.with.group$group)] = "unknown"
    # rownames(counts.with.group) = counts.with.group$Row.names
    # counts.with.group <- counts.with.group[, -which(colnames(counts.with.group) %in% c("Row.names"))]
}
## else {
##     counts.with.group <- data.frame(countsHM)
##     counts.with.group["group"] = NA
##     counts.with.group["Row.names"] = rownames(counts.with.group)
##     # Rearrange the columns... R is such a hassle...
##     counts.with.group = counts.with.group[c(ncol(counts.with.group), seq(1, ncol(counts.with.group) - 1))]
##     # counts.with.group$group <- list(rep(NA, times = nrow(countsHM)))
## }

library(ggplot2)

levels(counts.with.group) = c("1","2","3", "unknown")

# "HEIL <- id3"  "HEIL <- id4"  "JMAC <- id3"  "JMAC <- id4"  "JMAC <- id5"
# "JMAC <- id6"  "KNBD <- id3"  "KNBD <- id4"  "KNBD <- id5"  "KNBD <- id6"  "OQSU <- id1"
# "OQSU <- id2"  "OQSU <- id3"  "OQSU <- id4"  "PRTV <- id1"  "PRTV <- id2"  "PRTV <- id3"
# "PRTV <- id4"


transposon.plot <- function(interesting, plot.name, dir, lim) {
    for (i in 1:ncol(interesting)) {
        x.sample = interesting[,i][1]
        x.no = which(colnames(counts.with.group) == x.sample) - 1 # Get the index in GT, and get xlab
        x.label = as.character(GT[x.no])
        print(c("Sample on the x-axis:", x.sample, x.label))
        print(x.no)
        y.sample = interesting[,i][2]
        y.no = which(colnames(counts.with.group) == y.sample) - 1
        y.label = as.character(GT[y.no])
        print(c("Sample on the y-axis", y.sample, y.label))
        print(y.no)
        if (processing.transposon == TRUE) {
            p = ggplot(data=counts.with.group, aes_string(x=x.sample, y=y.sample, col="group"))
        } else {
            # If it is not transposon, we do not need the group info
            p = ggplot(data=counts.with.group, aes_string(x=x.sample, y=y.sample))
        }
        p = p + xlim(lim) + ylim(lim)

        # Unlabeled
        pdf(paste(dir, "/", plot.name, "_unlabeled_", x.sample, "_", y.sample ,".pdf", sep="") )
        print(p + geom_point() + geom_abline(intercept=0, slope=1) + xlab(x.label) + ylab(y.label))
    dev.off()
        # Labeled
        pdf(paste(dir, "/", plot.name, "_labeled_", x.sample, "_",  y.sample ,".pdf", sep=""))
        print(p + geom_point() + geom_abline(intercept=0, slope=1) + geom_text(aes(label=Row.names),hjust=0, vjust=0, size=2) + xlab(x.label) + ylab(y.label))
        dev.off()
    }
}

interesting.pairs.br = matrix( # br: biological replicates
    c("X10131619_id1", "OQSU_id1", # Ago3 Het
      "X10131619_id2", "OQSU_id2", # Ago3 Mut
      "X10131619_id3", "OQSU_id3", # Aub Het
      "X10131619_id4", "OQSU_id4", # Aub Mut
      "X12pmpmmp_id5", "JMAC_id5", # piwi Het
      "X12pmpmmp_id6", "JMAC_id6", # piwi Mut
      "HEIL_id3", "JMAC_id3", # XO
      "HEIL_id4", "JMAC_id4", # WT
      "OQSU_id1", "PRTV_id1", # Ago3 Het
      "OQSU_id2", "PRTV_id2", # Ago3 Het
      "OQSU_id3", "PRTV_id3", # Aub Het
      "OQSU_id4", "PRTV_id4", # Aub Mut
      "JMAC_id5", "KNBD_id5", # piwi Het
      "JMAC_id6", "KNBD_id6", # piwi Mut
      "JMAC_id3", "KNBD_id3", # XO
      "JMAC_id4", "KNBD_id4"  # WT
    ), nrow=2)

dir.create(dir)

transposon.plot(interesting.pairs.br, paste("between_replicates_", suf, sep=""), dir=dir, lim=c(0, 300))

interesting.pairs = matrix(
    c("X10131619_id1", "X10131619_id2", # Ago3 Het vs Ago3 Mut BR1
      "X10131619_id3", "X10131619_id4", # Aub Het vs Aub Mut BR1
      "X12pmpmmp_id5", "X12pmpmmp_id6", # piwi Het vs piwi Mut BR1
      "HEIL_id3", "HEIL_id4", # XO vs WT BR1

      "OQSU_id1", "OQSU_id2", # Ago3 Het vs Ago3 Mut BR2
      "OQSU_id3", "OQSU_id4", # Aub Het vs Aub Mut BR
      "JMAC_id5", "JMAC_id6", # piwi Het vs piwi Mut BR2
      "JMAC_id3", "JMAC_id4", # XO vs WT BR2
      
      "PRTV_id1", "PRTV_id2", # Ago3 Het vs Ago3 Mut BR3
      "PRTV_id3", "PRTV_id4", # Aub Het vs Aub Mut BR3
      "KNBD_id5", "KNBD_id6", # piwi Het vs piwi Mut BR3
      "KNBD_id3", "KNBD_id4", # XO vs WT BR3
      
      "OQSU_id1", "JMAC_id4", # ago3 Het vs WT
      "OQSU_id2", "JMAC_id4", # ago3 Mut vs WT
      "OQSU_id3", "JMAC_id4", # Aub Het vs WT
      "OQSU_id4", "JMAC_id4" # Aub Mut vs WT
      ), nrow = 2)

transposon.plot(interesting.pairs, paste("between_het_mut_wt_", suf, sep=""), dir=dir, lim=c(0, 300))

# Now is the time to calculate all kinds of p-values
replicates = factor(c(
    # 10131619
    "ago3Het", "ago3Mut",
    "aubHet", "aubMut",
    # 12pmpmmp
    "piwiHet", "piwiMut",
    # HEIL
    "XO", "WT",
    # JMAC
    "XO", "WT",
    "piwiHet", "piwiMut",
    # KNBD
    "XO", "WT",
    "piwiHet", "piwiMut",
    # OQSU
    "ago3Het", "ago3Mut",
    "aubHet", "aubMut",
    # PRTV
    "ago3Het", "ago3Mut",
    "aub3Het", "aubMut"
))

dds2 <- dds
colData(dds2)$replicates = factor(as.vector(replicates), levels=c("WT", "ago3Het", "ago3Mut", "aub3Het",
                                                             "aubHet", "aubMut", "piwiHet", "piwiMut", "XO"))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)
# ddsw1 <- nbinomWaldTest(ddsw1,cooksCutoff=FALSE)
dds2 <- nbinomWaldTest(dds2)
# Now calculate the p-values
res <- results(dds2)
res <- res[order(res$padj), ]
res <- as.matrix(res)
gene.rows <- sapply(row.names(res), function(x) grepl("CG", x))
res <- res[!gene.rows,] 
write.table(as.matrix(res), file="transposons.XO.v.WT.txt", quote=FALSE, sep="\t", col.names=NA)


dds2 <- dds
colData(dds2)$replicates = factor(as.vector(replicates), levels=c("ago3Het", "WT", "aub3Het",
                                                             "aubHet", "aubMut", "piwiHet", "piwiMut", "XO", "ago3Mut"))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)
# ddsw1 <- nbinomWaldTest(ddsw1,cooksCutoff=FALSE)
dds2 <- nbinomWaldTest(dds2)
# Now calculate the p-values
res <- results(dds2)
res <- res[order(res$padj), ]
res <- as.matrix(res)
gene.rows <- sapply(row.names(res), function(x) grepl("CG", x))
res <- res[!gene.rows,] 
write.table(as.matrix(res), file="transposons.ago3Het.v.ago3Mut.txt", quote=FALSE, sep="\t", col.names=NA)

dds2 <- dds
colData(dds2)$replicates = factor(as.vector(replicates), levels=c("aubHet", "ago3Het", "WT", "aub3Het",
                                                               "piwiHet", "piwiMut", "XO", "ago3Mut", "aubMut"))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)
# ddsw1 <- nbinomWaldTest(ddsw1,cooksCutoff=FALSE)
dds2 <- nbinomWaldTest(dds2)
# Now calculate the p-values
res <- results(dds2)
res <- res[order(res$padj), ]
res <- as.matrix(res)
gene.rows <- sapply(row.names(res), function(x) grepl("CG", x))
res <- res[!gene.rows,] 
write.table(as.matrix(res), file="transposons.aubHet.v.aubMut.txt", quote=FALSE, sep="\t", col.names=NA)

dds2 <- dds
colData(dds2)$replicates = factor(as.vector(replicates), levels=c("piwiHet", "aubHet", "ago3Het", "WT", "aub3Het",
                                                                "XO", "ago3Mut", "aubMut", "piwiMut"))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)
# ddsw1 <- nbinomWaldTest(ddsw1,cooksCutoff=FALSE)
dds2 <- nbinomWaldTest(dds2)
# Now calculate the p-values
res <- results(dds2)
res <- res[order(res$padj), ]
res <- as.matrix(res)
gene.rows <- sapply(row.names(res), function(x) grepl("CG", x))
res <- res[!gene.rows,] 
write.table(as.matrix(res), file="transposons.piwiHet.v.piwiMut.txt", quote=FALSE, sep="\t", col.names=NA)
