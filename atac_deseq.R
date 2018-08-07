library("Rsubread")

outdir=("~/Dropbox/timplab_data/fragillis/180725_revision")

##From jawara, but heavily edited by Timp:
if (TRUE) {
    ##this code sums up the reads from all the bam files using an saf file,
    ##which I made from the narrowPeak (narrowpeaktosaf.sh)

    bampath="/atium/Data/NGS/Aligned/170423_atacseq_JA/"
    
    feature = featureCounts(annot.ext="/atium/Data/NGS/Aligned/170423_atacseq_JA/allsamples.peaks.qfilt.saf",
                            nthreads=10,
                            
                      files=file.path(bampath, c("aligned_031/24hrbft2.sorted.031.bam",
                                                 "aligned_031/48hrbft2.sorted.031.bam",
                                                 "aligned_031/24hrblank.sorted.031.bam",
                                                 "aligned_031/48hrblank.sorted.031.bam",
                                                 "aligned_031/time0.sorted.031.bam",
                                                 "aligned_220/24hrbft2.sorted.220.bam",
                                                 "aligned_220/48hrbft2.sorted.220.bam",
                                                 "aligned_220/24hrblank.sorted.220.bam",
                                                 "aligned_220/48hrblank.sorted.220.bam",
                                                 "aligned_220/time0.sorted.220.bam",
                                                 "aligned_trial1/24hrbft2.sorted.trial1.bam",
                                                 "aligned_trial1/48hrbft2.sorted.trial1.bam",
                                                 "aligned_trial1/24hrblank.sorted.trial1.bam",
                                                 "aligned_trial1/48hrblank.sorted.trial1.bam",
                                                 "aligned_trial1/time0.sorted.trial1.bam")),
                                      isGTFAnnotationFile=FALSE)

    save(feature, file=file.path(outdir, "features.rda"))
    
}






##Jawara's stuff (deseq.R)
#read in data
counts = read.csv("~/Code/atacseq/resubmission/featurecounts.allsamples.qfilt.txt", header=TRUE, sep= "\t")
counts2=counts
colnames(counts2) = c("GeneID","Length","24hrbft2.031","48hrbft2.031","24hrblank.031","48hrblank.031","time0.031","24hrbft2.220","48hrbft2.220","24hrblank.220","48hrblank.220","time0.220","24hrbft2.trial1","48hrbft2.trial1","24hrblank.trial1","48hrblank.trial1","time0.trial1")
counts2 = counts2[ ,3:17]

#separate data by 24 and 48 hours
counts.24hrs = counts2[ , c("24hrbft2.031","24hrblank.031","24hrbft2.220","24hrblank.220","24hrbft2.trial1","24hrblank.trial1")]
counts.48hrs = counts2[ , c("48hrbft2.031","48hrblank.031","48hrbft2.220","48hrblank.220","48hrbft2.trial1","48hrblank.trial1")]

#read in phenotype data and separate by 24 and 48 hours
coldata = read.csv("~/Code/atacseq/resubmission/featurecounts.phenotypedata.csv", header=TRUE)
row.names(coldata) = coldata$ids
coldata.24hrs = coldata[ c("24hrbft2.031","24hrblank.031","24hrbft2.220","24hrblank.220","24hrbft2.trial1","24hrblank.trial1"),]
coldata.48hrs = coldata[ c("48hrbft2.031","48hrblank.031","48hrbft2.220","48hrblank.220","48hrbft2.trial1","48hrblank.trial1"),]
coldata.24hrs$toxin <- relevel(coldata.24hrs$toxin, "blank")
coldata.48hrs$toxin <- relevel(coldata.48hrs$toxin, "blank")


###run Deseq2
 #source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)

dds.24 <- DESeqDataSetFromMatrix(countData = counts.24hrs, colData = coldata.24hrs, design = ~ toxin)
dds.48 <- DESeqDataSetFromMatrix(countData = counts.48hrs, colData = coldata.48hrs, design = ~ toxin)

data.24 = DESeq(dds.24)
results.24 = results(data.24)
res.ordered.24 = results.24[order(results.24$pvalue),]

data.48 = DESeq(dds.48)
results.48 = results(data.48)
res.ordered.48 = results.48[order(results.48$pvalue),]

#code used to correct padj values because histogram of pvalues looks like a hill.
#install.packages("fdrtool")
library("fdrtool")


FDR.DESeq2Res <- fdrtool(res.ordered.24$stat, statistic= "normal", plot = T)
res.ordered.24[,"padjcorrect"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")

res.sig.24 = res.ordered.24[res.ordered.24$pvalue < .01,]
res.summary.24 = summary(res.ordered.24)
sum(res.ordered.24$pvalue < .01, na.rm = TRUE)
 

FDR.DESeq2Res.48 <- fdrtool(res.ordered.48$stat, statistic= "normal", plot = T)
res.ordered.48[,"padjcorrect"]  <- p.adjust(FDR.DESeq2Res.48$pval, method = "BH")

res.sig.48 = res.ordered.48[res.ordered.48$pvalue < .01,]
res.summary.48 = summary(results.48)
sum(res.ordered.48$padjcorrect < .25, na.rm = TRUE)




                                        #write out final files

write.table(res.sig.48, file="~/Code/atacseq/resubmission/48hrpeaks.qfilt.deseq.csv")


write.table(res.sig.24, file="~/Code/atacseq/resubmission/24hrpeaks.qfilt.deseq.csv")


#this code will print pdf of the pvalues
pdf(file="hist.48.corrected.pdf")
hist(FDR.DESeq2Res.48$pval, xlim = c(0,1))
dev.off()


pdf(file="hist.24.corrected.pdf")
hist(FDR.DESeq2Res$pval, xlim = c(0,1))
dev.off()

pdf(file="hist.48.pdf")
hist(res.ordered.48$pvalue, xlim = c(0,1))
dev.off()




    ##write.table(x=data.frame(b$annotation[,c("GeneID","Length")],b$counts,stringsAsFactors=FALSE),file="featurecounts.allsamples.qfilt.txt",quote=FALSE,sep="\t")

