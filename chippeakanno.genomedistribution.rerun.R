source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)

#load in my data and make granges object
peaks.unique24bft= read.table(file="~/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrbft2.narrowPeak",sep="\t",head=F,stringsAsFactors=F)
peaks.unique24bft.gr=GRanges(seqnames=peaks.unique24bft[,1],IRanges(start=peaks.unique24bft[,2],end=peaks.unique24bft[,3]),score=peaks.unique24bft[,5])
peaks.unique24bft.gr$score <- as.numeric(peaks.unique24bft.gr$score)

peaks.unique24blank= read.table(file="~/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrblank.narrowPeak",sep="\t",head=F,stringsAsFactors=F)
peaks.unique24blank.gr=GRanges(seqnames=peaks.unique24blank[,1],IRanges(start=peaks.unique24blank[,2],end=peaks.unique24blank[,3]),score=peaks.unique24blank[,5])
peaks.unique24blank.gr$score <- as.numeric(peaks.unique24blank.gr$score)

peaks.unique48bft= read.table(file="~/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrbft2.narrowPeak",sep="\t",head=F,stringsAsFactors=F)
peaks.unique48bft.gr=GRanges(seqnames=peaks.unique48bft[,1],IRanges(start=peaks.unique48bft[,2],end=peaks.unique48bft[,3]),score=peaks.unique48bft[,5])
peaks.unique48bft.gr$score <- as.numeric(peaks.unique48bft.gr$score)

peaks.unique48blank= read.table(file="~/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrblank.narrowPeak",sep="\t",head=F,stringsAsFactors=F)
peaks.unique48blank.gr=GRanges(seqnames=peaks.unique48blank[,1],IRanges(start=peaks.unique48blank[,2],end=peaks.unique48blank[,3]),score=peaks.unique48blank[,5])
peaks.unique48blank.gr$score <- as.numeric(peaks.unique48blank.gr$score)



peaks.24bft= read.table(file="/home/jallen66/atac_dnase_pipelines/rerun/out.24hrbft2/peak/macs2/idr/optimal_set/out.24hrbft2_ppr.IDR0.1.filt.narrowPeak.gz",sep="\t",head=F,stringsAsFactors=F)
peaks.24bft.gr=GRanges(seqnames=peaks.24bft[,1],IRanges(start=peaks.24bft[,2],end=peaks.24bft[,3]),score=peaks.24bft[,5])
peaks.24bft.gr$score <- as.numeric(peaks.24bft.gr$score)

peaks.24blank= read.table(file="/home/jallen66/atac_dnase_pipelines/rerun/out.24hrblank/peak/macs2/idr/optimal_set/out.24hrblank_ppr.IDR0.1.filt.narrowPeak.gz",sep="\t",head=F,stringsAsFactors=F)
peaks.24blank.gr=GRanges(seqnames=peaks.24blank[,1],IRanges(start=peaks.24blank[,2],end=peaks.24blank[,3]),score=peaks.24blank[,5])
peaks.24blank.gr$score <- as.numeric(peaks.24blank.gr$score)

peaks.48bft= read.table(file="/home/jallen66/atac_dnase_pipelines/rerun/out.48hrbft2/peak/macs2/idr/optimal_set/out.48hrbft2_ppr.IDR0.1.filt.narrowPeak.gz",sep="\t",head=F,stringsAsFactors=F)
peaks.48bft.gr=GRanges(seqnames=peaks.48bft[,1],IRanges(start=peaks.48bft[,2],end=peaks.48bft[,3]),score=peaks.48bft[,5])
peaks.48bft.gr$score <- as.numeric(peaks.48bft.gr$score)

peaks.48blank= read.table(file="/home/jallen66/atac_dnase_pipelines/rerun/out.48hrblank/peak/macs2/idr/optimal_set/out.48hrblank_ppr.IDR0.1.filt.narrowPeak.gz",sep="\t",head=F,stringsAsFactors=F)
peaks.48blank.gr=GRanges(seqnames=peaks.48blank[,1],IRanges(start=peaks.48blank[,2],end=peaks.48blank[,3]),score=peaks.48blank[,5])
peaks.48blank.gr$score <- as.numeric(peaks.48blank.gr$score)




#here, I am determing which region of the genome each peak is in. I am then creating an ouput table with this info foreach sample. I am using that table for the next step of the analysis.

biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

aCR.unique24bft <- assignChromosomeRegion(peaks.unique24bft.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.unique24blank <- assignChromosomeRegion(peaks.unique24blank.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.unique48bft <- assignChromosomeRegion(peaks.unique48bft.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.unique48blank <- assignChromosomeRegion(peaks.unique48blank.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)

aCR.24bft <- assignChromosomeRegion(peaks.24bft.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.24blank <- assignChromosomeRegion(peaks.24blank.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.48bft <- assignChromosomeRegion(peaks.48bft.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aCR.48blank <- assignChromosomeRegion(peaks.48blank.gr, nucleotideLevel=FALSE, precedence= c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Introns","Exons"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)


write.table(aCR.unique24bft, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbftopened", sep="\t")
write.table(aCR.unique24blank, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbftclosed", sep="\t")
write.table(aCR.unique48bft, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbftopened", sep="\t")
write.table(aCR.unique48blank, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbftclosed", sep="\t")

write.table(aCR.24bft, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbft2", sep="\t")
write.table(aCR.24blank, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrblank", sep="\t")
write.table(aCR.48bft, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbft2", sep="\t")
write.table(aCR.48blank, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrblank", sep="\t")


#now, I am creating objects for peak genome distribution for each sample and then merging them and printing them. For my input file, I am reading in the files that I create in the previous step. Basically, I am reading in each file with info about the location of peaks for each sample, and then creating one large data table with all of the data. I then put that data table into excel in order to makea  barchart, but I am sure I can make the same stacked bar chart with ggplot or some other R application


dis.unique24bft = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbftopened")
colnames(dis.unique24bft) = c("region", "percentage.unique24bft",'jaccardindex.unique24bft')

dis.unique24blank = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbftclosed")
colnames(dis.unique24blank) = c("region", "percentage.unique24blank",'jaccardindex.unique24blank')

dis.unique48bft = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbftopened")
colnames(dis.unique48bft) = c("region", "percentage.unique48bft",'jaccardindex.unique48bft')

dis.unique48blank = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbftclosed")
colnames(dis.unique48blank) = c("region", "percentage.unique48blank",'jaccardindex.unique48blank')



dis.24bft = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrbft2")
colnames(dis.24bft) = c("region", "percentage.24bft","jaccardindex.24bft")

dis.48bft = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrbft2")
colnames(dis.48bft) = c("region", "percentage.48bft","jaccardindex.48bft")

dis.24blank = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_24hrblank")
colnames(dis.24blank) = c("region", "percentage.24blank","jaccardindex.24blank")

dis.48blank = read.table("/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_48hrblank")
colnames(dis.48blank) = c("region", "percentage.48blank","jaccardindex.48blank")




peak.distribution = cbind(dis.unique24bft,dis.unique24blank, dis.unique48bft, dis.unique48blank, dis.24bft, dis.24blank, dis.48bft, dis.48blank)
peak.distribution = subset(peak.distribution, select = -c(region))
peak.distribution= subset(peak.distribution, select = -c(region.1,region.2,region.3,region.4,region.5,region.6))
write.table(peak.distribution, "/home/jallen66/Data/atac/final.analysis/rerun/chippeakanno_peaks/peak_genome_distribution_allsamples.csv", sep=",",row.names=F)
