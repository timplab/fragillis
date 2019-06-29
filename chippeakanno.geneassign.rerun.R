source("https://bioconductor.org/biocLite.R")
biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

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


                                        #load in annotation data and make granges object

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
annoData2 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
annoData2[1:2]


                                        #With this code, I am associating my peaks with genes by overlapping the peaks with genes or promoter regions of genes. Here, I define promoter regions as anything less than 1000bp upstream of a TSS. This specific code finds peaks located within promoters only... this is defined by the term bindingRegion=c(-1000,1), which means up to 1000bp upstream of the TSS and 1bp downstream of the TSS.


  
overlaps.unique24bft = annoPeaks(peaks.unique24bft.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.unique24bft.id = addGeneIDs(overlaps.unique24bft, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.unique24bft.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/24hrbft2openedpeaks.promoters.1000.1.chippeakanno.csv")

overlaps.unique24blank = annoPeaks(peaks.unique24blank.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.unique24blank.id = addGeneIDs(overlaps.unique24blank, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.unique24blank.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/24hrbft2closedpeaks.promoters.1000.1.chippeakanno.csv")

overlaps.unique48bft = annoPeaks(peaks.unique48bft.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.unique48bft.id = addGeneIDs(overlaps.unique48bft, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.unique48bft.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/48hrbft2openedpeaks.promoters.1000.1.chippeakanno.csv")

overlaps.unique48blank = annoPeaks(peaks.unique48blank.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.unique48blank.id = addGeneIDs(overlaps.unique48blank, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.unique48blank.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/48hrbft2closedpeaks.promoters.1000.1.chippeakanno.csv")




overlaps.24bft = annoPeaks(peaks.24bft.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.24bft.id = addGeneIDs(overlaps.24bft, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.24bft.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/24hrbft2peaks.promoters.1000.1.chippeakanno.csv")

overlaps.24blank = annoPeaks(peaks.24blank.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.24blank.id = addGeneIDs(overlaps.24blank, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.24blank.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/24hrblankpeaks.promoters.1000.1.chippeakanno.csv")

overlaps.48bft = annoPeaks(peaks.48bft.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.48bft.id = addGeneIDs(overlaps.48bft, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.48bft.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/48hrbft2peaks.promoters.1000.1.chippeakanno.csv")

overlaps.48blank = annoPeaks(peaks.48blank.gr, annoData=annoData2, bindingType="startSite", bindingRegion=c(-1000,1)) 
overlaps.48blank.id = addGeneIDs(overlaps.48blank, "org.Hs.eg.db", feature_id_type= "entrez_id", IDs2Add="symbol")
write.csv(overlaps.48blank.id, "~/Data/atac/final.analysis/rerun/chippeakanno_peaks/48hrblankpeaks.promoters.1000.1.chippeakanno.csv")



