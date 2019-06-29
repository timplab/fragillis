install.packages("bedr")
library("bedr")
library("VennDiagram")


#first, read in peak data for unique peak files and individual sample peak files
sample.peak.24hrbft2= read.csv("/home/jallen66/atac_dnase_pipelines/rerun/out.24hrbft2/peak/macs2/idr/optimal_set/out.24hrbft2_ppr.IDR0.1.filt.narrowPeak.gz", sep="\t",header=FALSE)
sample.peak.24hrblank= read.csv("/home/jallen66/atac_dnase_pipelines/rerun/out.24hrblank/peak/macs2/idr/optimal_set/out.24hrblank_ppr.IDR0.1.filt.narrowPeak.gz", sep="\t",header=FALSE)
sample.peak.48hrbft2= read.csv("/home/jallen66/atac_dnase_pipelines/rerun/out.48hrbft2/peak/macs2/idr/optimal_set/out.48hrbft2_ppr.IDR0.1.filt.narrowPeak.gz", sep="\t",header=FALSE)
sample.peak.48hrblank= read.csv("/home/jallen66/atac_dnase_pipelines/rerun/out.48hrblank/peak/macs2/idr/optimal_set/out.48hrblank_ppr.IDR0.1.filt.narrowPeak.gz", sep="\t",header=FALSE)



unique.peak.24hropened = read.csv("/home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrbft2.narrowPeak", sep="\t", header=FALSE)
unique.peak.24hrclosed = read.csv("/home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrblank.narrowPeak", sep="\t", header=FALSE)
unique.peak.48hropened = read.csv("/home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrbft2.narrowPeak", sep="\t", header=FALSE)
unique.peak.48hrclosed = read.csv("/home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrblank.narrowPeak", sep="\t", header=FALSE)
 

crc = c("AKT1", "AKT2","AKT3","APC","APC2","APPL1","ARAF","Axin1","AXIN2","BAD","BAX","BCL2","BIRC5","BRAF","casp3","CASP9","CCND1","CTNNB1","CYCS","DCC","FOS","GSK3B","JUN","KRAS","LEF1","MAP2K1","MAPK1","MAPK10","MAPK3","MAPK8","MAPK9","MLH1","MSH2","MSH3","MSH6","MYC","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","RAC1","RAC2","RAC3","RAF1","RALGDS","RHOA","SMAD2","SMAD3","SMAD4","TCF7","TCF7L1","TCF7L2","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TP53")

wnt = c("APC","ASCL1","ASCL2","AXIN1","AXIN2","CTNNB1","BCL9","BCL92","ABL1","CSNK1A1","CSNK1D","DIXDC1","CREB","DVL1","DVL2","DVL3","DRAXIN","GSK3A","GSK3B","TGFB1I1","CTNNBIP1","PYGO1","PYGO2","HNF1B","UBE2A","TCF7","TCF7L1","TLE2")


bcat = c("MYC","MYCN","CCND1","HNF1A","LEF1","PPARD","JUN","FOSL1","PLAUR","MMP7","AXIN2","NRCAM","TCF4","GAST","CD44","EFNB1","EFNB2","BMP4","CLDN1","BIRC5","VEGFA","FGF18","ATOH1","MET","EDN1","MYCBP","L1CAM","ID2","JAG1","MSL1","TIAM1","NOS2","TERT","DKK1","FGF9","LBH","FGF20","LGR5","SOX9","SOX17","RUNX2","GREM1","SALL4","TNFSF11","TNFRSF11B","CYR61","PTTG1","DLL1","FOXN1","MMP26","NANOG","POU5F1","SNAI1","FN1","FZD7","ISL1","MMP2","MMP9","FST","WNT3A","TWIST1","TBX3","GBX2","CACNA1G","CDC25","WISP1","WISP2","IGF2","EMP1","IGF1","VEGFC","PTGS2","IL6","PITX2","EGFR","CDH1","CDKN2A","CTLA4","CXCL8","VCAN","TNFRSF19") # INFO FROM https://web.stanford.edu/group/nusselab/cgi-bin/wnt/target_genes



###this code will give me the percent of my peaks that overlap with eachother
require(GenomicRanges)

sample.peak.24hrbft2.gr=GRanges(seqnames=sample.peak.24hrbft2[,1],IRanges(start=sample.peak.24hrbft2[,2],end=sample.peak.24hrbft2[,3]))
sample.peak.24hrblank.gr=GRanges(seqnames=sample.peak.24hrblank[,1],IRanges(start=sample.peak.24hrblank[,2],end=sample.peak.24hrblank[,3]))
sample.peak.48hrbft2.gr=GRanges(seqnames=sample.peak.48hrbft2[,1],IRanges(start=sample.peak.48hrbft2[,2],end=sample.peak.48hrbft2[,3]))
sample.peak.48hrblank.gr=GRanges(seqnames=sample.peak.48hrblank[,1],IRanges(start=sample.peak.48hrblank[,2],end=sample.peak.48hrblank[,3]))


unique.peak.24hropened.gr=GRanges(seqnames=unique.peak.24hropened[,1],IRanges(start=unique.peak.24hropened[,2],end=unique.peak.24hropened[,3]))
unique.peak.24hrclosed.gr=GRanges(seqnames=unique.peak.24hrclosed[,1],IRanges(start=unique.peak.24hrclosed[,2],end=unique.peak.24hrclosed[,3]))
unique.peak.48hropened.gr=GRanges(seqnames=unique.peak.48hropened[,1],IRanges(start=unique.peak.48hropened[,2],end=unique.peak.48hropened[,3]))
unique.peak.48hrclosed.gr=GRanges(seqnames=unique.peak.48hrclosed[,1],IRanges(start=unique.peak.48hrclosed[,2],end=unique.peak.48hrclosed[,3]))



#code for overlaps
overlaps.24hrbft2.48hrbft2= as.data.frame(findOverlaps(sample.peak.24hrbft2.gr, sample.peak.48hrbft2.gr))
overlaps.list.24hrbft2.48hrbft2 = overlaps.24hrbft2.48hrbft2 [ ,"subjectHits"]
overlaps.list.24hrbft2.48hrbft2 = unique(overlaps.list.24hrbft2.48hrbft2) #here, i am using subject hits, which is the second peak file in findoverlaps. So, i am getting a number that represents the number of peaks in the second file that overlap with peaks in the first file. I always want number of peaks in the smaller file that overlaps with peaks in the larger file.


sample.peak.24hrbft2.length = sample.peak.24hrbft2[ ,1]
sample.peak.48hrbft2.length = sample.peak.48hrbft2[ ,1]

length(overlaps.list.24hrbft2.48hrbft2)
length(sample.peak.24hrbft2.length)
length(sample.peak.48hrbft2.length)


overlaps.24hrblank.48hrblank= as.data.frame(findOverlaps(sample.peak.48hrblank.gr, sample.peak.24hrblank.gr))
overlaps.list.24hrblank.48hrblank = overlaps.24hrblank.48hrblank [ ,"subjectHits"]
overlaps.list.24hrblank.48hrblank = unique(overlaps.list.24hrblank.48hrblank) #here, i am using subject hits, which is the second peak file in findoverlaps. So, i am getting a number that represents the number of peaks in the second file that overlap with peaks in the first file. I always want number of peaks in the smaller file that overlaps with peaks in the larger file.


sample.peak.24hrblank.length = sample.peak.24hrblank[ ,1]
sample.peak.48hrblank.length = sample.peak.48hrblank[ ,1]

length(overlaps.list.24hrblank.48hrblank)
length(sample.peak.24hrblank.length)
length(sample.peak.48hrblank.length)



overlaps.24hropened.48hropened= as.data.frame(findOverlaps(unique.peak.24hropened.gr, unique.peak.48hropened.gr))
overlaps.list.24hropened.48hropened = overlaps.24hropened.48hropened [ ,"subjectHits"]
overlaps.list.24hropened.48hropened = unique(overlaps.list.24hropened.48hropened) #here, i am using subject hits, which is the second peak file in findoverlaps. So, i am getting a number that represents the number of peaks in the second file that overlap with peaks in the first file. I always want number of peaks in the smaller file that overlaps with peaks in the larger file.


unique.peak.24hropened.length = unique.peak.24hropened[ ,1]
unique.peak.48hropened.length = unique.peak.48hropened[ ,1]

length(overlaps.list.24hropened.48hropened)
length(unique.peak.24hropened.length)
length(unique.peak.48hropened.length)



overlaps.24hrclosed.48hrclosed= as.data.frame(findOverlaps(unique.peak.48hrclosed.gr, unique.peak.24hrclosed.gr))
overlaps.list.24hrclosed.48hrclosed = overlaps.24hrclosed.48hrclosed [ ,"subjectHits"]
overlaps.list.24hrclosed.48hrclosed = unique(overlaps.list.24hrclosed.48hrclosed) #here, i am using subject hits, which is the second peak file in findoverlaps. So, i am getting a number that represents the number of peaks in the second file that overlap with peaks in the first file. I always want number of peaks in the smaller file that overlaps with peaks in the larger file.

unique.peak.24hrclosed.length = unique.peak.24hrclosed[ ,1]
unique.peak.48hrclosed.length = unique.peak.48hrclosed[ ,1]

length(unique.peak.24hrclosed.length)
length(unique.peak.48hrclosed.length)
length(overlaps.list.24hrclosed.48hrclosed)


                                        #make venn diagrams

pdf("~/Data/atac/final.analysis/rerun/venn2448closed.pdf")
draw.pairwise.venn(length(unique.peak.24hrclosed.length),length(unique.peak.48hrclosed.length),length(overlaps.list.24hrclosed.48hrclosed), c("24hrclosed", "48hrclosed"),fill= c("blue","red"), cat.pos= c(0,0), cat.dist=c(0.1,0.1), cex=1.5, print.mode=c('raw','percent'), inverted=TRUE)
dev.off()

pdf("~/Data/atac/final.analysis/rerun/venn2448opened.pdf")
draw.pairwise.venn(length(unique.peak.24hropened.length),length(unique.peak.48hropened.length),length(overlaps.list.24hropened.48hropened), c("24hropened", "48hropened"),fill= c("blue","red"), cat.pos= c(0,0), cat.dist=c(0.1,0.1), cex=1.5, print.mode=c('raw','percent'))
dev.off()

pdf("~/Data/atac/final.analysis/rerun/venn2448bft2.pdf")
draw.pairwise.venn(length(sample.peak.24hrbft2.length),length(sample.peak.48hrbft2.length),length(overlaps.list.24hrbft2.48hrbft2), c("24hrbft2", "48hrbft2"),fill= c("blue","red"), cat.pos= c(0,0), cat.dist=c(0.1,0.1), cex=1.5, print.mode=c('raw','percent'))
dev.off()

pdf("~/Data/atac/final.analysis/rerun/venn2448blank.pdf")
draw.pairwise.venn(length(sample.peak.24hrblank.length),length(sample.peak.48hrblank.length),length(overlaps.list.24hrblank.48hrblank), c("24hrblank", "48hrblank"),fill= c("blue","red"), cat.pos= c(0,0), cat.dist=c(0.1,0.1), cex=1.5, print.mode=c('raw','percent'), inverted=TRUE)
dev.off()
