
  #Use fold change change data to look for trends in gene expression changes and make line plots
 
install.packages("ggplot2")
install.packages("reshape2")
install.packages("pheatmap")
install.packages("tidyr")
install.packages("reshape")
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyr)
library(reshape)


#Here, I am reading in data that describes the samples associated with each file name
phenodata2 = read.csv("~/Data/rnaseq_data/rnaseq.phenotypedata.csv", header = FALSE)

 #here, I am creating lists with specific genes I want to look at


wnt = c("APC","ASCL1","ASCL2","AXIN1","AXIN2","CTNNB1","BCL9","BCL92","ABL1","CSNK1A1","CSNK1D","DIXDC1","CREB","DVL1","DVL2","DVL3","DRAXIN","GSK3A","GSK3B","TGFB1I1","CTNNBIP1", "PYGO1","PYGO2","HNF1B","UBE2A","TCF7","TCF7L1","TLE2")


crc = c("AKT1", "AKT2","AKT3","APC","APC2","APPL1","ARAF","Axin1","AXIN2","BAD","BAX","BCL2","BIRC5","BRAF","casp3","CASP9","CCND1","CTNNB1","CYCS","DCC","FOS","GSK3B","JUN","KRAS","LEF1","MAP2K1","MAPK1","MAPK10","MAPK3","MAPK8","MAPK9","MLH1","MSH2","MSH3","MSH6","MYC","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","RAC1","RAC2","RAC3","RAF1","RALGDS","RHOA","SMAD2","SMAD3","SMAD4","TCF7","TCF7L1","TCF7L2","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TP53")

bcat = c("MYC","MYCN","CCND1","HNF1A","LEF1","PPARD","JUN","FOSL1","PLAUR","MMP7","AXIN2","NRCAM","TCF4","GAST","CD44","EFNB1","EFNB2","BMP4","CLDN1","BIRC5","VEGFA","FGF18","ATOH1","MET","EDN1","MYCBP","L1CAM","ID2","JAG1","MSL1","TIAM1","NOS2","TERT","DKK1","FGF9","LBH","FGF20","LGR5","SOX9","SOX17","RUNX2","GREM1","SALL4","TNFSF11","TNFRSF11B","CYR61","PTTG1","DLL1","FOXN1","MMP26","NANOG","POU5F1","SNAI1","FN1","FZD7","ISL1","MMP2","MMP9","FST","WNT3A","TWIST1","TBX3","GBX2","CACNA1G","CDC25","WISP1","WISP2","IGF2","EMP1","IGF1","VEGFC","PTGS2","IL6","PITX2","EGFR","CDH1","CDKN2A","CTLA4","CXCL8","VCAN","TNFRSF19") # INFO FROM https://web.stanford.edu/group/nusselab/cgi-bin/wnt/target_genes


                                        #read in 24hrbft2 expression data
expres.24 = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv", row.names = 1)
expres.24 = expres.24[expres.24$condition %in% c("24hrbft2"),]
expres.24 = expres.24[order(expres.24$target_id),]
expres.24$logtpm= log2(expres.24$tpm + 1)
#expres.24 = expres.24[expres.24$logtpm>1,]
expres.24 =  merge(expres.24, phenodata2, by.x= "sample", by.y= "V1")
expres.24 = expres.24[order(expres.24$target_id),]

                                        #read in 24hrblank data
expres.24b = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv", row.names = 1)
expres.24b = expres.24b[expres.24b$condition %in% c("24blank"),]
expres.24b = expres.24b[order(expres.24b$target_id),]
expres.24b$logtpm= log2(expres.24b$tpm + 1)
#expres.24b = expres.24b[expres.24b$logtpm>1,]
expres.24b =  merge(expres.24b, phenodata2, by.x= "sample", by.y= "V1")
expres.24b = expres.24b[order(expres.24b$target_id),]

                                        #split up data for each group
expres.24.220 = expres.24[expres.24$V2=="24hrbft2.220",]
expres.24.031 = expres.24[expres.24$V2=="24hrbft2.031",]
expres.24.trial1 = expres.24[expres.24$V2=="24hrbft2.trial1",]

expres.24b.220 = expres.24b[expres.24b$V2=="24blank.220",]
expres.24b.031 = expres.24b[expres.24b$V2=="24blank.031",]
expres.24b.trial1 = expres.24b[expres.24b$V2=="24blank.trial1",]


#calculate the logfc in expression between blank and bft2 samples then combine all the samples, then get rid of low expressing genes
expres.24.220$logfc = log2((1+expres.24.220$tpm)/(1+expres.24b.220$tpm))
expres.24.031$logfc = log2((1+expres.24.031$tpm)/(1+expres.24b.031$tpm))
expres.24.trial1$logfc = log2((1+expres.24.trial1$tpm)/(1+expres.24b.trial1$tpm))

expres.24.all=rbind(expres.24.220, expres.24.trial1, expres.24.031)
#expres.24.all = expres.24.all[expres.24.all$logtpm>1,]
 

                                        #read in 48hrbft2 expression data
expres.48 = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv", row.names = 1)
expres.48 = expres.48[expres.48$condition %in% c("48hrbft2"),]
expres.48 = expres.48[order(expres.48$target_id),]
expres.48$logtpm= log2(expres.48$tpm + 1)
#expres.48 = expres.48[expres.48$logtpm>1,]
expres.48 =  merge(expres.48, phenodata2, by.x= "sample", by.y= "V1")
expres.48 = expres.48[order(expres.48$target_id),]

                                        #read in 48hrblank data
expres.48b = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv", row.names = 1)
expres.48b = expres.48b[expres.48b$condition %in% c("48blank"),]
expres.48b = expres.48b[order(expres.48b$target_id),]
expres.48b$logtpm= log2(expres.48b$tpm + 1)
#expres.48b = expres.48b[expres.48b$logtpm>1,]
expres.48b =  merge(expres.48b, phenodata2, by.x= "sample", by.y= "V1")
expres.48b = expres.48b[order(expres.48b$target_id),]

                                        #separate out replicate groups
expres.48.220 = expres.48[expres.48$V2=="48hrbft2.220",]
expres.48.031 = expres.48[expres.48$V2=="48hrbft2.031",]
expres.48.trial1 = expres.48[expres.48$V2=="48hrbft2.trial1",]

expres.48b.220 = expres.48b[expres.48b$V2=="48blank.220",]
expres.48b.031 = expres.48b[expres.48b$V2=="48blank.031",]
expres.48b.trial1 = expres.48b[expres.48b$V2=="48blank.trial1",]


#calculate logfc for the differnce in expression between bft2 and blank, then combine samples, then get rid of low expressing genes
expres.48.220$logfc = log2((1+expres.48.220$tpm)/(1+expres.48b.220$tpm))
expres.48.031$logfc = log2((1+expres.48.031$tpm)/(1+expres.48b.031$tpm))
expres.48.trial1$logfc = log2((1+expres.48.trial1$tpm)/(1+expres.48b.trial1$tpm))

expres.48.all = rbind(expres.48.220, expres.48.trial1, expres.48.031)#do not remove 031 time point because it looks like an outlier
#expres.48.all = expres.48.all[expres.48.all$logtpm>1,]




####################



#set up data to be run with pheatmap function

expres.2448.wnt = rbind(expres.24.all, expres.48.all)


                                        #get average logfc
setup.expres.wnt =  expres.2448.wnt[ ,c("target_id","V2","logfc")]
setup.expres.wnt.melt = melt(setup.expres.wnt, target_id=c("target_id","logfc"))
setup.expres.wnt.melt.cast = cast(setup.expres.wnt.melt, target_id~V2, mean)

setup.expres.wnt.melt.names = setup.expres.wnt.melt.cast
row.names(setup.expres.wnt.melt.names) = setup.expres.wnt.melt.names[ ,1]
setup.expres.wnt.melt.names = setup.expres.wnt.melt.names[ ,-1] #these 3 commands make my sample names my row names
colnames(setup.expres.wnt.melt.names) = c("hr24bft2.031","hr24bft2.220","hr24bft2.trial1","hr48bft2.031","hr48bft2.220","hr48bft2.trial1")



expres.wnt.averages = setup.expres.wnt.melt.names
expres.wnt.averages$average24hrbft2 = (expres.wnt.averages$hr24bft2.031 + expres.wnt.averages$hr24bft2.220 + expres.wnt.averages$hr24bft2.trial1)/3
expres.wnt.averages$average48hrbft2 = (expres.wnt.averages$hr48bft2.031 + expres.wnt.averages$hr48bft2.220 + expres.wnt.averages$hr48bft2.trial1)/3
expres.wnt.averages = expres.wnt.averages [ , c("average24hrbft2","average48hrbft2")]
expres.wnt.averages.df = data.frame(expres.wnt.averages)
expres.wnt.averages.df$target_id= rownames(expres.wnt.averages.df)
expres.wnt.averages.melt = melt(expres.wnt.averages.df)
colnames(expres.wnt.averages.melt) = c("target_id","variable","average_logfc")

#get average tpm values
setup.expres.wnt.tpm =  expres.2448.wnt[ ,c("target_id","V2","tpm")]
setup.expres.wnt.tpm.melt = melt(setup.expres.wnt.tpm, target_id=c("target_id","tpm"))
setup.expres.wnt.tpm.melt.cast = cast(setup.expres.wnt.tpm.melt, target_id~V2, mean)

setup.expres.wnt.tpm.melt.names = setup.expres.wnt.tpm.melt.cast
row.names(setup.expres.wnt.tpm.melt.names) = setup.expres.wnt.tpm.melt.names[ ,1]
setup.expres.wnt.tpm.melt.names = setup.expres.wnt.tpm.melt.names[ ,-1] #these 3 commands make my sample names my row names
colnames(setup.expres.wnt.tpm.melt.names) = c("hr24bft2.031","hr24bft2.220","hr24bft2.trial1","hr48bft2.031","hr48bft2.220","hr48bft2.trial1")


expres.wnt.tpm.averages = setup.expres.wnt.tpm.melt.names
expres.wnt.tpm.averages$average24hrbft2 = (expres.wnt.tpm.averages$hr24bft2.031 + expres.wnt.tpm.averages$hr24bft2.220 + expres.wnt.tpm.averages$hr24bft2.trial1)/3
expres.wnt.tpm.averages$average48hrbft2 = (expres.wnt.tpm.averages$hr48bft2.031 + expres.wnt.tpm.averages$hr48bft2.220 + expres.wnt.tpm.averages$hr48bft2.trial1)/3
expres.wnt.tpm.averages = expres.wnt.tpm.averages [ , c("average24hrbft2","average48hrbft2")]
expres.wnt.tpm.averages.df = data.frame(expres.wnt.tpm.averages)
expres.wnt.tpm.averages.df$target_id= rownames(expres.wnt.tpm.averages.df)
expres.wnt.tpm.averages.melt = melt(expres.wnt.tpm.averages.df)
colnames(expres.wnt.tpm.averages.melt) = c("target_id","variable","average_tpm")

                                        #here, I am separating out my 24hr and 48hr data, then merging my logfc data with my average tpm data
expres.wnt.averages.melt.24 = expres.wnt.averages.melt[expres.wnt.averages.melt$variable=="average24hrbft2", ]
expres.wnt.averages.melt.48 = expres.wnt.averages.melt[expres.wnt.averages.melt$variable=="average48hrbft2", ]

expres.wnt.tpm.averages.melt.24 = expres.wnt.tpm.averages.melt[expres.wnt.tpm.averages.melt$variable=="average24hrbft2", ]
expres.wnt.tpm.averages.melt.48 = expres.wnt.tpm.averages.melt[expres.wnt.tpm.averages.melt$variable=="average48hrbft2", ]

expres.wnt.tpm.logfc.24 =  merge(expres.wnt.tpm.averages.melt.24, expres.wnt.averages.melt.24, by.x= "target_id", by.y= "target_id")
expres.wnt.tpm.logfc.48 =  merge(expres.wnt.tpm.averages.melt.48, expres.wnt.averages.melt.48, by.x= "target_id", by.y= "target_id")

                                        #here, I am making a new logtpm column
expres.wnt.tpm.logfc.24$logtpm = log2(expres.wnt.tpm.logfc.24$average_tpm +1)
expres.wnt.tpm.logfc.48$logtpm = log2(expres.wnt.tpm.logfc.48$average_tpm +1)







                                        #Add chromatin state data to  all of my samples

                                        #first, I will read in data with info on uniquepeak locations
upeaks.24hrbft2= read.csv("~/Data/atac/final.analysis/chippeakanno_peaks/unique24hrbft2peaks.promoters.1000.1.chippeakanno.csv") #use promoters only or promoters.andgenebodies file

upeaks.24hrblank= read.csv("~/Data/atac/final.analysis/chippeakanno_peaks/unique24hrblankpeaks.promoters.1000.1.chippeakanno.csv") #use promoters only or promoters.andgenebodies file

upeaks.48hrbft2= read.csv("~/Data/atac/final.analysis/chippeakanno_peaks/unique48hrbft2peaks.promoters.1000.1.chippeakanno.csv") #use promoters only or promoters.andgenebodies file

upeaks.48hrblank= read.csv("~/Data/atac/final.analysis/chippeakanno_peaks/unique48hrblankpeaks.promoters.1000.1.chippeakanno.csv") #use promoters only or promoters.andgenebodies file




                                        #genes in wnt pathway
expres.wnt.tpm.logfc.24$bft48hr.chromatin="nochange"
expres.wnt.tpm.logfc.24$bft48hr.chromatin[expres.wnt.tpm.logfc.24$target_id %in% upeaks.48hrbft2$symbol]="opened"
expres.wnt.tpm.logfc.24$bft48hr.chromatin[expres.wnt.tpm.logfc.24$target_id %in% upeaks.48hrblank$symbol]="closed"

expres.wnt.tpm.logfc.24$bft24hr.chromatin="nochange"
expres.wnt.tpm.logfc.24$bft24hr.chromatin[expres.wnt.tpm.logfc.24$target_id %in% upeaks.24hrbft2$symbol]="opened"
expres.wnt.tpm.logfc.24$bft24hr.chromatin[expres.wnt.tpm.logfc.24$target_id %in% upeaks.24hrblank$symbol]="closed"




expres.wnt.tpm.logfc.48$bft48hr.chromatin="nochange"
expres.wnt.tpm.logfc.48$bft48hr.chromatin[expres.wnt.tpm.logfc.48$target_id %in% upeaks.48hrbft2$symbol]="opened"
expres.wnt.tpm.logfc.48$bft48hr.chromatin[expres.wnt.tpm.logfc.48$target_id %in% upeaks.48hrblank$symbol]="closed"

expres.wnt.tpm.logfc.48$bft24hr.chromatin="nochange"
expres.wnt.tpm.logfc.48$bft24hr.chromatin[expres.wnt.tpm.logfc.48$target_id %in% upeaks.24hrbft2$symbol]="opened"
expres.wnt.tpm.logfc.48$bft24hr.chromatin[expres.wnt.tpm.logfc.48$target_id %in% upeaks.24hrblank$symbol]="closed"


table(expres.wnt.tpm.logfc.24$bft24hr.chromatin)
table(expres.wnt.tpm.logfc.24$bft48hr.chromatin)
table(expres.wnt.tpm.logfc.24$bftboth.chromatin)
table(expres.wnt.tpm.logfc.48$bft24hr.chromatin)
table(expres.wnt.tpm.logfc.48$bft48hr.chromatin)
table(expres.wnt.tpm.logfc.48$bftboth.chromatin)

expres.wnt.tpm.logfc.24.nono = expres.wnt.tpm.logfc.24[ expres.wnt.tpm.logfc.24$bft24hr.chromatin !="nochange",] #writes only genes that are opened or closed
expres.wnt.tpm.logfc.48.nono = expres.wnt.tpm.logfc.48[ expres.wnt.tpm.logfc.48$bft48hr.chromatin !="nochange",] #writes only genes that are opened or closed


   #perform stats
wilcox.test (expres.wnt.tpm.logfc.24[ expres.wnt.tpm.logfc.24$bft24hr.chromatin=="opened", "average_logfc"], expres.wnt.tpm.logfc.24[ expres.wnt.tpm.logfc.24$bft24hr.chromatin=="closed","average_logfc" ]) #testing whether baseline expression is different for open and closed genes

mean(expres.wnt.tpm.logfc.24[expres.wnt.tpm.logfc.24$bft24hr.chromatin=="opened","average_logfc"])
mean(expres.wnt.tpm.logfc.24[expres.wnt.tpm.logfc.24$bft24hr.chromatin=="closed","average_logfc"])

median(expres.wnt.tpm.logfc.24[expres.wnt.tpm.logfc.24$bft24hr.chromatin=="opened","average_logfc"])
median(expres.wnt.tpm.logfc.24[expres.wnt.tpm.logfc.24$bft24hr.chromatin=="closed","average_logfc"])




wilcox.test (expres.wnt.tpm.logfc.48[ expres.wnt.tpm.logfc.48$bft48hr.chromatin=="opened", "average_logfc"], expres.wnt.tpm.logfc.48[ expres.wnt.tpm.logfc.48$bft48hr.chromatin=="closed","average_logfc" ]) 

mean(expres.wnt.tpm.logfc.48[expres.wnt.tpm.logfc.48$bft48hr.chromatin=="opened","average_logfc"])
mean(expres.wnt.tpm.logfc.48[expres.wnt.tpm.logfc.48$bft48hr.chromatin=="closed","average_logfc"])

median(expres.wnt.tpm.logfc.48[expres.wnt.tpm.logfc.48$bft48hr.chromatin=="opened","average_logfc"])
median(expres.wnt.tpm.logfc.48[expres.wnt.tpm.logfc.48$bft48hr.chromatin=="closed","average_logfc"])



#get rid of low expression genes
expres.wnt.tpm.logfc.24.nolow = expres.wnt.tpm.logfc.24[expres.wnt.tpm.logfc.24$logtpm>1,]
expres.wnt.tpm.logfc.48.nolow = expres.wnt.tpm.logfc.48[expres.wnt.tpm.logfc.48$logtpm>1,]

expres.wnt.tpm.logfc.24.nono.nolow = expres.wnt.tpm.logfc.24.nolow[ expres.wnt.tpm.logfc.24.nolow$bft24hr.chromatin !="nochange",] #writes only genes that are opened or closed
expres.wnt.tpm.logfc.48.nono.nolow = expres.wnt.tpm.logfc.48.nolow[ expres.wnt.tpm.logfc.48.nolow$bft48hr.chromatin !="nochange",] #writes only genes that are opened or closed


#stats with low expression genes removed
wilcox.test (expres.wnt.tpm.logfc.24.nolow[ expres.wnt.tpm.logfc.24.nolow$bft24hr.chromatin=="opened", "average_logfc"], expres.wnt.tpm.logfc.24.nolow[ expres.wnt.tpm.logfc.24.nolow$bft24hr.chromatin=="closed","average_logfc" ]) 

wilcox.test (expres.wnt.tpm.logfc.48.nolow[ expres.wnt.tpm.logfc.48.nolow$bft48hr.chromatin=="opened", "average_logfc"], expres.wnt.tpm.logfc.48.nolow[ expres.wnt.tpm.logfc.48.nolow$bft48hr.chromatin=="closed","average_logfc" ]) 


#plot
pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.24hrs.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.24, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft24hr.chromatin), alpha=0.75) + geom_abline(intercept=0.5, slope=0, linetype=2)+ geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.015 \nwilcoxon rank sum test") + ylim(-1.5,1.5) + xlim(0,15)
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.24hrs.nonochange.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.24.nono, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft24hr.chromatin), alpha=0.75) + geom_abline(intercept=0.5, slope=0, linetype=2)+ geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.015 \nwilcoxon rank sum test") +  ylim(-1.5,1.5) + xlim(0,15) 
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.24hrs.nonochange.nolow.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.24.nono.nolow, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft24hr.chromatin), alpha=0.75) + geom_abline(intercept=0.5, slope=0, linetype=2)+ geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.009 \nwilcoxon rank sum test") +  ylim(-1.5,1.5) + xlim(0,15) 
dev.off()


pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.48hrs.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.48, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft48hr.chromatin), alpha=0.75) + geom_abline(intercept=0.5, slope=0, linetype=2) + geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.456 \nwilcoxon rank sum test") +  ylim(-1.5,1.5) + xlim(0,15)
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.48hrs.nonochange.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.48.nono, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft48hr.chromatin),alpha=0.75) + geom_abline(intercept=0.5, slope=0, linetype=2)+ geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.456 \nwilcoxon rank sum test") +  ylim(-1.5,1.5) + xlim(0,15)
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/maplot.48hrs.nonochange.nolow.rerun.pdf")
ggplot( expres.wnt.tpm.logfc.48.nono.nolow, aes(x=logtpm, y=average_logfc)) +  geom_point(aes(colour=bft48hr.chromatin),alpha=0.75) + geom_abline(intercept=0.5, slope=0)+ geom_abline(slope=0, intercept = -0.5, linetype=2) +labs(x= "mean (logTPM +1) ", y="logFC expression") + annotate("text", x = 10,y=1.3, label = "pvalue [opened ~ closed] =0.383 \nwilcoxon rank sum test") +  ylim(-1.5,1.5) + xlim(0,15)
dev.off()


##################

