###This code creates a heatmap  of gene expression from sleuth data using pheatmap

library(ggplot2)
library(reshape2)
library(pheatmap)
library(reshape)
library(tidyr)

#24hrs
#read in list of genes that are most differentially expressed at 24hrs
hr24genes= read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/24hr.expression.changes.genelevel.rerun.pval01.csv")
hr24genes= hr24genes[order(hr24genes$pval),]
hr24genes= hr24genes[, "target_id"] #gets top 20 differentially expressed genes and names

#read in gene expression data
sleuth.data = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv")
p.data = read.csv("~/Data/rnaseq_data/rnaseq.phenotypedata.csv", header=FALSE)
colnames(p.data) = c("seqnumber", "experiment", "treatment")
sleuth.data = merge(sleuth.data, p.data, by.x = "sample", by.y = "seqnumber")
sleuth.data = sleuth.data[ , c("sample","target_id","tpm","experiment","treatment")]
sleuth.data=sleuth.data[ sleuth.data$treatment %in% c("24hrbft2","24blank"),] #change to 48hrs if you want that data
sleuth.data$log_tpm= log2(sleuth.data$tpm +1)
sleuth.data.24hr= sleuth.data[ sleuth.data$target_id %in% hr24genes, ]


sleuth.data.24hr.tpm= sleuth.data.24hr[order(sleuth.data.24hr$sample),]
sleuth.data.24hr.tpm = subset(sleuth.data.24hr.tpm, select = -c(log_tpm))#removes  log tpm

sleuth.data.24hr.log=subset(sleuth.data.24hr, select = -c( tpm)) #removes tpm and leaves only logtpm
sleuth.data.24hr.log= sleuth.data.24hr.log[order(sleuth.data.24hr.log$sample),]

#set up data to be run with pheatmap function
setup.24= subset(sleuth.data.24hr.log, select = -c(treatment,sample)) #use the .log file or the .tpm file here. want to keep genes, logtpm, and sample identifier
setup.melt.24 = melt(setup.24,target_id=c("target_id","log_tpm")) #use log_tpm or tpm here
setup.melt.24 = cast(setup.melt.24,target_id~experiment, mean)

setup.melt.24.names= setup.melt.24
row.names(setup.melt.24.names)= setup.melt.24.names[ ,1]
setup.melt.24.names=setup.melt.24.names[ ,-1]#these three commands make my sample names my rownames



 #run pheatmap function and print pdf
pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/heatmap.24hrexp.24hrgenes.pdf")
pheatmap(setup.melt.24.names, fontsize=5) 
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/heatmap.24hrexp.24hrgenes.nocluster.pdf")
pheatmap(setup.melt.24.names, cluster_rows = F, cluster_cols=T) 

dev.off()





#48hrs
#read in list of genes that are most differentially expressed at 48hrs
hr48genes= read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/48hr.expression.changes.genelevel.rerun.pval01.csv")
hr48genes= hr48genes[order(hr48genes$pval),]
hr48genes= hr48genes[, "target_id"] #gets top 20 differentially expressed genes and names

#read in gene expression data
sleuth.data = read.csv("~/Data/rnaseq_data/gene.level.analysis/tpmcounts.sleuth.genelevel.allsamples.rerun.csv")
p.data = read.csv("~/Data/rnaseq_data/rnaseq.phenotypedata.csv", header=FALSE)
colnames(p.data) = c("seqnumber", "experiment", "treatment")
sleuth.data = merge(sleuth.data, p.data, by.x = "sample", by.y = "seqnumber")
sleuth.data = sleuth.data[ , c("sample","target_id","tpm","experiment","treatment")]
sleuth.data=sleuth.data[ sleuth.data$treatment %in% c("48hrbft2","48blank"),] #change to 48hrs if you want that data
sleuth.data$log_tpm= log2(sleuth.data$tpm +1)
sleuth.data.48hr= sleuth.data[ sleuth.data$target_id %in% hr48genes, ]


sleuth.data.48hr.tpm= sleuth.data.48hr[order(sleuth.data.48hr$sample),]
sleuth.data.48hr.tpm = subset(sleuth.data.48hr.tpm, select = -c(log_tpm))#removes  log tpm

sleuth.data.48hr.log=subset(sleuth.data.48hr, select = -c( tpm)) #removes tpm and leaves only logtpm
sleuth.data.48hr.log= sleuth.data.48hr.log[order(sleuth.data.48hr.log$sample),]




#set up data to be run with pheatmap function
setup.48= subset(sleuth.data.48hr.log, select = -c(treatment,sample)) #use the .log file or the .tpm file here. want to keep genes, logtpm, and sample identifier
setup.melt.48 = melt(setup.48,target_id=c("target_id","log_tpm")) #use log_tpm or tpm here
setup.melt.48 = cast(setup.melt.48,target_id~experiment, mean)

setup.melt.48.names= setup.melt.48
row.names(setup.melt.48.names)= setup.melt.48.names[ ,1]
setup.melt.48.names=setup.melt.48.names[ ,-1]#these three commands make my sample names my rownames



#run pheatmap function and print pdf

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/heatmap.48hrexp.48hrgenes.pdf")
pheatmap(setup.melt.48.names)
dev.off()

pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/heatmap.48hrexp.48hrgenes.nocluster.pdf")
pheatmap(setup.melt.48.names, cluster_rows=F, cluster_cols=T) 
dev.off()

