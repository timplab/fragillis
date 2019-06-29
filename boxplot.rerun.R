library(ggplot2)



###code used to calculate baseline differences in expresion by peak status

 #first, read in location of peaks and data describing file names
regpeaks = read.csv("~/Data/atac/final.analysis/rerun/chippeakanno_peaks/24hrblankpeaks.promoters.1000.1.chippeakanno.csv")
regpeaks = regpeaks[ ,"symbol"]
phenodata2 = read.csv("~/Data/rnaseq_data/rnaseq.phenotypedata.csv", header = FALSE)

  # then, read in expresion data and edit it
expres.baseline = read.csv("~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv", row.names = 1)
expres.baseline = expres.baseline[expres.baseline$condition=="24blank",]
expres.baseline$logtpm= log2(expres.baseline$tpm + 1)
expres.baseline= expres.baseline[expres.baseline$logtpm>1,]
expres.baseline =  merge(expres.baseline, phenodata2, by.x= "sample", by.y= "V1")


  #add chromatin structure data to gene expression data
expres.baseline$bft24hr.chromatin="nopeak"
expres.baseline$bft24hr.chromatin[expres.baseline$target_id %in% regpeaks]="peak"


table(expres.baseline$bft24hr.chromatin)



  #perform stats
wilcox.test (expres.baseline[ expres.baseline$bft24hr.chromatin=="nopeak", "logtpm"], expres.baseline[ expres.baseline$bft24hr.chromatin=="peak","logtpm" ]) #testing whether baseline expression is different for open and closed genes

mean(expres.baseline[expres.baseline$bft24hr.chromatin=="nopeak","logtpm"])
mean(expres.baseline[expres.baseline$bft24hr.chromatin=="peak","logtpm"])

median(expres.baseline[expres.baseline$bft24hr.chromatin=="nopeak","logtpm"])
median(expres.baseline[expres.baseline$bft24hr.chromatin=="peak","logtpm"])



                                        #plot boxplot of data
pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/baseline_atac_rna_comparison__pandgb_boxplot.pdf")
ggplot(expres.baseline, aes( x =bft24hr.chromatin, y=logtpm)) + geom_boxplot() +  labs(title="correlation between expression and  \npromoter/genebody chromatin structure at baseline", y="logtpm expression", x = "chromatin state at baseline") + annotate("text", x = 2, y=10, label = "pvalue [opened ~ closed] = 2e^-16 \nwilcoxon ranksum test")

dev.off()



pdf("~/Data/rnaseq_data/gene.level.analysis/rerun/baseline_atac_rna_comparison_histo.pdf")
ggplot(expres.baseline, aes(x=logtpm, fill=bft24hr.chromatin)) + geom_density(alpha=0.3) + annotate("text", x = 10, y=0.20, label = "pvalue [opened ~ closed] = 2e^-16 \nwilcoxon ranksum test")
dev.off()


