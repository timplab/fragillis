library(tidyverse)
library(googlesheets)
library(GenomicRanges)
#library(GenomicAlignements)
library(Rsubread)
library(DESeq2)

#Code dir is ~/Code/

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_title("Fragillis Data")
dataloc=gs_read(fullsheet, ws="180608_dat")

workdir="/atium/Data/NGS/Aligned/180806_fragillis"
outdir="~/Dropbox/timplab_data/fragillis/180806_revision"

if (!dir.exists(workdir)) {
    dir.create(workdir)
}

if (!dir.exists(outdir)) {
    dir.create(outdir)
}


setwd(workdir)


peak2gr <- function(peak) {    
    mygr=GRanges(seqnames=peak$chrom, ranges=IRanges(start=peak$chromStart, end=peak$chromEnd), name=peak$name)
}




##Load peak files

narrowpeak=c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

merged=read_tsv(file.path(workdir, "sharpmerged_peaks.narrowPeak"), col_names=narrowpeak)
merged.gr=peak2gr(merged)

ori.24open=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrbft2.narrowPeak"), col_names=narrowpeak)
ori.24open.gr=peak2gr(ori.24open)

ori.24close=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrblank.narrowPeak"), col_names=narrowpeak)
ori.24close.gr=peak2gr(ori.24close)

ori.48open=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.48hrbft2.narrowPeak"), col_names=narrowpeak)
ori.48open.gr=peak2gr(ori.48open)

ori.48close=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrblank.narrowPeak"), col_names=narrowpeak)
ori.48close.gr=peak2gr(ori.48close)


sum(overlapsAny(ori.24open.gr, merged.gr))/length(ori.24open.gr)
sum(overlapsAny(ori.24close.gr, merged.gr))/length(ori.24close.gr)
sum(overlapsAny(ori.48open.gr, merged.gr))/length(ori.48open.gr)
sum(overlapsAny(ori.48close.gr, merged.gr))/length(ori.48close.gr)


##Jawara filters here based on qval, I'm going to skip it for now


##Use featureCounts to count reads in features
##first convert/extract SAF from merged.gr
##Filter merged peaks to keep only autosomes
merged.saf=merged %>%
    filter(!(chrom %in% c('chrX', 'chrY'))) %>%
    select(chrom, chromStart, chromEnd, name, strand) %>%
    rename_(GeneID='name', Chr='chrom', Start='chromStart', End='chromEnd', Strand='strand')

merged.saf=as.data.frame(merged.saf[,c(4, 1, 2, 3, 5)])

feature = featureCounts(annot.ext=merged.saf,
                        nthreads=10,                            
                        files=file.path(workdir, paste0(dataloc$label, "nodup.bam")))

counts=feature$counts

colnames(counts)=dataloc$label

pdata=data.frame(time=dataloc$time, bft=dataloc$bft=="y", row.names=dataloc$label)


##Just do 24 hrs
timedo=24
dds=DESeqDataSetFromMatrix(countData = counts[,pdata$time==timedo], colData=pdata[pdata$time==timedo,], design= ~ bft)

##pre-filter less than 10 reads
dds=dds[rowSums(counts(dds)) >= 10,]
##Tell it control is untreated
dds$bft = factor(dds$bft, levels=c(FALSE, TRUE))
##run
dds=DESeq(dds)

res=as.tibble(as.data.frame(results(dds)), rownames="name") %>%
    arrange(pvalue)

pdf(file.path(outdir, "24hrMA.pdf"))
plotMA(results(dds))
dev.off()


res.filt=res %>%
    filter(pvalue<.01) %>%
    left_join(merged)

sig..gr=peak2gr(res.filt)

sum(overlapsAny(sig.24.gr, ori.24open.gr))/length(sig.24.gr)
sum(overlapsAny(sig.24.gr, ori.24close.gr))/length(sig.24.gr)

