library(googlesheets)
library(GenomicRanges)
#library(GenomicAlignements)
library(Rsubread)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BiocParallel)
library(tidyverse)

register(MulticoreParam(8))

##Code dir is ~/Code/

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

gr2peak <- function(mygr) {
    peak=tibble(chrom=as.character(seqnames(mygr)), chromStart=start(mygr), chromEnd=end(mygr), strand=as.character(strand(mygr)))
    peak$name=as.character(1:length(mygr))
    return(peak)
}

call.sig <- function(peaks, dataloc, timedo=24) {

    just.these=dataloc[dataloc$time==timedo,]
    
    peaks.saf=peaks %>%
        filter(!(chrom %in% c('chrX', 'chrY'))) %>%
        select(chrom, chromStart, chromEnd, name, strand) %>%
        rename_(GeneID='name', Chr='chrom', Start='chromStart', End='chromEnd', Strand='strand')
    
    peaks.saf=as.data.frame(peaks.saf[,c(4, 1, 2, 3, 5)])
    
    feature = featureCounts(annot.ext=peaks.saf,
                            nthreads=10,                            
                            files=file.path(workdir, paste0(just.these$label, "nodup.bam")))
    
    counts=feature$counts
    
    colnames(counts)=just.these$label
    
    pdata=data.frame(time=just.these$time, bft=just.these$bft=="y", row.names=just.these$label)
    
    
    ##Just do 24 hrs
    ##timedo=24
    dds=DESeqDataSetFromMatrix(countData = counts, colData=pdata, design= ~ bft)
    
    ##pre-filter less than 10 reads
    dds=dds[rowSums(counts(dds)) >= 10,]
    ##Tell it control is untreated
    dds$bft = factor(dds$bft, levels=c(FALSE, TRUE))
    ##run
    dds=DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(8))
    
    res=as.tibble(as.data.frame(results(dds)), rownames="name") %>%
        arrange(pvalue)
    
    res.filt=res %>%
        filter(pvalue<.01) %>%
        left_join(peaks)

    return(res.filt)
}    



##Load peak files

narrowpeak=c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

merged=read_tsv(file.path(workdir, "sharpmerged_peaks.narrowPeak"), col_names=narrowpeak)
merged.gr=peak2gr(merged)

ori.24open=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrbft2.narrowPeak"), col_names=narrowpeak)
ori.24open.gr=peak2gr(ori.24open)

ori.24close=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrblank.narrowPeak"), col_names=narrowpeak)
ori.24close.gr=peak2gr(ori.24close)

ori.24full.gr=GenomicRanges::reduce(c(ori.24open.gr, ori.24close.gr))
ori.24full=gr2peak(ori.24full.gr)


ori.48open=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.48hrbft2.narrowPeak"), col_names=narrowpeak)
ori.48open.gr=peak2gr(ori.48open)

ori.48close=read_tsv(file.path(workdir, "old", "uniquepeaks.optimalset.24hrblank.narrowPeak"), col_names=narrowpeak)
ori.48close.gr=peak2gr(ori.48close)

ori.48full.gr=GenomicRanges::reduce(c(ori.48open.gr, ori.48close.gr))
ori.48full=gr2peak(ori.48full.gr)

sum(overlapsAny(ori.24open.gr, merged.gr))/length(ori.24open.gr)
sum(overlapsAny(ori.24close.gr, merged.gr))/length(ori.24close.gr)
sum(overlapsAny(ori.48open.gr, merged.gr))/length(ori.48open.gr)
sum(overlapsAny(ori.48close.gr, merged.gr))/length(ori.48close.gr)





##Jawara filters here based on qval, I'm going to skip it for now

res.24.merged.filt=call.sig(merged, dataloc, timedo=24)
res.48.merged.filt=call.sig(merged, dataloc, timedo=48)

res.24.old.filt=call.sig(ori.24full, dataloc, timedo=24)
res.48.old.filt=call.sig(ori.48full, dataloc, timedo=48)




write_csv(res.24.merged.filt, file.path(outdir, "24hrsig.merged.csv"))
write_csv(res.48.merged.filt, file.path(outdir, "48hrsig.merged.csv"))

write_csv(res.24.old.filt, file.path(outdir, "24hrsig.oldpeak.csv"))

res.open=filter(res.filt, stat>0)
res.close=filter(res.filt, stat<0)

sig.24.gr=peak2gr(res.filt)
sig.24.op.gr=peak2gr(res.open)
sig.24.cl.gr=peak2gr(res.close)

sum(overlapsAny(sig.24.op.gr, ori.24open.gr))/length(sig.24.op.gr)
sum(overlapsAny(sig.24.cl.gr, ori.24close.gr))/length(sig.24.cl.gr)

##make MA plot

##look at promoter regions
if (FALSE) {
    promoters=promoters(TxDb.Hsapiens.UCSC.hg38.knownGene)
    
    sum(overlapsAny(sig.24.op.gr, promoters))
    sum(overlapsAny(sig.24.cl.gr, promoters))

    sig.op.prom=subsetByOverlaps(promoters, sig.24.op.gr)
    sig.cl.prom=subsetByOverlaps(promoters, sig.24.cl.gr)
    

    sig.op.genes=mapIds(org.Hs.eg.db, keys=as.character(sig.op.prom$tx_id), column="SYMBOL", keytype="ENTREZID")
    sig.cl.genes=mapIds(org.Hs.eg.db, keys=as.character(sig.cl.prom$tx_id), column="SYMBOL", keytype="ENTREZID")
    

}




