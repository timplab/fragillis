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
##No 48 hr sig peaks


##Look for overlap
sig.24.op.gr=peak2gr(filter(res.24.merged.filt, stat>0))
sig.24.cl.gr=peak2gr(filter(res.24.merged.filt, stat<0))

sum(overlapsAny(sig.24.op.gr, ori.24open.gr))/length(sig.24.op.gr)
sum(overlapsAny(sig.24.cl.gr, ori.24close.gr))/length(sig.24.cl.gr)

sig.48.op.gr=peak2gr(filter(res.48.merged.filt, stat>0))
sig.48.cl.gr=peak2gr(filter(res.48.merged.filt, stat<0))

sum(overlapsAny(sig.48.op.gr, ori.48open.gr))/length(sig.48.op.gr)
sum(overlapsAny(sig.48.cl.gr, ori.48close.gr))/length(sig.48.cl.gr)



##look at promoter regions
if (FALSE) {
    promoters=promoters(TxDb.Hsapiens.UCSC.hg38.knownGene)
    
    sum(overlapsAny(sig.24.op.gr, promoters))
    sum(overlapsAny(sig.24.cl.gr, promoters))

    sig.op.prom=subsetByOverlaps(promoters, sig.24.op.gr)
    sig.cl.prom=subsetByOverlaps(promoters, sig.24.cl.gr)
    
    sig.op.genes=mapIds(org.Hs.eg.db, keys=as.character(sig.op.prom$tx_name), column="SYMBOL", keytype="UCSCKG")
    sig.cl.genes=mapIds(org.Hs.eg.db, keys=as.character(sig.cl.prom$tx_name), column="SYMBOL", keytype="UCSCKG")
    
}

##Check Cancer

if (FALSE) {

    ##First load in files
    cosmic.mut=read_tsv(file.path(outdir, "CosmicGenomeScreensMutantExport.tsv.gz")) %>%
        rename_all(funs(gsub(" ", ".", .))) %>%
        filter(Primary.site=="large_intestine") %>%

    cosmic.mut=cosmic.mut %>%
        separate(Mutation.genome.position, into=c("chrom", "chromStart", "chromEnd"), sep="[:-]", convert=T) %>%
        mutate(chrom=paste0("chr", chrom)) %>%
        mutate(chrom=ifelse(chrom=="chr23", "chrX", chrom)) %>%
        mutate(chrom=ifelse(chrom=="chr24", "chrY", chrom)) %>%
        mutate(strand=".", name=1:dim(cosmic.mut)[1])

    cosmic.mut=cosmic.mut %>%
        filter(!is.na(chromStart))
        

    

    crc = c("AKT1", "AKT2","AKT3","APC","APC2","APPL1","ARAF","Axin1","AXIN2","BAD","BAX","BCL2","BIRC5","BRAF",
            "casp3","CASP9","CCND1","CTNNB1","CYCS","DCC","FOS","GSK3B","JUN","KRAS","LEF1","MAP2K1","MAPK1",
            "MAPK10","MAPK3","MAPK8","MAPK9","MLH1","MSH2","MSH3","MSH6","MYC","PIK3CA","PIK3CB","PIK3CD",
            "PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","RAC1","RAC2","RAC3","RAF1","RALGDS","RHOA","SMAD2",
            "SMAD3","SMAD4","TCF7","TCF7L1","TCF7L2","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TP53")
    
    wnt = c("APC","ASCL1","ASCL2","AXIN1","AXIN2","CTNNB1","BCL9","BCL92","ABL1","CSNK1A1","CSNK1D","DIXDC1",
            "CREB","DVL1","DVL2","DVL3","DRAXIN","GSK3A","GSK3B","TGFB1I1","CTNNBIP1","PYGO1","PYGO2","HNF1B",
            "UBE2A","TCF7","TCF7L1","TLE2")


    bcat = c("MYC","MYCN","CCND1","HNF1A","LEF1","PPARD","JUN","FOSL1","PLAUR","MMP7","AXIN2","NRCAM","TCF4",
             "GAST","CD44","EFNB1","EFNB2","BMP4","CLDN1","BIRC5","VEGFA","FGF18","ATOH1","MET","EDN1","MYCBP",
             "L1CAM","ID2","JAG1","MSL1","TIAM1","NOS2","TERT","DKK1","FGF9","LBH","FGF20","LGR5","SOX9","SOX17",
             "RUNX2","GREM1","SALL4","TNFSF11","TNFRSF11B","CYR61","PTTG1","DLL1","FOXN1","MMP26","NANOG",
             "POU5F1","SNAI1","FN1","FZD7","ISL1","MMP2","MMP9","FST","WNT3A","TWIST1","TBX3","GBX2","CACNA1G",
             "CDC25","WISP1","WISP2","IGF2","EMP1","IGF1","VEGFC","PTGS2","IL6","PITX2","EGFR","CDH1","CDKN2A",
             "CTLA4","CXCL8","VCAN","TNFRSF19")
    ## INFO FROM https://web.stanford.edu/group/nusselab/cgi-bin/wnt/target_genes


    snvs.crc = filter(cosmic.mut, Gene.name %in% crc)  #snvs present in genes found most frequently mutated  in CRC
    snvs.wnt = filter(cosmic.mut, Gene.name %in% wnt)  # snvs present in genes in the wnt pathway
    snvs.bcat= filter(cosmic.mut, Gene.name %in% bcat) # snvs present in genes in the wnt pathway

    snvs.crc.gr=peak2gr(snvs.crc)
    snvs.wnt.gr=peak2gr(snvs.wnt)
    snvs.bcat.gr=peak2gr(snvs.bcat)

}    

    

##Check DMRs    
 
 
##Check CTCF




