library(tidyverse)

##Ok - first pull in csv files giving locations of rna data
##


##First one doesn't have info, but I pulled an excel from genesifter, made a csv
meta=read_csv("/atium/Data/NGS/Aligned/170613_fragillis/meta/128239.csv", skip=2) %>%
    rename(idx=`Index Sequence`) %>%
    rename(samp=Template) %>%
    select(samp, idx) %>%
    mutate(run=1, paired=T, rawdir="/mithril/Data/NGS/Raw/160708_fragillis/wtimp1_128239/FASTQ/")

##Ok - Jawara renamed files, so manual editing here of idx
##First, time0 is "Time 0" - remove whitespace
meta$samp=gsub(" ", "", meta$samp)

##Now make idex from samp name instead of actual barcode idx
meta$idx=meta$samp

##CHECK IF THIS WORKS


##continue to build full filename list
##Ok - genesifter table for 2 doesn't have idx codes, using this file instead

temp=read_csv("/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135818/wtimp1_135818.csv") %>%
    rename(idx=Index) %>%
    rename(samp=SM_Tag) %>%
    filter(Lane==1) %>%
    select(samp, idx) %>%
    mutate(run=2, paired=T, rawdir="/atium/Data/NGS/Raw/170419_fragillis_wtimp1_135818/FASTQ/")

meta=bind_rows(meta, temp)

##3rd run
temp=read_csv("/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135820/wtimp1_135820.csv") %>%
    rename(idx=Index) %>%
    rename(samp=SM_Tag) %>%
    filter(Lane==1) %>%
    select(samp, idx) %>%
    mutate(run=3, paired=T, rawdir="/atium/Data/NGS/Raw/170419_fragillis_wtimp1_135820/FASTQ/")

meta=bind_rows(meta, temp)

##Set parameters for each run, first init
meta=meta %>%
    mutate(time=0, drug="blank")

##Then set times based on searching for 24hr or 48hr in samp name.  Nxt time pls have clearer samp names so that this is easier
meta$time[grepl("24hr", meta$samp)]=24
meta$time[grepl("48hr", meta$samp)]=48

##Drugs
meta$drug[grepl("bft1", meta$samp)]="bft1"
meta$drug[grepl("bft2", meta$samp)]="bft2"

##Ok - get atac seq fastq packaged up?





##Ok - start atac-seq pipeline

if (FALSE) {

    ##loop through all conditions
    for (t in c(0, 24, 48)) {
        for (d in c("blank", "bft2")) {
            condition=filter(meta, time==t, drug==d)           

            ##ok - 3 replicates, so
            for (r in 1:3) {
                read1=list.files(path=condition$rawdir[r], pattern=paste0("*", condition$idx[r], "_1.fastq.gz"), full.names=T)
                system(paste0("cat ", paste(read1, collapse=" "), " >/tmp/rep", r, "_read1.fq.gz"))
                read2=list.files(path=condition$rawdir[r], pattern=paste0("*", condition$idx[r], "_2.fastq.gz"), full.names=T)
                system(paste0("cat ", paste(read1, collapse=" "), " >/tmp/rep", r, "_read2.fq.gz"))
            }
                
            


            
        }
    }
}


