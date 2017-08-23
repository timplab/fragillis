library(tidyverse)

##Ok - first pull in csv files giving locations of rna data


##First one doesn't have info, but I pulled an excel from genesifter, made a csv
meta=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/129658.csv", skip=2) %>%
    rename(idx=`Index Sequence`) %>%
    rename(samp=Template) %>%
    select(samp, idx) %>%
    mutate(run=1, paired=T, rawdir="/mithril/Data/NGS/Raw/160823_bft_rna/wtimp1_129656/FASTQ/")

##continue to build full filename list

temp=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/135814.csv", skip=2) %>%
    rename(idx=`Index Sequence`) %>%
    rename(samp=Template) %>%
    select(samp, idx) %>%
    mutate(run=2, paired=F, rawdir="/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135814/FASTQ/") 

meta=bind_rows(meta, temp)

##3rd run
temp=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/135816.csv", skip=2) %>%
    rename(idx=`Index Sequence`) %>%
    rename(samp=Template) %>%
    select(samp, idx) %>%
    mutate(run=3, paired=F, rawdir="/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135816/FASTQ/") 

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


##Setup reference locations
ref="/mithril/Data/NGS/Reference/human_hisat38/genome_tran"
gtfref="/mithril/Data/NGS/Reference/human_hisat38/hg38_ucsc.annotated.gtf"

##Setup output locations
outdir="/atium/Data/NGS/Aligned/170613_fragillis/"
meta$bamout=NA


if (TRUE) {

    ##loop through all samples
    #for (i in 1:dim(meta)[1]) {
    for (i in 1:1) { 

        meta$bamout[i]=file.path(outdir, "aligned", paste0(meta$samp[i], ".sorted.bam"))
        ##If paired, do paired alignment, otherwise, single end
        if (meta$paired[i]) {
            ##Get reads, cat, and trim in tmp dir          
            read1=list.files(path=meta$rawdir[i], pattern=paste0("*", meta$idx[i], "_1.fastq.gz"), full.names=T)
            read2=list.files(path=meta$rawdir[i], pattern=paste0("*", meta$idx[i], "_2.fastq.gz"), full.names=T)

            system(paste0("cat ", paste(read1, collapse=" "), " >/tmp/read1.fq.gz"))
            system(paste0("cat ", paste(read2, collapse=" "), " >/tmp/read2.fq.gz"))
            
            system(paste0("atropos --aligner insert -T 6 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ",
                          "-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ",
                          "-q 20 -o /tmp/read1.trim.fq.gz -p /tmp/read2.trim.fq.gz -pe1 /tmp/read1.fq.gz -pe2 /tmp/read2.fq.gz"))

            system(paste0("hisat2 --threads 10 --tmo -x ", ref ," -1 ", "/tmp/read1.trim.fq.gz", " -2 ", "/tmp/read2.trim.fq.gz",
                          " | samtools view -b - | samtools sort -o ", meta$bamout[i]))

            
#            system(paste0("hisat2 --threads 10 --tmo -x ", ref ," -1 ", read1, " -2 ", read2,
#                          " | samtools view -b - | samtools sort -o ", meta$bamout[i]))
            
        } else {
            read1=paste(list.files(path=meta$rawdir[i], pattern=paste0("*", meta$idx[i], "_1.fastq.gz"), full.names=T), collapse=",")
            system(paste0("hisat2 --threads 10 --tmo -x ", ref ," -U ", read1, 
                          " | samtools view -b - | samtools sort -o ", meta$bamout[i]))
        }
    }

}

##Stringtie

meta$ballout=NA

if (TRUE) {

    ##for (i in 1:dim(meta)[1]) {
    for (i in 1:1) {
        meta$ballout[i]=file.path(outdir, "ballgown", paste0(meta$samp[i], ".gtf"))
        
        system(paste0("stringtie -e -B -p 10 -G ", gtfref, " -o ", meta$ballout[i], " ", meta$bamout[i]))
    }

}


