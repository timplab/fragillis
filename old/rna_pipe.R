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
##Get gencode reference


ref="/mithril/Data/NGS/Reference/human_genes/gencode.v27.salmon"
##Make index
if (FALSE) {
    ##wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz
    system(paste0("~/Code/salmon/bin/salmon index --gencode -p 10 ",
                  "-t /mithril/Data/NGS/Reference/human_genes/gencode.v27.transcripts.fa.gz -i ",
                  ref))
}
    

    
##Setup output locations
outdir="/atium/Data/NGS/Aligned/170613_fragillis/salmon/"


if (TRUE) {

    meta$salmon.out=file.path(outdir, paste0(meta$samp, ".quant"))
    
    ##loop through all samples
    for (i in 1:dim(meta)[1]) {
    #for (i in 1:1) { 


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


            system(paste0("~/Code/salmon/bin/salmon quant -i ", ref," -l A",
                          " -1 /tmp/read1.trim.fq.gz ",
                          " -2 /tmp/read2.trim.fq.gz ",
                          "-p 8 -o ", meta$salmon.out[i]))

            system("rm /tmp/read*.gz")
            
        } else {
            read1=list.files(path=meta$rawdir[i], pattern=paste0("*", meta$idx[i], "_1.fastq.gz"), full.names=T)                        

            system(paste0("cat ", paste(read1, collapse=" "), " >/tmp/read1.fq.gz"))

            system(paste0("atropos -T 6 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ",
                          "-q 20 -o /tmp/read1.trim.fq.gz -se /tmp/read1.fq.gz"))

            system(paste0("~/Code/salmon/bin/salmon quant -i ", ref," -l A",
                          " -r /tmp/read1.trim.fq.gz ",
                          "-p 8 -o ", meta$salmon.out[i]))

            system("rm /tmp/read*.gz")
        }
    }

}

