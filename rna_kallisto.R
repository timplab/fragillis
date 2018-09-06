library(tidyverse)
library(googlesheets)

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

if (FALSE) {

    ##Put file locations into gsheet
    
    ##First one doesn't have info, but I pulled an excel from genesifter, made a csv
    meta=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/129658.csv", skip=2) %>%
        rename(idx=`Index Sequence`) %>%
        rename(samp=Template) %>%
        select(samp, idx) %>%
        mutate(replicate=1, paired=T, rawdir="/mithril/Data/NGS/Raw/160823_bft_rna/wtimp1_129656/FASTQ/")
    
    ##continue to build full filename list
    
    temp=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/135814.csv", skip=2) %>%
        rename(idx=`Index Sequence`) %>%
        rename(samp=Template) %>%
        select(samp, idx) %>%
        mutate(replicate=2, paired=F, rawdir="/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135814/FASTQ/") 
    
    meta=bind_rows(meta, temp)
    
    ##3rd run
    temp=read_csv("~/Dropbox/Data/Genetics/Fragillis/170716_pipe/135816.csv", skip=2) %>%
        rename(idx=`Index Sequence`) %>%
        rename(samp=Template) %>%
        select(samp, idx) %>%
        mutate(replicate=3, paired=F, rawdir="/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135816/FASTQ/") 
    
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

    meta=meta %>%
        filter(drug!="bft1") %>%
        mutate(bft=ifelse(drug=="bft2", "y", "n")) %>%
        select(-paired, -samp)

    dataloc=left_join(dataloc, meta)

    for (i in 1:dim(dataloc)[2]) {       
        ##Loop this
        
        rna.files=list.files(path=dataloc$rawdir[i], pattern=paste0("*", dataloc$idx[i], "_1.fastq.gz"), full.names=T)
        
        system(paste0("zcat ", paste(rna.files, collapse=" "), " | gzip >/tmp/rna.fq.gz"))

        dataloc$rna.fq[i]=file.path(workdir, "rna", paste0(dataloc$label[i], ".fq.gz"))
        
        system(paste0("java -jar ~/Code/Trimmomatic-0.38/trimmomatic-0.38.jar SE -phred33 -threads 6 ",
                      "/tmp/rna.fq.gz ", dataloc$rna.fq[i],
                      " ILLUMINACLIP:/home/timp/Code/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10",
                      " LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:40"))
        
    }
    
    dataloc=select(dataloc, -idx, -rawdir, -drug)

    colnum=which(colnames(dataloc)=="rna.fq")

    gs_edit_cells(fullsheet, ws="180608_dat", input=c("rna.fq", dataloc$rna.fq), anchor=paste0("R1C", colnum))

}


##Run kallisto

if (FALSE) {

    ref="/mithril/Data/NGS/Reference/human_ensembl93/ensembl93.cdna"
    
    for (i in 1:dim(dataloc)[2]) {       

        system(paste0("/bin/bash -c ", shQuote(paste0("source activate fragillis; kallisto quant -i ",
                                                      ref, " -o ", dataloc$label[i], ".kallisto -b 100 ",
                                                      "--single -l 200 -s 20 -t 8 ", dataloc$rna.fq[i]))))

    }

}
