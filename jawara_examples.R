library(tidyverse)
library(readxl)

workdir="~/Dropbox/Data/Genetics/Fragillis/170725_jallen/"

upeaks=read_excel(file.path(workdir, "unique.peaks.all.dedup.xlsx"))
express=read_csv(file.path(workdir, "24hr.expression.changes_final.csv")) %>%
    rename(gene.names=`gene names`) 




express=express %>%
    mutate(bft24hr.chromatin=ifelse( (gene.names %in% z), yes="opened", no="nochange")) %>%

express$bft24hr.chromatin="nochange"
express$bft24hr.chromatin[express$gene.names %in% upeaks$`unique 24hrbft2 peaks`]="opened"
express$bft24hr.chromatin[express$gene.names %in% upeaks$`unique 24hrblank peaks`]="closed"

table(express$bft24hr.chromatin)

express$bft48hr.chromatin="nochange"
express$bft48hr.chromatin[express$gene.names %in% upeaks$`unique 48hrbft2 peaks`]="opened"
express$bft48hr.chromatin[express$gene.names %in% upeaks$`unique 48hrblank peaks`]="closed"

table(express$bft48hr.chromatin)

