library(tidyverse)
library(googlesheets)

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

whichcomp=3

for (i in seq(whichcomp, dim(dataloc)[1], 3)) {

    system(paste0("java -jar ~/Code/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 6 ",
                  dataloc$ill.r1[i], " ", dataloc$ill.r2[i],
                  " /tmp/read1.trim.fq.gz /tmp/read1.unpaired.fq.gz",
                  " /tmp/read2.trim.fq.gz /tmp/read2.unpaired.fq.gz",
                  " ILLUMINACLIP:/home/timp/Code/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10",
                  " LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40"))
    

    system(paste0("bowtie2 -X 2000 -p 10 -t -x /mithril/Data/NGS/Reference/human38/GRCH38 ",
                  "-1 /tmp/read1.trim.fq.gz -2 /tmp/read2.trim.fq.gz ",
                  "2> ", workdir, "/", dataloc$label[i], "_bowtie2.log ",
                  "| samtools view -bS - | samtools sort - -o ", dataloc$label[i], ".bam"))

    system(paste0("samtools index ", dataloc$label[i], ".bam"))

}
