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
    ##Distribute which set to 3 different computers by changing the starting index
    for (i in seq(1, dim(dataloc)[1], 3)) {
        
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
        ##started at 8:44AM
        system(paste0("samtools index ", dataloc$label[i], ".bam"))
        
    }
}

nucchr=paste(c(paste0("chr", 1:22), "chrX", "chrY"), collapse=" ")

if (TRUE) {
    ##remove chrM
    for (i in 1:dim(dataloc)[1]) {
        #system(paste0("samtools idxstats ", dataloc$bam[i], " | cut -f 1 | grep -v chrM | xargs samtools view -b ",
        #              dataloc$bam[i], " >", dataloc$label[i], "noM.bam"))
        
        system(paste0("samtools view -b ", dataloc$bam[i], " ", nucchr, " >", dataloc$label[i], "noM.bam"), wait=F)
    }

    ##remove dups
    for (i in 1:dim(dataloc)[1]) {
        system(paste0("java -Xmx8G -jar ~/Code/picard/picard.jar MarkDuplicates ",
                      "REMOVE_DUPLICATES=true INPUT=", dataloc$label[i], "noM.bam", 
                      " OUTPUT=", dataloc$label[i], "nodup.bam",
                      " METRICS_FILE=", dataloc$label[i], "_duplicate_metrics.txt"))
    }

    for (i in 1:dim(dataloc)[1]) {

        system(paste0("samtools index ", dataloc$label[i], "nodup.bam"), wait=F)
    }
    
}

##Generate and share coverage plots

if (TRUE) {

    for (i in 1:dim(dataloc)[1]) {

        system(paste0("/bin/bash -c ", shQuote(paste0("source activate fragillis; bamCoverage -b ",
                                                      dataloc$label[i], "nodup.bam -bs 100 -o ",
                                                      dataloc$label[i], ".bw"))))
    }


    
    for (cond in unique(dataloc$condition)) {
        just.these=dataloc %>%
            filter(condition==cond)
        
        system(paste0("samtools merge -r ", cond, ".merged.bam ", paste(paste0(just.these$label, "nodup.bam"), collapse=" ")), wait=F)
        system(paste0("samtools index ", cond, ".merged.bam"))
        system(paste0("/bin/bash -c ", shQuote(paste0("source activate fragillis; bamCoverage -b ",
                                                      cond, ".merged.bam -bs 100 -o ", cond, ".bw"))))

    }
    
}


if (TRUE) {

    ##Merge bam and call peaks on it
    system(paste0("samtools merge merged.bam ", paste(paste0(dataloc$label, "nodup.bam"), collapse=" ")))


    
    ##call peaks
    system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                  "--nomodel -g hs -t merged.bam --broad -n broadmerged --keep-dup all"))),
           wait=F)


    system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                  "--nomodel -g hs -t merged.bam -n sharpmerged --keep-dup all"))),
           wait=F)


    ##per condition
    for (cond in unique(dataloc$condition)) {
        just.these=dataloc %>%
            filter(condition==cond)

        system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                      "--nomodel -g hs -t ", cond, ".merged.bam --broad -n ",
                                                      cond, ".broadmerged", " --keep-dup all"))),
               wait=F)

        system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                      "--nomodel -g hs -t ", cond, ".merged.bam -n ",
                                                      cond, ".sharpmerged", " --keep-dup all"))),
               wait=F)

    }

    ##individ
    for (i in 1:dim(dataloc)[1]) {
        system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                      "--nomodel -g hs -t ", dataloc$label[i], "nodup.bam --broad -n ",
                                                      dataloc$label[i], ".broadmerged", " --keep-dup all"))),
               wait=F)

        system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27.fragillis; macs2 callpeak ",
                                                      "--nomodel -g hs -t ", dataloc$label[i], "nodup.bam -n ",
                                                      dataloc$label[i], ".sharpmerged", " --keep-dup all"))),
               wait=F)

    }
    ##Denny et al command
    ##macs2 callpeak --nomodle --broad --keep-dup all

    ##Jawara command
    ##macs2 callpeak --nomodel -t ${cleaned} -n ${peakout} --nolambda --keep-dup all --shift -75 --extsize 150 -p 0.1 -B -SPMR

}

