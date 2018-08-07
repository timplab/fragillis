#!/bin/bash

rawdir=/atium/Data/NGS/Aligned/170423_atacseq_JA/raw_data_renamed
outdir=/atium/Data/NGS/Aligned/170423_atacseq_JA/aligned_220

for samp in "24hrbft2"  "24hrblank" "48hrbft2" "48hrblank" "time0" 
do
    echo $samp

    read1=${rawdir}/${samp}_220_read1.fastq.gz
    read2=${rawdir}/${samp}_220_read2.fastq.gz
    
    echo ${read1}
    echo ${read2}
    
    ##Bowtie2 alignment
    sorted=${outdir}/${samp}.sorted

    # in the paper, they noted that they used -X 2000 because of the possibility of long fragments (no shearing)
    bowtie2 -X 2000 -p 10 -t  \
        -x /mithril/Data/NGS/Reference/human38/GRCH38 \
    -1 ${read1} -2 ${read2} \
   2> ${outdir}/${samp}_bowtie2.log\
   | samtools view -bS -| samtools sort - -o ${sorted}.bam
    echo ${sorted}
  ##samtools: 1) convert sam to bam, 2) sort, 3) index
  samtools index ${sorted}.bam

done
