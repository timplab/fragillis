#!/bin/bash

##Move all jawara raw concat fq to Aligned dir to make life easier.

rawdir1=/mithril/Data/NGS/Raw/160823_bft_rna/wtimp1_129656/FASTQ/
rawdir2=/atium/Data/NGS/Raw/170419_fragillis/wtimp1_135816/FASTQ/

##Get gencode reference
##wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz


ref=/mithril/Data/NGS/Reference/human_genes/gencode.v27.transcripts.fa.gz
gtfref=/mithril/Data/NGS/Reference/human_hisat38/hg38_ucsc.annotated.gtf


outdir=/atium/Data/NGS/Aligned/170613_fragillis/



if [ 0 -eq 0 ]; then
    #for samp in "Time0"
    for samp in "Time0" "24hrbft1" "24hrbft2" "24hrblank" "48hrbft1" "48hrbft2" "48hrblank"
    do
    
       hisat2 --threads 10 --tmo -x ${ref} -1 ${rawdir1}/${samp}_1.fastq.gz -2 ${rawdir1}/${samp}_2.fastq.gz |
	    samtools view -b - | samtools sort -o ${outdir}/${samp}.sorted.bam
       
    done
    
fi


if [ 0 -eq 0 ]; then
    for samp in "Time0"
    do
	
	echo "Stringtie"
	stringtie -e -B -p 10 -G ${gtfref} -o ${outdir}/ballgown/${samp}.gtf ${outdir}/${samp}.sorted.bam

    done

fi


##for samp in "24hrbft2.031" "24hrblank.031"  "48hrbft2.031" "48hrblank.031" "time0.031" "24hrbft2.220" "24hrblank.220"\
##			   "48hrbft2.220" "48hrblank.220" "time0.220" "24hrbft2.trial1" "24hrblank.trial1" "48hrbft2.trial1" "48hrblank.trial1" "time0.trial1"




##for samp in "24hrbft2.031" "24hrblank.031"  "48hrbft2.031" "48hrblank.031" "time0.031" "24hrbft2.220" "24hrblank.220"\
##			   "48hrbft2.220" "48hrblank.220" "time0.220" "24hrbft2.trial1" "24hrblank.trial1" "48hrbft2.trial1" "48hrblank.trial1" "time0.trial1"


