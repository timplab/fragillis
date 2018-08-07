library(tidyverse)

cosmic=read_tsv("/atium/Data/NGS/Aligned/170613_fragillis/cosmic/CosmicGenomeScreensMutantExport.tsv") %>%
    separate(`Mutation genome position`, into=c("chr", "start", "end"), sep=[:-])

