
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("dplyr")
library("dplyr")

install.packages("devtools")
devtools::install_github("pachterlab/sleuth")
library("sleuth")

###First, I must tell the sleuth program where to find my data and prepare my annotation data

#here, I am telling the program where to find my files
sample.id <- dir(file.path("~/kallisto-0.43.1","output/rerun/sleuthformat"))
kal.dirs <- file.path("~/kallisto-0.43.1","output/rerun/sleuthformat",sample.id,"abundance.h5")

#here, I am reading in a txt file that tells the program  what experimental condition each file corresponds to 
s2c <- read.table("~/kallisto-0.43.1/rnaseq_sampleinfo.txt", header=TRUE, stringsAsFactors=FALSE)
s2c <- s2c[ ,c(1,2)]
s2c <- dplyr::mutate(s2c, path= kal.dirs)

#subset data to just get 24hr expression data
s2c.24 <- s2c[s2c$condition != "48hrbft2", ] #this is how i subset my data to compare certain groups
s2c.24 <- s2c.24[s2c.24$condition != "48blank", ] #this is how i subset my data to compare certain groups
s2c.24 <- s2c.24[s2c.24$condition != "time0", ] #this is how i subset my data to compare certain groups

#subset data to just get 48hr expression data
s2c.48 <- s2c[s2c$condition != "24hrbft2", ] #this is how i subset my data to compare certain groups
s2c.48 <- s2c.48[s2c.48$condition != "24blank", ] #this is how i subset my data to compare certain groups
s2c.48 <- s2c.48[s2c.48$condition != "time0", ] #this is how i subset my data to compare certain groups

s2c.48.no031 <- s2c.48[s2c$sample != "CGTACG", ] #REMOVING 48HRBFT2 031 SAMPLE BECAUSE IT IS AN OUTLIER
s2c.48.no031 <- s2c.48[s2c.48$sample != "GTCCGC", ] # REMOVING 48HRBLANK 031 SAMPLE BECAUSE IT IS AN OUTLIER

                                    
#this lets me annotate my transcript data by adding gene names

#source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library("biomaRt")
library("sleuth")
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl" , host = 'ensembl.org')#MUST LOAD SLEUTH BEFORE LOADING THIS... for some reason, this is currently not working
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)



###### Now I will analyze my data



##This code will give me gene names with my gene level analysis, but will not allow me to plot my expression data via plot_bootstraps 
so9 <- sleuth_prep(s2c.24 ~condition, target_mapping=t2g, aggregation_column= 'ext_gene', extra_bootstrap_summary=TRUE)#the ext gene here is what gives you gene names.  
so9<- sleuth_fit(so5,~condition,'full')
so9 <- sleuth_wt(so5, which_beta="condition48hrbft2")

sleuth_table9 <- sleuth_results(so9,"condition24hrbft2", show_all=FALSE)
sleuth_significant9 <- dplyr::filter(sleuth_table9, pval <= 0.01)
sleuth_significant9 <- sleuth_significant9[order(sleuth_significant9$pval), ]
head(sleuth_significant9)
write.csv(sleuth_significant9, "~/Data/rnaseq_data/gene.level.analysis/rerun/24hr.expression.changes.genelevel.rerun.pval01.csv")


so5 <- sleuth_prep(s2c.48, ~condition, target_mapping=t2g, aggregation_column= 'ext_gene', extra_bootstrap_summary=TRUE)#the ext gene here is what gives you gene names.  
so5<- sleuth_fit(so5,~condition,'full')
so5 <- sleuth_wt(so5, which_beta="condition48hrbft2")

sleuth_table5 <- sleuth_results(so5,"condition48hrbft2", show_all=FALSE)
sleuth_significant5 <- dplyr::filter(sleuth_table5, pval <= 0.01)
sleuth_significant5 <- sleuth_significant5[order(sleuth_significant5$pval), ]
head(sleuth_significant5)
write.csv(sleuth_significant5, "~/Data/rnaseq_data/gene.level.analysis/rerun/48hr.expression.changes.genelevel.rerun.pval01.csv")


#this code will give me gene level expression data that I can use to create a table of tpm values for all samples
so6 <- sleuth_prep(s2c, extra_bootstrap_summary=TRUE, aggregation_column="ext_gene", target_mapping=t2g)
so6<- sleuth_fit(so6, ~condition, 'full')
so6 <- sleuth_fit(so6, ~1, 'reduced')
so6 <- sleuth_lrt(so6, 'reduced', 'full')

sleuth_table <- sleuth_results(so6, 'reduced:full', 'lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
sleuth_significant <- sleuth_significant[order(sleuth_significant$pval), ]
head(sleuth_significant, 50)
head(sleuth_significant)


counts=kallisto_table(so6) #creates count table from my data

write.csv(counts, "~/Data/rnaseq_data/gene.level.analysis/rerun/tpmcounts.sleuth.genelevel.allsamples.rerun.csv")



###################################



# This will perform transcript level differential expression and give you gene names associated with those transcripts
so3 <- sleuth_prep(s2c.24, ~condition, target_mapping=t2g, extra_boostrap_summary=TRUE) 
so3<- sleuth_fit(so3, ~condition, 'full')
so3 <- sleuth_wt(so3, which_beta="condition24hrbft2")

sleuth_table3 <- sleuth_results(so3,"condition24hrbft2", show_all=FALSE)
sleuth_significant3 <- dplyr::filter(sleuth_table3, pval <= 0.05)
sleuth_significant3 <- sleuth_significant3[order(sleuth_significant3$pval), ]
head(sleuth_significant3)


#this is the code used to get a value called 'b' which is a biased logfc (the sleuth_wt function does this) and to perform gene level diferential expression analysis. With this, I can also look at expression data for each gene. The key to getting that data is to have extra_bootstrap_summary in the code. This will not give me gene names though, only gene IDs
so2 <- sleuth_prep(s2c.24, ~condition, target_mapping=t2g, aggregation_column= 'ens_gene', extra_bootstrap_summary=TRUE)#here, use s2c.24 or s2c.48 
so2<- sleuth_fit(so2,~condition,'full')
so2 <- sleuth_wt(so2, which_beta="condition24hrbft2")

sleuth_table2 <- sleuth_results(so2,"condition24hrbft2", show_all=FALSE)
sleuth_significant2 <- dplyr::filter(sleuth_table2, pval <= 0.05)
sleuth_significant2 <- sleuth_significant2[order(sleuth_significant2$qval), ]
head(sleuth_significant2)

#this lets me visualize the rna data for each group and sample. When run on Rstudio, I Can run sleuth_live(so2) in order to look at the data in an interactive way. This does not work on Bender.
pdf("sleuth.test2.pdf")
plot_bootstrap(so2,"ENSG00000167964", units = "scaled_reads_per_base", color_by = "condition")

dev.off()

