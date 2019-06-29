
###unique bft2 peaks (BFT2 opened peaks)

bedtools intersect -v -bed -a ~/atac_dnase_pipelines/rerun/out.24hrbft2/peak/macs2/idr/optimal_set/out.24hrbft2_ppr.IDR0.1.filt.narrowPeak.gz -b ~/atac_dnase_pipelines/rerun/out.24hrblank/peak/macs2/idr/optimal_set/out.24hrblank_ppr.IDR0.1.filt.narrowPeak.gz > /home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrbft2.narrowPeak



bedtools intersect -v -bed -a ~/atac_dnase_pipelines/rerun/out.48hrbft2/peak/macs2/idr/optimal_set/out.48hrbft2_ppr.IDR0.1.filt.narrowPeak.gz -b ~/atac_dnase_pipelines/rerun/out.48hrblank/peak/macs2/idr/optimal_set/out.48hrblank_ppr.IDR0.1.filt.narrowPeak.gz > /home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrbft2.narrowPeak


### unique blank peaks (BFT2 closed peaks)


bedtools intersect -v -bed -a ~/atac_dnase_pipelines/rerun/out.24hrblank/peak/macs2/idr/optimal_set/out.24hrblank_ppr.IDR0.1.filt.narrowPeak.gz -b ~/atac_dnase_pipelines/rerun/out.24hrbft2/peak/macs2/idr/optimal_set/out.24hrbft2_ppr.IDR0.1.filt.narrowPeak.gz > /home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.24hrblank.narrowPeak



bedtools intersect -v -bed -a ~/atac_dnase_pipelines/rerun/out.48hrblank/peak/macs2/idr/optimal_set/out.48hrblank_ppr.IDR0.1.filt.narrowPeak.gz -b ~/atac_dnase_pipelines/rerun/out.48hrbft2/peak/macs2/idr/optimal_set/out.48hrbft2_ppr.IDR0.1.filt.narrowPeak.gz > /home/jallen66/Data/atac/final.analysis/rerun/uniquepeaks.optimalset.48hrblank.narrowPeak


	 
