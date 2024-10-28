ID=NA12878_rmdup_srt12
./mkpe.pl $ID.bedpe | sort -S 30% -k1,1 -k2,2n --parallel=8 > ${ID}_mkpe.bed
bedtools genomecov -i NA12878_rmdup_srt12_mkpe.bed -bg -g chrom.sizes > NA12878_rmdup_srt12_mkpe.bedgraph
# bedSort NA12878.bedGraph NA12878s.bedGraph
bedGraphToBigWig ${ID}_mkpe.bedgraph chrom.sizes ${ID}_mkpe.bw
