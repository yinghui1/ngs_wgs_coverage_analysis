ID=NA12878_rmdup
bedtools genomecov -ibam $ID.bam -bg > $ID.bedgraph
LC_ALL=C sort -S 30% --parallel=8 -k1,1 -k2,2n $ID.bedgraph > ${ID}_srt.bedgraph
[bli@silencer fg]$ bedGraphToBigWig ${ID}_srt.bedgraph chrom.sizes $ID.bw
