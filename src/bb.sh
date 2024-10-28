awk '{print $1"\t"$2"\t"$3}' NA12878_rmdup_srt12_mkpe.bed | LC_ALL=C sort -S 30% --parallel=8 -k1,1 -k2,2n > tmp.bed
bedToBigBed tmp.bed chrom.sizes NA12878_rmdup_srt12_mkpe.bb
