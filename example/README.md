Run:

```
../src/random.pl > random.bed
python ../src/coverage_dev.py --bam NA12878_rmdup.bam --bed random.bed --output random --reference genome.fa
```

* outputs (top 3 lines):

```
==> example/gc_bias.csv <==
,gc_content,mean_coverage
0,0.4875,4.9391
1,0.5257,4.2162

==> example/global_stats.csv <==
,mean_coverage,median_coverage,std_coverage
0,4.9771412882787756,5.0,2.3473035118324352

==> example/percentiles.csv <==
,p10,p25,p50,p75,p90
0,2.0,3.0,5.0,7.0,8.0

==> example/random.bed <==
chr1    860368  870368
chr1    4510135 4520135
chr1    8312956 8322956

==> example/regions.csv <==
,chrom,start,end,mean_coverage,median_coverage,std_coverage
0,chr1,860368,870368,4.9391,5.0,2.428660369421793
1,chr1,4510135,4520135,4.2162,4.0,2.2004221322282684

==> example/regions_high.csv <==
,chrom,start,end,mean_coverage,median_coverage,std_coverage
28,chr1,85434652,85444652,6.234,6.0,2.4378769452127806
92,chr2,32529632,32539632,6.1144,6.0,2.3733336554306894

==> example/regions_low.csv <==
,chrom,start,end,mean_coverage,median_coverage,std_coverage
2,chr1,8312956,8322956,3.715,3.0,2.3888438626247637
37,chr1,113064858,113074858,3.9336,4.0,2.841512104496477
```

