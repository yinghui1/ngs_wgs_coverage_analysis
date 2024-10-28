bwa mem -t 8 Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/genome.fa fgs1.fastq.gz fgs2.fastq.gz | samtools view - -Sb | samtools sort - > fg.bam
