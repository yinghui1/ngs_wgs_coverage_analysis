#bwa mem -t 8 Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0/genome.fa fgs1.fastq.gz fgs2.fastq.gz | samtools view - -Sb > NA12878_pe.bam

#bedtools bamtobed -bedpe -i NA12878_pe.bam > NA12878.bedpe

#grep -v "\-1" NA12878.bedpe | sort -k1,1 -k2,2n | bedtools merge -i - > NA12878_pe.bed

# Calculate coverage over merged insert intervals
## bedtools coverage -a NA12878_pe.bed -b NA12878.bam > NA12878_coverage.txt

mkfifo NA12878_rmdup_srt12.bam

samtools view -h NA12878_rmdup.bam | sort -S 30% -k1,1 -k2,2n --parallel=8 | samtools view -Sb - > NA12878_rmdup_srt12.bam &

bedtools bamtobed -bedpe -i NA12878_rmdup_srt12.bam > NA12878_rmdup_srt12.bedpe

grep -v "\-1" NA12878_rmdup_srt12.bedpe | sort -S 30% -k1,1 -k2,2n --parallel=8 | bedtools merge -i - > NA12878_rmdup_srt12_pe.bed
