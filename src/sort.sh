zcat fg2.fastq.gz | LC_ALL=C sort --parallel=4 | sed 's/\t/\n/g' | gzip > fgs2.fastq.gz &
zcat fg1.fastq.gz | LC_ALL=C sort --parallel=4 | sed 's/\t/\n/g' | gzip > fgs1.fastq.gz &
