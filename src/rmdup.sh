# https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar

samtools addreplacerg -r "@RG\tID:RG1\tSM:NA12878\tPL:Illumina\tLB:NA12878.fa" -o NA12878rg.sam NA12878.bam
java -jar picard.jar MarkDuplicates     I=NA12878rg.sam     O=NA12878_rmdup.bam     M=NA12878_duplication_metrics.txt     REMOVE_DUPLICATES=true
