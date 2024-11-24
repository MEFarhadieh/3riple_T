#Histogram of ATAC-seq signal
module load picard/2.23.4

java -Xmx32G -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics \
 -I ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam \
 -O ENCFF045OAB.chr14.proc.fraglen.stats \
 -H ENCFF045OAB.chr14.proc.fraglen.pdf -M 0.5
