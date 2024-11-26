# Histogram of ATAC-seq signal
module load picard/2.23.4

java -Xmx32G -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics \
 -I ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam \
 -O ENCFF045OAB.chr14.proc.fraglen.stats \
 -H ENCFF045OAB.chr14.proc.fraglen.pdf -M 0.5

# Fraction of Reads in Peaks
awk -F $'\t' 'BEGIN {OFS = FS}{ $2=$2+1; peakid="macs3Peak_"++nr;  print peakid,$1,$2,$3,"."}' ENCFF045OAB.macs3.default.summits.bampe_peaks.narrowPeak > ENCFF045OAB.chr14.macs3.saf
module load subread/2.0.0
featureCounts -p -F SAF -a ENCFF045OAB.chr14.macs3.saf --fracOverlap 0.2 -o ENCFF045OAB.peaks_macs3.counts ENCFF045OAB.chr14.proc_rh.nsort.bam
head ENCFF045OAB.chr14.macs3.saf

# Summarise reads
module load subread/2.0.0
featureCounts -p -F SAF -a ENCFF045OAB.chr14.macs3.saf --fracOverlap 0.2 -o ENCFF045OAB.peaks_macs3.counts ENCFF045OAB.chr14.proc_rh.nsort.bam

featureCounts -p -F SAF -a nk_merged_peaks.saf --fracOverlap 0.2 -o nk_merged_peaks_macs3.counts ENCFF363HBZ.chr14.proc.bam ENCFF398QLV.chr14.proc.bam ENCFF828ZPN.chr14.proc.bam ENCFF045OAB.chr14.proc.bam
head nk_merged_peaks_macs3.counts

# Prepare for R
awk '(NR>1)' nk_merged_peaks_macs3.counts > nk_merged_peaks_macs3.counts.tsv
