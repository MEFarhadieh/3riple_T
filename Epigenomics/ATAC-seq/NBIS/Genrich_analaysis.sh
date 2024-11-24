
module load bioinfo-tools
module load samtools/1.8


#subset bam and change header
samtools view -h ENCFF045OAB.chr14.proc.bam chr14 | grep -P "@HD|@PG|chr14" | samtools view -Shbo ENCFF045OAB.chr14.proc_rh.bam
samtools view -h  ENCFF828ZPN.chr14.proc.bam chr14 | grep -P "@HD|@PG|chr14" | samtools view -Shbo  ENCFF828ZPN.chr14.proc_rh.bam


# sort the bam file by read name (required by Genrich)
samtools sort -n -o ENCFF045OAB.chr14.proc_rh.nsort.bam -T sort.tmp  ENCFF045OAB.chr14.proc_rh.bam
samtools sort -n -o ENCFF828ZPN.chr14.proc_rh.nsort.bam -T sort.tmp  ENCFF828ZPN.chr14.proc_rh.bam

# Detect peaks
Genrich -j -t ENCFF045OAB.chr14.proc_rh.nsort.bam  -o ENCFF045OAB.chr14.genrich.narrowPeak
wc -l ENCFF045OAB.chr14.genrich.narrowPeak

Genrich -j -t ENCFF828ZPN.chr14.proc_rh.nsort.bam -o ENCFF828ZPN.chr14.genrich.narrowPeak
# or detect more peaks in the joint mode
Genrich -j -t ENCFF045OAB.chr14.proc_rh.nsort.bam,ENCFF828ZPN.chr14.proc_rh.nsort.bam  -o nk_stim.chr14.genrich.narrowPeak
wc -l *narrowPeak
