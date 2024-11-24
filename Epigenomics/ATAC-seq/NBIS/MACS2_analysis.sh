module load MACS/2.2.6
module load BEDTools/2.25.0

bedtools bamtobed -i ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam >ENCFF045OAB.chr14.bed

macs2 callpeak -t ENCFF045OAB.chr14.bed \
-n ENCFF045OAB.chr14.macs2.encode -f BED \
-g 107043718 -q 0.05 --nomodel --shift -75 --extsize 150 \
--call-summits --keep-dup all
