module load conda
conda activate macs3

# Peak calling with callpeak
macs3 callpeak -f BAMPE --call-summits -t ENCFF045OAB.chr14.proc_rh.nsort.bam -g 107043718 -n ENCFF045OAB.macs3.default.summits.bampe -B -q 0.05

# Peak calling with hmmratac
macs3 hmmratac -b ../genrich/ENCFF045OAB.chr14.proc_rh.nsort.bam -n ENCFF045OAB.macs3.hmmratac.bampe





















