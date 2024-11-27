library(ATACseqQC)
library(MotifDb)

library(BSgenome.Hsapiens.UCSC.hg38)
genome <- Hsapiens

# Load data from shifted bam file
shifted.bamFile="shifted.bam"

# Limit the analysis to chr14
seqlev <- "chr14"

# Summarise the signal in the vicinity of CTCF motifs
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits=2)

ctcf <- factorFootprints(shifted.bamFile, pfm=CTCF[[1]],
                         genome=genome,
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)

ctcf$spearman.correlation

pdf("ctcf_footprnt.pdf")
sigs <- factorFootprints(shifted.bamFile, pfm=CTCF[[1]],
                         genome=genome,
                         min.score="90%", seqlev=seqlev,
                         upstream=100, downstream=100)
dev.off()

