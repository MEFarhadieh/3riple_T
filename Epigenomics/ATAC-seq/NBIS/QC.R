
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(Rsamtools)

# Read data
bamFile="ENCFF045OAB.chr14.blacklist_M_filt.mapq5.dedup.bam"
bamFileLabels <- "ENCFF045OAB"

### Shiftig and Splitting Aligned Reads
# Collect library statistics made by QC.sh
bam_qc=bamQC(bamFile, outPath = NULL)
bam_qc[1:10]

outPath <- "splitBam"
dir.create(outPath)

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

bamTop100 <- scanBam(BamFile(bamFile, yieldSize = 100),
                     param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]

# shift the coordinates only for alignments on chr14
seqlev <- "chr14"
which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
gal <- readBamFile(bamFile, tag=tags, which=which,asMates=TRUE, bigFile=TRUE)

shiftedBamFile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamFile)

saveRDS(gal1, file = "gal1.rds", ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)

# Find genomic locations of TSS
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
txs <- txs[seqnames(txs) %in% "chr14"]
genome <- Hsapiens

# Split the alignments
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)

saveRDS(objs, file = "objs.rds", ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)

### Signal in NFR and Mononucleosome Fractions
# Use TSS coordinates
bamFiles <- file.path(outPath,
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

# log2 transform the signal around TSS
librarySize <- estLibSize(bamFiles)

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")],
                          TSS=TSS,
                          librarySize=librarySize,
                          seqlev=seqlev,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

# Heatmap of ATAC-seq signal in NFR and mononculeosome fractions                  
pdf("Heatmap_splitbam.pdf")
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)

dev.off()

### Signal at TSS
out <- featureAlignedDistribution(sigs,
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l",
                                  ylab="Averaged coverage")

# rescale the nucleosome-free and nucleosome signals to 0~1 for plotting
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)

# Plot of ATAC-seq signal in NFR and mononculeosome fractions
pdf("TSSprofile_splitbam.pdf")
        matplot(out, type="l", xaxt="n",
        xlab="Position (bp)",
        ylab="Fraction of signal")
        axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
        abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
                  
