library(edgeR)
library(EDASeq)

library(GenomicAlignments)
library(GenomicFeatures)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(wesanderson)

library(Hmisc)
library(dplyr)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
library(ReactomePA)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

ff = FaFile("/references/hg38.fa")

# read data 
cnt_table = read.table("nk_merged_peaks_macs3.counts", sep="\t", header=TRUE, blank.lines.skip=TRUE)
rownames(cnt_table)=cnt_table$Geneid

# update colnames of this count table
colnames(cnt_table)=c("Geneid","Chr","Start","End","Strand","Length","ENCFF363HBZ","ENCFF398QLV","ENCFF045OAB","ENCFF828ZPN")


groups = factor(c(rep("NK",2),rep("NKstim",2)))

# this data frame contains only read counts to peaks on assembled chromosomes
reads.peak = cnt_table[,c(7:10)]

# GC-aware normalisation
gr = GRanges(seqnames=cnt_table$Chr, ranges=IRanges(cnt_table$Start, cnt_table$End), strand="*", mcols=data.frame(peakID=cnt_table$Geneid))

peakSeqs = getSeq(x=ff, gr)

gcContentPeaks = letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]

# divide into 20 bins by GC content
gcGroups = Hmisc::cut2(gcContentPeaks, g=20)
mcols(gr)$gc = gcContentPeaks

# visualise GC bias in peaks
lowListGC = list()
for(kk in 1:ncol(reads.peak)){
  set.seed(kk)
  lowListGC[[kk]] = lowess(x=gcContentPeaks, y=log1p(reads.peak[,kk]), f=1/10)
}

names(lowListGC)=colnames(reads.peak)

dfList = list()
for(ss in 1:length(lowListGC)){
  oox = order(lowListGC[[ss]]$x)
  dfList[[ss]] = data.frame(x=lowListGC[[ss]]$x[oox], y=lowListGC[[ss]]$y[oox], sample=names(lowListGC)[[ss]])
}
dfAll = do.call(rbind, dfList)
dfAll$sample = factor(dfAll$sample)

p1.1 = ggplot(dfAll, aes(x=x, y=y, group=sample, color=sample)) +
  geom_line(size = 1) +
  xlab("GC-content") +
  ylab("log(count + 1)") +
  theme_classic()

pdf("GCcontent_peaks.pdf")
p1.1
dev.off()


reads.peak=as.matrix(reads.peak)

dataOffset = withinLaneNormalization(reads.peak,y=gcContentPeaks,num.bins=20,which="full",offset=TRUE)
dataOffset = betweenLaneNormalization(reads.peak,which="full",offset=TRUE)

# DA peaks statistics with edgeR
design = model.matrix(~groups)

d = DGEList(counts=reads.peak, group=groups)

keep = filterByExpr(d)
summary(keep)

d=d[keep,,keep.lib.sizes=FALSE]

d$offset = -dataOffset[keep,]
d.eda = estimateGLMCommonDisp(d, design = design)
d.eda = estimateGLMCommonDisp(d, design = design)
fit = glmFit(d.eda, design = design)
lrt.EDASeq = glmLRT(fit, coef = 2)

DA_res=as.data.frame(topTags(lrt.EDASeq, nrow(lrt.EDASeq$table)))
head(DA_res)

# add more peak information
DA_res$Geneid = rownames(DA_res)
DA.res.coords = left_join(DA_res,cnt_table[1:4],by="Geneid")
head(DA.res.coords)

write.table(DA.res.coords, "nk_DA_stim_vs_ctrl.tsv", quote = FALSE, sep = "\t",
    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
    col.names = TRUE, fileEncoding = "")

# check the dependency of log2FC on GC content
gcGroups.sub=gcGroups[keep]
dfEdgeR = data.frame(logFC=log(2^lrt.EDASeq$table$logFC), gc=gcGroups.sub)

pedgeR = ggplot(dfEdgeR) +
  aes(x=gc, y=logFC, color=gc) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  scale_color_manual(values=wesanderson::wes_palette("Zissou1", nlevels(gcGroups), "continuous")) +
  geom_abline(intercept = 0, slope = 0, col="black", lty=2) +
  ylim(c(-1,1)) +
  ggtitle("log2FCs in bins by GC content, FQ-FQ normalisation") +
  xlab("GC-content bin") +
  theme_bw()+
  theme(aspect.ratio = 1)+
  theme(axis.text.x = element_text(angle = 45, vjust = .5),
        legend.position = "none",
        axis.title = element_text(size=16))

ggsave(filename="log2FC_vs_GCcontent.EDAseq.pdf",plot=pedgeR ,path=".",device="pdf")

# Peaks Coverage Plot
pth2peaks_bed="nk_merged_peaksid.bed"

peaks.bed=read.table(pth2peaks_bed, sep="\t", header=FALSE, blank.lines.skip=TRUE)
rownames(peaks.bed)=peaks.bed[,4]

peaks.gr <- GRanges(seqnames=peaks.bed[,1], ranges=IRanges(peaks.bed[,2], peaks.bed[,3]), strand="*", mcols=data.frame(peakID=peaks.bed[,4]))

pdf("PeakCoverage.pdf")
covplot(peaks.gr, chrs=c("chr14", "chr15"))
dev.off()

# annotate peaks with closest genomic features
bed.annot = annotatePeak(peaks.gr, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
bed.annot
annot_peaks=as.data.frame(bed.annot)
head(annot_peaks)

write.table(annot_peaks, "nk_merged_annotated.txt",
        append = FALSE,
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        fileEncoding = "")

pdf("AnnotVis.pdf")
upsetplot(bed.annot, vennpie=TRUE)
dev.off()

# distribution of loci with respect to TSS
pdf("TSSdist.pdf")
plotDistToTSS(bed.annot, title="Distribution of ATAC-seq peaks loci\nrelative to TSS")
dev.off()

# finding enriched Reactome pathways using chromosome 1 and 2 genes as a background
pathway.reac <- enrichPathway(as.data.frame(annot_peaks)$geneId)

# previewing enriched Reactome pathways
head(pathway.reac)

colnames(as.data.frame(pathway.reac))
pathway.reac[1:10,c(1:7,9)]

# enriched GO terms
pathway.GO <- enrichGO(as.data.frame(annot_peaks)$geneId, org.Hs.eg.db, ont = "MF")
head(pathway.GO)











