# Load required library
library(ggplot2)

# Example data: TH expression levels across conditions
data <- read.csv(../AD23_gene_loci/ch112_185000_2195000.csv)

# Create bar plot
ggplot(data, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Tyrosine Hydroxylase (TH) Gene Expression",
    x = "Experimental Condition",
    y = "Expression Level (FPKM)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend if redundant

# Bioconductor packages (install if needed)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Gviz")

library(Gviz)

# Define TH gene coordinates (example: human, chr11)
genome <- "hg38"  # or "mm10" for mouse
chr <- "chr11"
start <- 2_185_000  # Update with actual coordinates
end <- 2_195_000

# Create axis track
axis_track <- GenomeAxisTrack()

# Gene region track (simplified)
gene_track <- GeneRegionTrack(
  genome = genome,
  chromosome = chr,
  start = start,
  end = end,
  name = "TH",
  transcriptAnnotation = "symbol"
)

# Plot
plotTracks(list(axis_track, gene_track), 
           main = "TH Gene Locus (chr11:2,185,000-2,195,000)")

# Install if needed: install.packages("pheatmap")
library(pheatmap)

# Example matrix (rows = samples, cols = genes)
expr_matrix <- read.delim(../AD23_expr/expr_sets.tsv)

# Plot
pheatmap(
  expr_matrix[, "TH", drop = FALSE],  # Focus on TH
  main = "TH Expression Heatmap",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = TRUE,
  show_colnames = TRUE
)


# Install: BiocManager::install("karyoploteR")
library(karyoploteR)

# Define TH region (human hg38)
th_region <- toGRanges("chr11:2,185,000-2,195,000")

# Plot
kp <- plotKaryotype(zoom = th_region, genome = "hg38")
kpAddBaseNumbers(kp)
kpPlotGenes(kp, data = th_region, r0 = 0, r1 = 0.3, gene.name.cex = 1.2)
kpAddMainTitle(kp, "TH Gene Region (hg38)", cex = 1.5)
