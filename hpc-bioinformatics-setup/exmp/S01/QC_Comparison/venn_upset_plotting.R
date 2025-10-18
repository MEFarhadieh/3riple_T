# ====================================================================
# R script for creating Venn diagrams and UpSet plot
# 5-way comparison
# ====================================================================

library(VennDiagram)
library(UpSetR)
library(tidyverse)
library(grid)

# ====================================================================
# Set paths
# ====================================================================

results_dir <- "/work/archive/farhadie/public_studies/Hill/QC_comparison_analysis/results"
figures_dir <- "/work/archive/farhadie/public_studies/Hill/QC_comparison_analysis/figures"

cat("Working directories:\n")
cat("  Results:", results_dir, "\n")
cat("  Figures:", figures_dir, "\n\n")

# ====================================================================
# Load data
# ====================================================================

cat("Loading data...\n")

venn_data <- read.csv(file.path(results_dir, "venn_pairwise_data.csv"))
upset_data_raw <- read.csv(file.path(results_dir, "upset_5way_data.csv"), stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("  Venn data rows:", nrow(venn_data), "\n")
cat("  UpSet data rows:", nrow(upset_data_raw), "\n\n")

# ====================================================================
# Create pairwise Venn diagrams (QClus vs each Combined)
# ====================================================================

cat("Creating Venn diagrams...\n\n")

comparisons <- c("Combined_5pct", "Combined_10pct", "Combined_20pct", "Combined_25pct")
labels <- c("5pct", "10pct", "20pct", "25pct")
display_labels <- c("5%", "10%", "20%", "25%")

for(i in 1:length(comparisons)) {
  
  cat("Processing:", comparisons[i], "\n")
  
  row <- venn_data[venn_data$Method2 == comparisons[i], ]
  
  venn.plot <- draw.pairwise.venn(
    area1 = row$Method1_Total,
    area2 = row$Method2_Total,
    cross.area = row$Overlap,
    category = c("QClus", paste0("Combined QC (", display_labels[i], " SoupX)")),
    fill = c("#E74C3C", "#3498DB"),
    alpha = 0.5,
    cat.fontface = "bold",
    cat.cex = 1.3,
    cex = 1.6,
    cat.pos = c(-20, 20),
    cat.dist = 0.05,
    euler.d = TRUE,
    scaled = TRUE
  )
  
  # Save PNG
  png_file <- file.path(figures_dir, paste0("venn_QClus_vs_Combined_", labels[i], ".png"))
  png(png_file, width = 900, height = 900, res = 150)
  grid.draw(venn.plot)
  dev.off()
  cat("  Saved PNG:", png_file, "\n")
  
  # Save PDF
  pdf_file <- file.path(figures_dir, paste0("venn_QClus_vs_Combined_", labels[i], ".pdf"))
  pdf(pdf_file, width = 9, height = 9)
  grid.draw(venn.plot)
  dev.off()
  cat("  Saved PDF:", pdf_file, "\n\n")
}

cat("All Venn diagrams created\n\n")

# ====================================================================
# Create 5-way UpSet plot - FIXED DATA CONVERSION
# ====================================================================

cat("Creating UpSet plot...\n")

# Check data structure
cat("  Checking data structure...\n")
cat("  First few rows of raw data:\n")
print(head(upset_data_raw))
cat("\n  Column types:\n")
print(sapply(upset_data_raw, class))

# Convert string "True"/"False" to logical, then to integer
upset_matrix <- upset_data_raw %>%
  select(QClus, Combined_5pct, Combined_10pct, Combined_20pct, Combined_25pct) %>%
  mutate(
    QClus = as.integer(QClus == "True" | QClus == TRUE | QClus == "TRUE" | QClus == 1),
    Combined_5pct = as.integer(Combined_5pct == "True" | Combined_5pct == TRUE | Combined_5pct == "TRUE" | Combined_5pct == 1),
    Combined_10pct = as.integer(Combined_10pct == "True" | Combined_10pct == TRUE | Combined_10pct == "TRUE" | Combined_10pct == 1),
    Combined_20pct = as.integer(Combined_20pct == "True" | Combined_20pct == TRUE | Combined_20pct == "TRUE" | Combined_20pct == 1),
    Combined_25pct = as.integer(Combined_25pct == "True" | Combined_25pct == TRUE | Combined_25pct == "TRUE" | Combined_25pct == 1)
  )

cat("  After conversion:\n")
print(head(upset_matrix))
cat("\n  Column sums (total flagged per method):\n")
print(colSums(upset_matrix, na.rm = TRUE))

# Check for NAs
if(any(is.na(upset_matrix))) {
  cat("  WARNING: NAs detected, removing rows with NAs\n")
  upset_matrix <- na.omit(upset_matrix)
}

cat("  Final matrix dimensions:", nrow(upset_matrix), "x", ncol(upset_matrix), "\n")

# Create UpSet plot - PNG
png_file <- file.path(figures_dir, "upset_5way_QClus_vs_CombinedQC.png")
png(png_file, width = 1400, height = 900, res = 150)

upset(
  as.data.frame(upset_matrix),
  nsets = 5,
  nintersects = 30,
  order.by = "freq",
  mainbar.y.label = "Number of Cells",
  sets.x.label = "Total Cells Flagged",
  text.scale = 1.3,
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

dev.off()
cat("  Saved PNG:", png_file, "\n")

# Create UpSet plot - PDF
pdf_file <- file.path(figures_dir, "upset_5way_QClus_vs_CombinedQC.pdf")
pdf(pdf_file, width = 14, height = 9)

upset(
  as.data.frame(upset_matrix),
  nsets = 5,
  nintersects = 30,
  order.by = "freq",
  mainbar.y.label = "Number of Cells",
  sets.x.label = "Total Cells Flagged",
  text.scale = 1.3,
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.6, 0.4)
)

dev.off()
cat("  Saved PDF:", pdf_file, "\n\n")

# ====================================================================
# Create alternative intersection bar plot
# ====================================================================

cat("Creating intersection summary bar plot...\n")

# Calculate intersections
intersection_counts <- upset_matrix %>%
  group_by(QClus, Combined_5pct, Combined_10pct, Combined_20pct, Combined_25pct) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

# Create readable labels
intersection_counts$Label <- apply(intersection_counts[,1:5], 1, function(x) {
  sets <- c("QC", "5%", "10%", "20%", "25%")[as.logical(x)]
  if(length(sets) == 0) return("None")
  paste(sets, collapse = " & ")
})

# Show top intersections
cat("\nTop 15 intersections:\n")
print(head(intersection_counts, 15))

# Create bar plot
png_file <- file.path(figures_dir, "intersection_barplot_top20.png")
png(png_file, width = 1600, height = 900, res = 150)

par(mar = c(12, 5, 4, 2))
top20 <- head(intersection_counts, 20)
barplot(
  top20$Count,
  names.arg = top20$Label,
  las = 2,
  col = rainbow(20, alpha = 0.7),
  border = "black",
  main = "Top 20 QC Method Intersections",
  ylab = "Number of Cells",
  cex.names = 0.9,
  cex.axis = 1.2,
  cex.lab = 1.2,
  cex.main = 1.5
)

# Add count labels on bars
text(
  x = 1:nrow(top20) * 1.2 - 0.5,
  y = top20$Count,
  labels = format(top20$Count, big.mark = ","),
  pos = 3,
  cex = 0.8
)

dev.off()
cat("  Saved intersection bar plot:", png_file, "\n")

# Save intersection counts
csv_file <- file.path(results_dir, "intersection_summary.csv")
write.csv(intersection_counts, csv_file, row.names = FALSE)
cat("  Saved intersection summary:", csv_file, "\n\n")

# ====================================================================
# Summary
# ====================================================================

cat("====================================\n")
cat("All visualizations complete!\n")
cat("====================================\n\n")

cat("Created files:\n")
cat("  - 4 Venn diagrams (PNG + PDF each)\n")
cat("  - 1 UpSet plot (PNG + PDF)\n")
cat("  - 1 Intersection bar plot (PNG)\n")
cat("  - 1 Intersection summary (CSV)\n")
cat("  Total: 12 files\n\n")

cat("Files saved in:", figures_dir, "\n")

# List all created visualization files
all_files <- list.files(figures_dir, pattern = "\\.(png|pdf)$", full.names = FALSE)
cat("\nAll visualization files:\n")
for(f in sort(all_files)) {
  cat("  -", f, "\n")
}

cat("\nDone!\n")
