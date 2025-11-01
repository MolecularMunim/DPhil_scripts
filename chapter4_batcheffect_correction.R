# --- Load required libraries ---
library(limma)
library(ggplot2)

# --- Step 1: Load metadata ---
sample_table_path <- "PCA/Brain_transcriptome/sample_table.csv"
sampleTable <- read.csv(sample_table_path, row.names = 1)

# --- Step 2: Build count file paths ---
count_dir <- "PCA/Brain_transcriptome/downsampled_brain_counts/"
countFiles <- file.path(count_dir, paste0(rownames(sampleTable), "fb.txt"))

# --- Step 3: Read count data ---
countList <- lapply(countFiles, function(f) read.table(f, header = TRUE, row.names = 1))
countMatrix <- do.call(cbind, lapply(countList, function(x) x[, 1]))  # Assuming counts are in the 1st column
colnames(countMatrix) <- rownames(sampleTable)

# --- Step 4: Define experimental variables ---

batch <- sampleTable$Dataset
condition <- sampleTable$condition

# --- Step 5: Log-transform counts ---

logCounts <- log2(countMatrix + 1)

# --- Step 6: Remove batch effects ---

# Build design matrix to preserve biological condition
design <- model.matrix(~ condition)
# Remove batch effects while keeping condition effects intact
logCounts_corrected <- removeBatchEffect(logCounts, batch = batch, design = design)

# ============================================================
# PCA Analysis
# ============================================================

# --- Helper function for PCA computation ---
perform_pca <- function(mat) {
  var_genes <- apply(mat, 1, var)
  mat_filtered <- mat[var_genes > 0, ]
  prcomp(t(mat_filtered), scale. = TRUE)
}

# Perform PCA on uncorrected and corrected data
pca_uncorrected <- perform_pca(logCounts)
pca_corrected   <- perform_pca(logCounts_corrected)

# --- Prepare PCA dataframes ---

pcaData_uncorrected <- as.data.frame(pca_uncorrected$x)
pcaData_uncorrected$condition <- condition
pcaData_uncorrected$Dataset <- batch

pcaData_corrected <- as.data.frame(pca_corrected$x)
pcaData_corrected$condition <- condition
pcaData_corrected$Dataset <- batch

# --- Helper to extract explained variance ---
expl_var <- function(pca_obj, pc_num) {
  round(summary(pca_obj)$importance[2, pc_num] * 100, 1)
}

# ============================================================
# Visualization
# ============================================================

# Common ggplot theme
square_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = c(0.80, 0.20),
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20, face = "bold"),
    aspect.ratio = 1
  )

# --- Plot function ---
plot_pca <- function(data, pca_obj, color_by, title_suffix) {
  ggplot(data, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 3) +
    coord_equal() +
    square_theme +
    labs(
      title = paste0("PCA (", title_suffix, ")"),
      x = paste0("PC1 (", expl_var(pca_obj, 1), "%)"),
      y = paste0("PC2 (", expl_var(pca_obj, 2), "%)")
    )
}

# --- Generate plots ---
p1 <- plot_pca(pcaData_uncorrected, pca_uncorrected, "condition", "Colored by Tissue, Uncorrected")
p2 <- plot_pca(pcaData_uncorrected, pca_uncorrected, "Dataset", "Colored by Dataset, Uncorrected")
p3 <- plot_pca(pcaData_corrected, pca_corrected, "condition", "Colored by Tissue, Batch-Corrected") +
  theme(legend.position = c(0.20, 0.80))  # top-left area

p4 <- plot_pca(pcaData_corrected, pca_corrected, "Dataset", "Colored by Dataset, Batch-Corrected") +
  theme(legend.position = c(0.20, 0.80))

# --- Display plots ---
p1
p2
p3
p4

# Save all PCA plots ----
ggsave("p1_tissue_uncorrected.png", p1, width = 10, height = 10, dpi = 300)
ggsave("p2_dataset_uncorrected.png", p2, width = 10, height = 10, dpi = 300)
ggsave("p3_tissue_corrected.png", p3, width = 10, height = 10, dpi = 300)
ggsave("p4_dataset_corrected.png", p4, width = 10, height = 10, dpi = 300)
