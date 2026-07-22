# Phase 1: Quality Control (QC) & Preprocessing
library(ggplot2)

# Load the raw data saved from Phase 1
load("raw_luad_data.RData")

print("--- Data Summary ---")
print(paste("Probes:", nrow(data_matrix)))
print(paste("Samples:", ncol(data_matrix)))

# 1. Check for missing values (NA)
na_count <- sum(is.na(data_matrix))
print(paste("Total missing values (NA):", na_count))

# If NAs exist, impute or drop them (microarrays usually shouldn't have NAs)
if (na_count > 0) {
  data_matrix <- na.omit(data_matrix)
}

# 2. Check if data is log2 transformed
max_val <- max(data_matrix, na.rm = TRUE)
print(paste("Max expression value:", max_val))

if (max_val > 100) {
  print("Data appears to be un-logged. Applying log2 transformation...")
  data_matrix <- log2(data_matrix + 1)
} else {
  print("Data is already log2 transformed.")
}

# 3. Principal Component Analysis (PCA) for Outlier Detection
pca <- prcomp(t(data_matrix), scale. = TRUE)
pca_data <- as.data.frame(pca$x)

# Plot PCA to visually inspect sample clustering and outliers
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "steelblue", size = 3) +
  theme_minimal() +
  labs(
    title = "PCA Plot - LUAD Discovery Dataset (GSE115002)",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% Variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% Variance)")
  )

# Save cleaned matrix and metadata
save(data_matrix, metadata, file = "clean_luad_data.RData")
print("QC Complete! Cleaned dataset saved as clean_luad_data.RData")