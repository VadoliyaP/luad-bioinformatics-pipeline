# Phase 2: Differential Expression Analysis (limma)
library(limma)

# 1. Load cleaned data
load("clean_luad_data.RData")

# 2. Extract and format sample groups from source_name_ch1
groups <- metadata$source_name_ch1

# Convert to clean factor names ("Normal" and "Tumor")
group_factor <- factor(ifelse(grepl("Tumor", groups), "Tumor", "Normal"), 
                       levels = c("Normal", "Tumor"))

print("Sample Group Summary:")
print(table(group_factor))

# 3. Create Design Matrix
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

# 4. Fit Linear Model
fit <- lmFit(data_matrix, design)

# Create contrast: Tumor vs Normal
contrast_matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 5. Extract Ranked Results
top_genes <- topTable(fit2, number = Inf, adjust.method = "BH")

# Filter significant Differentially Expressed Genes (DEGs)
# Adjusted P-value < 0.05 and |log2 Fold Change| > 1
sig_genes <- subset(top_genes, adj.P.Val < 0.05 & abs(logFC) > 1)

# 6. Print Summary Statistics
cat("\n--- Phase 2: Differential Expression Summary ---\n")
cat("Total genes tested:", nrow(top_genes), "\n")
cat("Significant DEGs (adj P < 0.05, |logFC| > 1):", nrow(sig_genes), "\n")
cat("  - Upregulated in Tumor:", sum(sig_genes$logFC > 1), "\n")
cat("  - Downregulated in Tumor:", sum(sig_genes$logFC < -1), "\n")

# 7. Save Outputs
write.csv(sig_genes, "significant_deg_luad.csv")
save(top_genes, sig_genes, file = "dea_results.RData")

cat("\nPhase 2 Complete! Results saved to 'significant_deg_luad.csv' and 'dea_results.RData'.\n")