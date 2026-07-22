# Phase 3 (Refined): WGCNA Module Analysis & Hub Gene Identification
library(WGCNA)
options(stringsAsFactors = FALSE)

# 1. Load Data and DEA Results
load("clean_luad_data.RData")
load("dea_results.RData")

# Top 4,000 DEGs
top_deg_genes <- head(rownames(sig_genes[order(sig_genes$adj.P.Val), ]), 4000)
datExpr <- t(data_matrix[top_deg_genes, ])

# 2. Re-construct Network with Proper Beta Power (14)
chosen_power <- 14
cat("Using Optimized Soft-Thresholding Power (beta):", chosen_power, "\n")

net <- blockwiseModules(
  datExpr, 
  power = 14,
  TOMType = "unsigned", 
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  remergeDefault = TRUE,
  numericLabels = FALSE, 
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# Extract module colors and Module Eigengenes (MEs)
moduleColors <- net$colors
MEs <- net$MEs

# Print detected modules
cat("\nModule Sizes:\n")
print(table(moduleColors))

# 3. Module-Trait Correlation
# Binary trait encoding: Normal = 0, Tumor = 1
trait_data <- data.frame(
  Tumor_Status = ifelse(grepl("Tumor", metadata$source_name_ch1), 1, 0)
)
rownames(trait_data) <- rownames(metadata)

# Calculate correlations between MEs and Tumor Status
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Display Module-Trait Correlations
trait_summary <- data.frame(
  Module = colnames(MEs),
  Correlation = moduleTraitCor[, 1],
  PValue = moduleTraitPvalue[, 1]
)
trait_summary <- trait_summary[order(-abs(trait_summary$Correlation)), ]

cat("\n--- Top Modules Correlated with Tumor Status ---\n")
print(head(trait_summary, 10))

# 4. Hub Gene Extraction for Top Key Module
top_module <- gsub("ME", "", trait_summary$Module[1])
cat("\nTop Module Identified:", top_module, "\n")

# Calculate Module Membership (kME)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
top_module_kME <- geneModuleMembership[, paste0("ME", top_module)]
names(top_module_kME) <- colnames(datExpr)

# Extract Top 20 Hub Genes for the module
hub_genes <- head(sort(abs(top_module_kME), decreasing = TRUE), 20)

cat("\nTop 10 Hub Genes in Module '", top_module, "':\n", sep = "")
print(head(hub_genes, 10))

# 5. Save Final WGCNA Results
save(datExpr, net, MEs, moduleColors, trait_summary, hub_genes, file = "wgcna_final_results.RData")
write.csv(trait_summary, "wgcna_module_trait_correlations.csv", row.names = FALSE)

cat("\nPhase 3 Analysis Complete! Results saved to 'wgcna_module_trait_correlations.csv' and 'wgcna_final_results.RData'.\n")