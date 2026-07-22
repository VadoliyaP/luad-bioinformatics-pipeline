# Phase 3: Co-Expression Network & Module Detection (WGCNA)
library(WGCNA)
options(stringsAsFactors = FALSE)

# Enable multithreading if supported
enableWGCNAThreads()

# 1. Load Clean Data and DEA Results
load("clean_luad_data.RData")
load("dea_results.RData")

# Filter Top DEGs for WGCNA (Top 4000 by Adjusted P-value)
top_deg_genes <- head(rownames(sig_genes[order(sig_genes$adj.P.Val), ]), 4000)
datExpr <- t(data_matrix[top_deg_genes, ])

# 2. Pick Soft-Thresholding Power (beta)
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot Soft Thresholding Results
pdf("wgcna_soft_threshold.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlabel="Soft Threshold (power)", ylabel="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.80, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlabel="Soft Threshold (power)", ylabel="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

# Automatically choose power where R^2 >= 0.8, default to 6 if none
chosen_power <- sft$powerEstimate
if (is.na(chosen_power)) chosen_power <- 6
cat("Chosen Soft-Thresholding Power (beta):", chosen_power, "\n")

# 3. Construct Co-Expression Network and Detect Modules
net <- blockwiseModules(
  datExpr, 
  power = chosen_power,
  TOMType = "unsigned", 
  minModuleSize = 30,
  remergeDefault = TRUE,
  numericLabels = TRUE, 
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# 4. Save WGCNA Results
moduleColors <- labels2colors(net$colors)
save(datExpr, net, moduleColors, chosen_power, file = "wgcna_results.RData")

cat("\nPhase 3 Network Construction Complete! Saved to 'wgcna_results.RData'.\n")