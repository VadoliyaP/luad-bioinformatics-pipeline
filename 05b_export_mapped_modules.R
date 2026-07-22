# ==============================================================================
# Step 05b: Export WGCNA Gene Modules with Probe-to-Symbol Mapping
# ==============================================================================

cat("Checking and loading required packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("hgug4112a.db", quietly = TRUE)) BiocManager::install("hgug4112a.db", update = FALSE)

library(hgug4112a.db)
library(AnnotationDbi)

# 1. Load WGCNA workspace
rdata_file <- "wgcna_final_results.RData"
if (!file.exists(rdata_file)) {
  stop("Error: 'wgcna_final_results.RData' not found.")
}

load(rdata_file)
cat("Loaded 'wgcna_final_results.RData' successfully.\n")

# 2. Extract Probe IDs and Module Colors
probes <- colnames(datExpr)
if (is.null(probes)) {
  probes <- names(moduleColors)
}

# 3. Map Agilent Probe IDs to Official HGNC Gene Symbols
cat("Mapping microarray probes to HGNC Gene Symbols...\n")
mapped_symbols <- mapIds(
  hgug4112a.db,
  keys = probes,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# 4. Extract kME (Module Membership) values
if (exists("geneModuleMembership") && !is.null(geneModuleMembership)) {
  kme_values <- apply(geneModuleMembership, 1, max, na.rm = TRUE)
} else {
  kme_values <- rep(1.0, length(probes))
}

# 5. Build Final Data Frame
df_out <- data.frame(
  ProbeID = probes,
  GeneSymbol = as.character(mapped_symbols),
  ModuleColor = moduleColors,
  kME = kme_values,
  stringsAsFactors = FALSE
)

# CRITICAL FIX: Keep ONLY probes that mapped to real HGNC Gene Symbols
df_out <- df_out[!is.na(df_out$GeneSymbol) & df_out$GeneSymbol != "" & !grepl("^A_", df_out$GeneSymbol), ]

# Deduplicate if multiple probes map to the same Gene Symbol (keep highest kME)
df_out <- df_out[order(df_out$GeneSymbol, -df_out$kME), ]
df_out <- df_out[!duplicated(df_out$GeneSymbol), ]

# Export to CSV
output_csv <- "wgcna_gene_modules.csv"
write.csv(df_out, file = output_csv, row.names = FALSE)

cat("\n=========================================================\n")
cat(" SUCCESS: Created 'wgcna_gene_modules.csv'!\n")
cat(" Total valid gene targets exported:", nrow(df_out), "\n")
cat(" Unique Gene Symbols mapped:", length(unique(df_out$GeneSymbol)), "\n")
cat("=========================================================\n")