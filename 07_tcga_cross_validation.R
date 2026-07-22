# =========================================================
# Phase 4b: TCGA-LUAD Module Cross-Validation
# =========================================================
library(WGCNA)
library(TCGAbiolinks)
library(DESeq2)
library(dplyr)

enableWGCNAThreads()

# 1. Load Discovery Workspace (GEO Data)
load("wgcna_final_results.RData") 
# Expecting: datExpr (GEO expression matrix: samples x genes), moduleColors

cat("Downloading TCGA-LUAD Validation Cohort...\n")
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, method = "api", files.per.chunk = 10)
tcga_data <- GDCprepare(query)

# Extract counts and convert to TPM / Log2 Normalized
tcga_counts <- assay(tcga_data, "unstranded")
# Convert Ensembl IDs to Gene Symbols
gene_info <- rowData(tcga_data)
rownames(tcga_counts) <- gene_info$gene_name

# Normalize TCGA Data
tcga_dds <- DESeqDataSetFromMatrix(countData = tcga_counts, 
                                  colData = colData(tcga_data), 
                                  design = ~ 1)
tcga_vsd <- vst(tcga_dds, blind = TRUE)
datExpr_tcga <- t(assay(tcga_vsd))

# Find overlapping genes between GEO and TCGA datasets
common_genes <- intersect(colnames(datExpr), colnames(datExpr_tcga))
cat(sprintf("Common genes across GEO and TCGA: %d\n", length(common_genes)))

datExpr_ref  <- datExpr[, common_genes]
datExpr_test <- datExpr_tcga[, common_genes]
color_ref    <- moduleColors[colnames(datExpr) %in% common_genes]

# 2. Run WGCNA Module Preservation Test
set.seed(123)
multiExpr <- list(GEO = list(data = datExpr_ref), TCGA = list(data = datExpr_test))
multiColor <- list(GEO = color_ref)

cat("Calculating Module Preservation Statistics (Permutations = 100)...\n")
mp <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  nPermutations = 100,
  randomSeed = 123,
  verbose = 3
)

# Extract Zsummary scores
ref <- 1
test <- 2
statsObs <- mp$preservation$observed[[ref]][[test]]
statsZ   <- mp$preservation$Z[[ref]][[test]]

preservation_df <- data.frame(
  Module = rownames(statsZ),
  ModuleSize = statsObs$moduleSize,
  Zsummary = statsZ$Zsummary.pres
)

# Filter for key modules
key_modules <- preservation_df %>% filter(Module %in% c("turquoise", "blue", "brown"))
print("=========================================================")
print("TCGA Cross-Validation Preservation Results:")
print(key_modules)
print("=========================================================")

write.csv(preservation_df, "tcga_module_preservation_results.csv", row.names = FALSE)
cat("Validation complete! Results saved to 'tcga_module_preservation_results.csv'.\n")