# Phase 1: Data Acquisition
# Target: Lung Adenocarcinoma (LUAD) Discovery Dataset
library(GEOquery)
library(DESeq2)

# Define the GEO accession
gse_id <- "GSE115002" 

print(paste("Starting download for:", gse_id))

# Fetch the dataset from GEO
gse <- getGEO(gse_id, destdir = ".", getGPL = TRUE)

# Check how many elements were returned (usually length is 1)
print(paste("Number of platforms/series in this GEO accession:", length(gse)))

# Extract expression matrix and metadata from the FIRST element [[1]]
gse_obj <- gse[[1]]

data_matrix <- exprs(gse_obj)
metadata <- pData(gse_obj)

# Inspect the dimensions to ensure it loaded properly
print(paste("Matrix dimensions (Genes x Samples):", paste(dim(data_matrix), collapse = " x ")))

# Save the raw objects for Quality Control
save(data_matrix, metadata, file = "raw_luad_data.RData")

print("Phase 1: Raw data successfully acquired and saved.")