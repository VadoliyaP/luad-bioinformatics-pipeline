## 🧬 Understanding Lung Adenocarcinoma (LUAD) Bioinformatics Pipeline

**Lung Adenocarcinoma (LUAD)** is the most prevalent histological subtype of non-small cell lung cancer (NSCLC). This pipeline provides an automated systems biology workflow designed to bridge raw genomics data with clinical translation, utilizing both R and Python to process, analyze, and visualize multi-omic profiles.



### 🔄 Multi-Phase Workflow Architecture

The project systematically moves from raw data extraction to interactive clinical deployment across five distinct phases:


[GEO Discovery Dataset] ➔ [DESeq2 / Limma DEG Analysis] ➔ [WGCNA Network Clustering]
                                                                     │
[Streamlit Interactive App] ◄── [DGIdb / Open Targets Annotation] ◄── [Hub Gene / kME Extraction]
             ▲
             │
[TCGA-LUAD Clinical Cross-Validation (Z-summary Preservation)]

* **Phase 1: Data Acquisition & Quality Control (QC)** Retrieval of high-throughput gene expression datasets from public repositories like NCBI GEO and TCGA, followed by sample metadata cleaning and outlier removal.
* **Phase 2: Differential Expression Analysis (DEGs)** Comparing tumor tissues against normal adjacent tissues to identify significantly upregulated or downregulated genes ($\log_2\text{FC}$, $p_{\text{adj}} < 0.05$).
* **Phase 3: Weighted Gene Co-Expression Network Analysis (WGCNA)** Constructing scale-free biological networks ($\beta = 14, R^2 \ge 0.80$) to cluster co-expressed genes into color-coded modules (e.g., turquoise, blue, brown) and extract high-centrality hub genes using Module Membership ($kME$).
* **Phase 4: Target Druggability & Annotation** Querying pharmacological resources—specifically DGIdb and the Open Targets GraphQL API—to map discovered gene targets to existing drugs, clinical trials, and actionable therapeutic families.
* **Phase 5: Clinical Cross-Validation & Dashboard** Validating the robustness of co-expression modules against independent patient cohorts using module preservation statistics ($Z_{\text{summary}}$) and deploying an interactive web application.

### 🛠️ Technology Stack

* **Statistical Computing & Bioinformatics (R 4.3+):** Utilizes packages such as `DESeq2` and `limma` for differential expression, `WGCNA` for network construction, and `GEOquery` / `TCGAbiolinks` for data retrieval.
* **API Integration & Dashboard (Python 3.10+):** Utilizes `Pandas` for data management, `Requests` for querying external bioinformatics APIs, and `Streamlit` alongside `Plotly` for real-time frontend visualization.
