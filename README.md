## 🧬 Understanding Lung Adenocarcinoma (LUAD) Bioinformatics Pipeline

This pipeline is a complete computational workflow built to study lung adenocarcinoma (LUAD)—the most common form of lung cancer. Instead of just looking at single genes one by one, it looks at the big picture of how genes interact with each other inside tumor cells.

By combining R for statistical modeling and Python for connecting with online medical databases, the pipeline takes raw genetic data from real cancer patients, cleans it up, groups related genes into networks to find the most important "hub" drivers, and checks whether those targets can actually be treated with existing drugs. Finally, everything comes together in an interactive web dashboard so you can easily search for genes and explore the results visually.



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
**Statistical Computing & Bioinformatics (R 4.3+):**

* DESeq2 / limma — Differential gene expression.

* WGCNA / fastcluster / dynamicTreeCut — Co-expression network construction and module detection.

 * GEOquery & TCGAbiolinks — Data retrieval from GEO and TCGA repositories.

 * AnnotationDbi & platform-specific annotation packages (hgug4112a.db) — Probe-to-gene symbol mapping.

**API Integration, Networks & Dashboard (Python 3.10+):**

 * Pandas / NumPy — Data manipulation and table preprocessing.

 * Requests / GraphQL — Programmatic fetching from DGIdb and Open Targets APIs.

* Streamlit — Interactive frontend dashboard design.

 * Plotly / Matplotlib / Seaborn — Network and scatter visualizations.
