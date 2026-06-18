# Genetic Expression Analysis
### Segal Lab Collaboration: UC Davis

A computational genomics study conducted in collaboration with the **Segal Lab at UC Davis**, investigating the off-target safety profile of engineered gene regulators designed to suppress **SNHG14**, a long non-coding RNA implicated in Angelman syndorme. The project integrates RNA-seq expression data, dimensionality reduction, bootstrapped machine learning, and publicly available disease association databases to rank candidate treatments by their SNHG14 knockdown efficacy relative to off-target transcriptional disruption.

---

## Research Context

The Segal Lab develops engineered transcription factors, including **Artificial Transcription Factors (hATFs)** and **Zinc Finger proteins (nZFs)**, capable of targeting and silencing specific genes. This study evaluates a panel of 22 such regulators in a human cell, all designed to suppress SNHG14. The central question is: *which candidates most effectively reduce SNHG14 expression while minimizing off-target effects?*

To answer this, the pipeline characterizes each treatment's expression signature across ~79,000 transcripts (CPM-normalized RNA-seq), identifies a Pareto-optimal frontier balancing strength against off-target disruptions, and then performs a multi-layer risk assessment of those selected treatments using clinical disease association evidence from Open Targets.

---

## Data Sources

| Source | Description |
|---|---|
| **CPM Expression Matrix** (`CPM-matrix_Tony.xlsx`) | CPM-normalized RNA-seq counts for 48 samples × ~79,000 genes. 22 treatment conditions (2 replicates each) plus 2 HSR6 control replicates. Produced by Tony. |
| **Open Targets Platform** | Public BigQuery dataset (`open-targets-prod.platform`) providing gene–disease association scores across 705,000+ associations. Queried via the GCP project `gene-expression-big-data`. |
| **Ensembl Gene Annotations** | Gene IDs (ENSG…) and biotype classifications embedded in the expression matrix. |

---

## Repository Structure & File Attribution

### Data Preparation & ETL

**`load_data.py`** *(Kai Barker)*
ETL module with two primary functions. `raw_data()` performs a one-time preprocessing step on the original Excel CPM matrix: drops raw count columns, renames the target gene to `Target(SNHG14)`, and renames control columns to `Control_#(HSR6)`, then saves the result as `data/gene_expression_matrix.parquet`. `pipeline_pull_opentargets()` transposes the matrix into sample-by-gene format (`data/sample_by_gene.parquet`), then queries the Open Targets BigQuery platform for disease associations across all expressed genes, resolves nested therapeutic area IDs to human-readable names, and saves the result as `data/disease_associations.parquet`.

**`generate_sample_by_gene.py`** *(Tony)*
Companion preprocessing script that produces the `sample_by_gene.parquet` format used as input across all downstream analyses.

---

### Exploratory Analysis

**`EDA.ipynb`**
Initial exploration of the expression matrix: sample-level summaries, gene biotype distributions, and data shape verification across the 22 treatment conditions and ~79,000 gene features.

---

### Dimensionality Reduction & Treatment Ranking

**`PCA_Analysis.ipynb`** *(Tony)*
Runs PCA on standardized gene expression values and produces biplots of the first two principal components to visualize separation between treatments and control. Computes log2 fold change of SNHG14 per treatment to rank effectiveness, and Euclidean distance from each treatment centroid to the control centroid in PCA space as a proxy for transcriptome-wide safety. Combines both metrics into a **Pareto optimality scatter plot** to identify the four candidates: `hATF567`, `hATF561`, `nZF105`, `nZF139` that best balance target knockdown with minimal off-target disruption. Outputs `data/pca_summary.csv` used by all downstream analyses.

**`multi_dim_PCA.ipynb`**
Extended PCA exploration across additional principal components.

---

### Disease Association Analysis

**`opentargets_analysis.ipynb`** *(Kai)*
Loads the disease association table from `data/disease_associations.parquet` and characterizes the therapeutic area footprint of each Pareto-optimal treatment. For each candidate, identifies genes with |log2FC| >= 1.0 relative to control, joins them to Open Targets associations, and aggregates disease impact by therapeutic area. Produces an interactive heatmap (raw counts and normalized %) of off-target therapeutic area disruption across the four finalists. Each treatment disturbs 25–64 genes at this threshold, touching 20–24 distinct therapeutic areas.

---

### Bootstrapped Classifier & Risk Scoring

**`NSC_risk_analysis.ipynb`** *(Kai)*
Core risk characterization notebook, executed on a Google Cloud Dataproc cluster. Runs **Nearest Shrunken Centroid (NSC) bootstrap stability analysis** across 1,800 Spark-parallelized model fits (50 bootstrap replicates across 6 shrinkage thresholds) to identify genes that reliably distinguish Pareto-optimal treatments from control. A permutation null (100 shuffles per threshold) validates statistical significance. Threshold refinement between 2.0 and 2.5 at 0.1 increments identifies the most stringent statistically significant cutoff (p < 0.001).

Three gene tiers are defined by stability across thresholds:
- **7-gene most-robust set**: 100% bootstrap stability across all thresholds <= 2.5
- **22-gene primary set**: >= 50% stability at threshold 2.4
- **1,484-gene broad signature**: 100% stability at threshold 0.5

Disease associations are mapped onto each tier and a **composite off-target risk score** is computed per treatment (disturbed gene count, evidence-weighted burden, high-confidence disruptions, therapeutic areas touched). Burden is then decomposed into shared (NSC-selected) vs. treatment-specific components to identify the safest Pareto-optimal candidate. `hATF561` emerges as the lowest-risk treatment with 25 disturbed genes and the smallest evidence-weighted burden.

---

### Differential Expression

**`DESeq2_analysis.ipynb`**
Runs **DESeq2** (via `pydeseq2`) on raw count data for the four Pareto-optimal treatments using Spark to parallelize across treatments. Identifies statistically significant differentially expressed genes (padj < 0.05): `hATF561` (7 genes), `nZF139` (31 genes), `hATF567` (35 genes), `nZF105` (45 genes), reinforcing the risk ordering from the NSC analysis.

**`DESeq2_pipeline.py`**
Supporting script for the DESeq2 workflow.

---

### Functional Enrichment

**`DAVID_process_extraction.ipynb`**
Functional enrichment analysis using the DAVID bioinformatics tool to characterize biological processes over-represented among differentially expressed genes.

---

### Network & Co-expression Analyses *(Lucas)*

**`Lucas_PPI/`**
Protein–protein interaction network analysis on gene sets identified in the expression study.

**`Lucas_WGCNA/`**
Weighted Gene Co-expression Network Analysis (WGCNA) to identify co-expressed gene modules and their relationship to treatment conditions.

---

## Cloud Infrastructure

All computationally intensive analyses were executed on **Google Cloud Platform** using the following services:

| Service | Role |
|---|---|
| **GCP Project** | `gene-expression-big-data` |
| **Google Cloud Dataproc** | Managed Spark cluster (`mycluster`, region `us-central1`, `e2-highmem-4` master node) used to parallelize NSC bootstrap fits (1,800 jobs) and DESeq2 runs across Pareto-optimal treatments |
| **Apache Spark (PySpark)** | Distributed computation framework running on Dataproc; workloads broadcast via `SparkContext` and parallelized with `sc.parallelize()` |
| **Google Cloud Storage** | GCS bucket `gene_datasets` serving as the data lake for all intermediate and final parquet/CSV outputs shared across cluster sessions |
| **BigQuery** | Used to query the Open Targets public dataset (`open-targets-prod.platform`) for gene–disease association scores; parameterized queries issued via the `google-cloud-bigquery` Python client |
| **JupyterLab on Dataproc** | Cluster notebooks accessed via SSH tunnel from local VS Code; cluster configured with auto-deletion after 2 hours idle or 5 hours total |

---


## Contributors

| Contributor | Affiliation | Contributions |
|---|---|---|
| **Kai Barker** | UC Santa Barbara | `load_data.py`, `opentargets_analysis.ipynb`, `NSC_risk_analysis.ipynb` |
| **Tony Segal** | UC Santa Barbara | `load_data.py`, `generate_sample_by_gene.py`, `PCA_Analysis.ipynb` |
| **Lucas Childs** | UC Santa Barbara | `Lucas_PPI/`, `Lucas_WGCNA/`
