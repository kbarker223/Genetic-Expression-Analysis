# Genetic-Expression-Analysis

## Files

**generate_sample_by_gene.py** — Preprocesses the raw Excel CPM matrix (`CPM-matrix_Tony.xlsx`) by dropping count columns, and transposing so samples are rows and genes (by `gene_id`) are columns. Also renames control samples as `Control_#(HSR6)` with `#` being the replicate number and specifies SNHG14 as `Target(SNHG14)` to mark the gene of interest we want to lower. Finally, the data is saved as `sample_by_gene.parquet`.

**PCA_Analysis.ipynb** — Runs PCA on standardized gene expression values using `sample_by_gene.parquet`, and produces plots of the first two principal components to measure gene variance between treatments and control. Computes log2 fold change of the target gene (SNHG14) per treatment to rank effectiveness, and measures Euclidean distance from each treatment to the Control centroid in PCA space as a proxy for safety (minimum variation from other genes). Combines both metrics into a Pareto optimality scatter plot to identify treatments that best balance SNHG14 knockdown with minimal off-target disruption.

**opentargets_analysis.ipynb** - Queries the Open Targets public BigQuery platform to retrieve gene-disease association scores for all genes present in the expression matrix. Joins the `association_overall_direct` and `disease` tables, filtering to associations with a score above 0.1, and resolves therapeutic area IDs to human-readable names. Saves the resulting gene-disease association table as `disease_associations.parquet` for later analysis. Visualizes the therapeutic areas across the Pareto-optimal treatments.


**NSC_treatments.ipynb** - Runs Nearest Shrunken Centroid (NSC) bootstrap stability analysis on the Pareto-optimal treatments identified in `PCA_Analysis.ipynb`, using Apache Spark on a Google Cloud Dataproc cluster to parallelize 1,800 model fits across bootstrap replicates, permutation nulls, and shrinkage thresholds. Identifies three tiers of stably selected genes and characterizes their disease associations using `disease_associations.parquet`. Computes per-treatment off-target burden scores weighted by Open Targets association evidence, and decomposes burden into shared (NSC-selected) and treatment-specific components to identify the safest Pareto-optimal candidate. Visualizations and tables included for each step of the process.

