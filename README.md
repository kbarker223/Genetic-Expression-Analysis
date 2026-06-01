# Genetic-Expression-Analysis

## Files

**generate_sample_by_gene.py** — Preprocesses the raw Excel CPM matrix (`CPM-matrix_Tony.xlsx`) by dropping count columns, and transposing so samples are rows and genes (by `gene_id`) are columns. Also renames control samples as `Control_#(HSR6)` with `#` being the replicate number and specifies SNHG14 as `Target(SNHG14)` to mark the gene of interest we want to lower. Finally, the data is saved as `sample_by_gene.parquet`.

**PCA_Analysis.ipynb** — Runs PCA on standardized gene expression values using `sample_by_gene.parquet`, and produces plots of the first two principal components to measure gene variance between treatments and control. Computes log2 fold change of the target gene (SNHG14) per treatment to rank effectiveness, and measures Euclidean distance from each treatment to the Control centroid in PCA space as a proxy for safety (minimum variation from other genes). Combines both metrics into a Pareto optimality scatter plot to identify treatments that best balance SNHG14 knockdown with minimal off-target disruption.

**NCS_treatments.ipynb** - 