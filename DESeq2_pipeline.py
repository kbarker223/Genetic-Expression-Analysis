from pyspark.sql import SparkSession
import pandas as pd
import os

# ── Config ────────────────────────────────────────────────────────────────────
COUNTS_PATH  = "gs://gene_datasets/cpm_matrix_raw.xlsx"
PCA_PATH     = "gs://gene_datasets/pca_summary.csv"
OUTPUT_PATH  = "gs://gene_datasets/deseq2_results.parquet"

# ── Spark ─────────────────────────────────────────────────────────────────────
os.environ['PYSPARK_PYTHON'] = '/opt/conda/miniconda3/bin/python3'

spark = SparkSession.builder \
    .appName("DESeq2-GProfiler") \
    .config("spark.driver.memory", "20g") \
    .getOrCreate()
sc = spark.sparkContext

# ── Load and prepare counts ───────────────────────────────────────────────────
print("Reading data...")
counts = pd.read_excel(COUNTS_PATH)

target_gene_id = counts[counts['gene_name'] == 'SNHG14']['gene_id'].iloc[0]

transposed = (
    counts
    .set_index('gene_id')
    .rename(index={target_gene_id: 'Target(SNHG14)'},
            columns={'hSR6_1_cpm': 'Control_1(HSR6)', 'hSR6_2_cpm': 'Control_2(HSR6)'})
    .reset_index()
    .drop(columns=['gene_name', 'gene_biotype'])
    .set_index('gene_id')
    .T
)
transposed.index.name = 'sample'

t_counts = (
    transposed[transposed.index.str.endswith('_count')]
    .copy()
    .round()
    .astype(int)
)
t_counts = t_counts.loc[:, t_counts.sum(axis=0) > 0]
print(f"Genes retained: {t_counts.shape[1]}")

# ── Metadata ──────────────────────────────────────────────────────────────────
def sample_to_condition(sample):
    group = sample.replace('_count', '').rsplit('_', 1)[0]
    return 'control' if group == 'hSR6' else group

metadata = pd.DataFrame(
    {"condition": [sample_to_condition(s) for s in t_counts.index]},
    index=t_counts.index
)

# ── Pareto-optimal treatments ─────────────────────────────────────────────────
all_treatments = [c for c in metadata['condition'].unique() if c != 'control']
pca_summary    = pd.read_csv(PCA_PATH)
pareto         = pca_summary[pca_summary['pareto_optimal']]['sample_type'].tolist()
treatments     = [t for t in pareto if t in all_treatments]
print(f"{len(treatments)} pareto-optimal treatments: {treatments}")

# ── Broadcast to workers ──────────────────────────────────────────────────────
t_counts_bc = sc.broadcast(t_counts)
metadata_bc  = sc.broadcast(metadata)

# ── DESeq2 per treatment ──────────────────────────────────────────────────────
def run_treatment(trt):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    tc   = t_counts_bc.value
    md   = metadata_bc.value
    mask = md["condition"].isin(["control", trt])

    dds = DeseqDataSet(counts=tc[mask], metadata=md[mask], design_factors="condition")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", trt, "control"])
    stat_res.summary()

    results = stat_res.results_df
    results = results[results['baseMean'] > 0].copy()
    results["treatment"] = trt

    print(f"Finished {trt}: {(results['padj'] < 0.05).sum()} DEGs")
    return results.reset_index().to_dict("records")

outputs = sc.parallelize(treatments, len(treatments)).map(run_treatment).collect()

# ── Save results ──────────────────────────────────────────────────────────────
combined = pd.concat([pd.DataFrame(o) for o in outputs])
combined.to_parquet(OUTPUT_PATH)

print(f"DESeq2: {len(combined)} rows → {OUTPUT_PATH}")
print(combined[combined['padj'] < 0.05].groupby('treatment').size())
