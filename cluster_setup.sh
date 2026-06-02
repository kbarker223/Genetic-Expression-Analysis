#!/bin/bash
# =============================================================
#  Dataproc Cluster — Connection Cheatsheet
#  Project: gene-expression-big-data
#  Cluster:  mycluster  |  Region: us-central1  |  Zone: us-central1-a
# =============================================================

# -------------------------------------------------------------
# 1. CREATE CLUSTER (run in Google Cloud Shell)
#    Only needed when cluster has been deleted/expired
# -------------------------------------------------------------
gcloud dataproc clusters create mycluster \
    --project gene-expression-big-data \
    --bucket gene_datasets \
    --region us-central1 \
    --single-node \
    --master-machine-type e2-highmem-4 \
    --master-boot-disk-size 50 \
    --optional-components=JUPYTER \
    --image-version 2.1-debian11 \
    --initialization-actions \
    gs://dataproc-initialization-actions/python/pip-install.sh \
    --metadata 'PIP_PACKAGES=google-cloud-storage scikit-learn pandas plotly pydeseq2 openpyxl numpy==1.26.4 fsspec gprofiler-official pyarrow' \
    --enable-component-gateway \
    --max-idle 7200s --max-age 18000s

#for kai on windows
gcloud dataproc clusters create mycluster-2 --project gene-expression-big-data --bucket gene_datasets --region us-central1 --single-node --master-machine-type e2-highmem-4 --master-boot-disk-size 50 --optional-components=JUPYTER --image-version 2.1-debian11 --initialization-actions gs://dataproc-initialization-actions/python/pip-install.sh --metadata "PIP_PACKAGES=google-cloud-storage scikit-learn pandas" --enable-component-gateway --max-idle 7200s --max-age 18000s


# -------------------------------------------------------------
# 2. OPEN SSH TUNNEL (run in a LOCAL terminal, keep it open)
#    Forwards cluster JupyterLab (port 8123) → localhost:8890
# -------------------------------------------------------------
gcloud compute ssh mycluster-m \
    --project=gene-expression-big-data \
    --zone=us-central1-a \
    -- -L 127.0.0.1:8890:localhost:8123 -N

## For Kai on Windows
gcloud compute ssh mycluster-2-m --project=gene-expression-big-data --zone=us-central1-a -- -L 127.0.0.1:8890:localhost:8123 -N


# -------------------------------------------------------------
# 3. VERIFY TUNNEL (run in a second LOCAL terminal)
#    Should return HTML — if "connection refused", tunnel isn't up
# -------------------------------------------------------------
curl http://127.0.0.1:8890/gateway/default/jupyter

# for kai
http://localhost:8890/gateway/default/jupyter/lab


# -------------------------------------------------------------
# 4. CONNECT VS CODE
#    - Open your .ipynb in VS Code
#    - Click kernel selector (top right) → Select Another Kernel
#      → Existing Jupyter Server
#    - Enter: http://127.0.0.1:8890/gateway/default/jupyter
#    - Leave token/password blank
#    - Select the PySpark kernel
# -------------------------------------------------------------


# -------------------------------------------------------------
# 5. VERIFY SPARK IS CONNECTED (run in a notebook cell)
# -------------------------------------------------------------
# spark
# Expected output:
#   SparkSession - hive
#   Master: yarn
#   AppName: PySparkShell


# -------------------------------------------------------------
# 6. LOAD DATA FROM GCS (example)
# -------------------------------------------------------------
# from pyspark.sql import SparkSession
# spark = SparkSession.builder.appName("GTEx Preview").getOrCreate()
# df = spark.read.parquet("gs://gene_datasets/GTEx_tissue_expression.parquet")
# df.select(df.columns[:5]).limit(5).toPandas()


# -------------------------------------------------------------
# NOTES
# - Cluster auto-deletes after 2hrs idle or 5hrs total
# - Always save your .ipynb locally before the cluster expires
# - Keep the tunnel terminal (step 2) open while working
# - To back up notebook to GCS:
#     gsutil cp my_notebook.ipynb gs://gene_datasets/
# -------------------------------------------------------------


