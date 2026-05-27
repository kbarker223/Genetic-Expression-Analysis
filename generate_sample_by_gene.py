# Import packages
import pandas as pd

# Read data and display head
raw_data = pd.read_excel("CPM-matrix_Tony.xlsx")

# Transpose: genes become columns, samples become rows
clean_data = (
    raw_data
    .drop(columns=[*raw_data.filter(like='count').columns, 'gene_biotype', 'gene_name'])
    .set_index('gene_id')
    .rename(index={'SNHG14': 'Target(SNHG14)'}, columns={'hSR6_1_cpm': 'Control_1(HSR6)', 'hSR6_2_cpm': 'Control_2(HSR6)'})
    .T
)
clean_data.index.name = 'sample'

# Add categorical sample_type column
clean_data.insert(0, 'sample_type', pd.Categorical(clean_data.index.str.split('_').str[0]))

# Save new data as parquet
clean_data.to_parquet("sample_by_gene.parquet")