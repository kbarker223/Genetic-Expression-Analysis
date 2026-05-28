# Import packages
import pandas as pd

def raw_data():
    """
    Used only once to create parquet file with dropped count columns. 
    Also renamed target gene and control treatments for clarity. 
    Raw data deleted to conserve space.
    """
    # Read raw data
    print("Reading raw data...")
    raw_data = pd.read_excel("data/CPM-matrix_Tony.xlsx")

    print("Cleaning...")
    # Find gene_id for target gene SNHG14
    target_gene_id = raw_data[raw_data['gene_name'] == 'SNHG14']['gene_id'].iloc[0]

    # Transpose: genes become columns, samples become rows
    clean_data = (
        raw_data
        .drop(columns=[*raw_data.filter(like='count').columns])
        .set_index('gene_id')
        .rename(index={target_gene_id: 'Target(SNHG14)'}, columns={'hSR6_1_cpm': 'Control_1(HSR6)', 'hSR6_2_cpm': 'Control_2(HSR6)'})
        .reset_index()
    )

    print("Saving parquet...")
    clean_data.to_parquet("data/gene_expression_matrix.parquet")
    print("Successfully saved as data/gene_expression_matrix.parquet")


if __name__ == "__main__":
    raw_data() # only run once to create original parquet file


