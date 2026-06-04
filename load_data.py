# Import packages
import pandas as pd
from google.cloud import bigquery

PROJECT_ID = 'gene-expression-big-data'
client = bigquery.Client(project=PROJECT_ID)

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


#helper for pipeline_pull_opentargets
def unwrap_ids(cell):
    '''handles checking whta format the data comes in as from Big Query
    three cases are BQ nested dict, already in a list or tuple so convert 
    directly no unwrapping needed, and none will return empty list 
    so nothing breaks later'''
    if cell is None:
        return []
    if isinstance(cell, dict) and 'list' in cell:
        return [item['element'] for item in cell['list'] if 'element' in item]
    if isinstance(cell, (list, tuple)):
        return list(cell)
    return []

def pipeline_pull_opentargets():
    # pull it in and transpose since the original format was genes as cols
    gene_expression_matrix = pd.read_parquet('data/gene_expression_matrix.parquet')

    cpm = (gene_expression_matrix
    .drop(columns=['gene_name', 'gene_biotype'])
    .set_index('gene_id')
    .T)
    
    cpm.index.name = 'sample'

    # Add categorical sample_type column
    cpm.insert(0, 'sample_type', pd.Categorical(cpm.index.str.split('_').str[0]))

    # save cpm for analysis notebook
    cpm.to_parquet('data/sample_by_gene.parquet')

    gene_list = [c for c in cpm.columns if c.startswith('ENSG')]


    # SQL for bigquery
    ot_query = """
    SELECT
    a.targetId AS gene_id,
    dis.name AS disease_name,
    dis.id AS disease_id,
    dis.therapeuticAreas AS therapeutic_area_ids,
    a.associationScore AS association_score
    FROM `open-targets-prod.platform.association_overall_direct` AS a
    JOIN `open-targets-prod.platform.disease` AS dis
    ON a.diseaseId = dis.id
    WHERE a.targetId IN UNNEST(@gene_list)
    AND a.associationScore > 0.1
    """

    job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ArrayQueryParameter("gene_list", "STRING", gene_list)]
    )
    disease_associations = client.query(ot_query, job_config=job_config).to_dataframe()

    # therapeutic area IDs are nested lists, so we need to unwrap them and get unique IDs
    # Collect unique therapeutic area IDs
    all_ta_ids = set()
    for cell in disease_associations['therapeutic_area_ids'].dropna():
        for ta in unwrap_ids(cell):
            all_ta_ids.add(ta)
    all_ta_ids = list(all_ta_ids)

    # Look up readable names as therapeutic area query
    ta_query = """
    SELECT id, name
    FROM `open-targets-prod.platform.disease`
    WHERE id IN UNNEST(@ta_ids)
    """
    ta_job_config = bigquery.QueryJobConfig(
        query_parameters=[bigquery.ArrayQueryParameter("ta_ids", "STRING", all_ta_ids)]
    )
    ta_lookup = client.query(ta_query, ta_job_config).to_dataframe()
    ta_dict = dict(zip(ta_lookup['id'], ta_lookup['name']))

    # Map IDs to names
    disease_associations['therapeutic_areas'] = disease_associations['therapeutic_area_ids'].apply(
        lambda cell: [ta_dict.get(i, i) for i in unwrap_ids(cell)]
    )

    # Save as parquet 
    disease_associations.to_parquet('data/disease_associations.parquet')



if __name__ == "__main__":
    #raw_data() # only run once to create original parquet file
    pipeline_pull_opentargets() # run for grabbing the opentargets data

