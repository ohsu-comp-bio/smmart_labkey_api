import argparse
import warnings
import pandas as pd
import numpy as np
import labkey
from labkey_client import SMMARTLabkey
import json

# Command line arguments:
#     --id  Create a summary table for a specific patient in addition to the summary of all the patients. Default is 101.
#     --out_dir  Where to save the output file. Default is Desktop.
# Additional markers and genes can be included by adding them to the appropriate json config dictionary. Any value that exists in the LabKey tables should work.
# Two important options are 'cnv_filter' and 'mutation_filter'. Both are set to 'true' by default, so the LabKey query will only select records when the 'Reported' has a true value.
# It is assumed the combination of 'Participant ID' and 'Date/Collection Date' can be used as a key when merging/joining records from different tables. 
# 'Bems IDs' would be useful for this, but unfortunately the 'Bems ID' is often missing for a record. 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', default = 0, type = int,  help = 'Create a summary table for a specific patient in addition to the summary of all the patients.')
    parser.add_argument('--study', default = 'amtec', type = str,  help = 'Create a summary table for patients in a spcefic study.')
    parser.add_argument('--out_dir', default = '~/Desktop', type = str, help = 'Where to save the output file.')
    args = parser.parse_args()

    # Use to join tables
    INDEX = ['Participant ID', 'Date']

    with open ('patient_summaries_config.json') as json_config:
        config = json.load(json_config)
        # Query filters.
        id_range = config['filters']['id_range']
        # TODO (hocumj 8/14/2019): Handle more studies than AMTEC.
        study_ids = config['studies']['amtec']
        cnv_reported_filter = config['filters']['cnv_reported_filter']
        mutation_reported_filter = config['filters']['mutation_reported_filter']
        # Columns of interest for summary output.
        # The first element in these lists should be the column to be summarized. e.g., 'Gene' or 'Marker'.
        summary_ihc_columns = config['summary_columns']['ihc_columns']
        summary_cnv_columns = config['summary_columns']['cnv_columns']
        summary_mutation_columns = config['summary_columns']['mutation_columns']
        # Genes/proteins of interest. 
        summary_ihc_markers = config['summary_columns']['ihc_markers']
        summary_cnv_genes = config['summary_columns']['cnv_genes']
        summary_mutation_genes = config['summary_columns']['mutation_genes']
        # Mapping values.
        positive_marker_values = config['marker_values']['positive']
        negative_marker_values = config['marker_values']['negative']
        # nan_marker_values = config['marker_values']['nan']

    # SMMARTLabkey object for queries.
    labkey_client = SMMARTLabkey(config['labkey']['server'], config['labkey']['project'], config['labkey']['context'])

    # Convert genes and markers of interest into a string suitable for the labkey queries.
    query_markers = ';'.join(summary_ihc_markers)
    query_cnvs = ';'.join(summary_cnv_genes)
    query_mutations = ';'.join(summary_mutation_genes)

    # Get tables.
    patient_markers = labkey_client.query_study(query = 'clinical_receptor_status', filter_array = [labkey.query.QueryFilter('ParticipantID', id_range[0], 'gte'), labkey.query.QueryFilter('ParticipantID', id_range[1], 'lt'), labkey.query.QueryFilter('marker_id/display_name', query_markers, 'in')])
    if cnv_reported_filter == 'true':
        patient_cnv = labkey_client.query_study(query = 'sample_genetrails_copy_number_variant', filter_array = [labkey.query.QueryFilter('ParticipantID', id_range[0], 'gte'), labkey.query.QueryFilter('ParticipantID', id_range[1], 'lt'), labkey.query.QueryFilter('reported', 'true', 'eq'), labkey.query.QueryFilter('gene_id/display_name', query_cnvs, 'in')])
    else:
        patient_cnv = labkey_client.query_study(query = 'sample_genetrails_copy_number_variant', filter_array = [labkey.query.QueryFilter('ParticipantID', id_range[0], 'gte'), labkey.query.QueryFilter('ParticipantID', id_range[1], 'lt'), labkey.query.QueryFilter('gene_id/display_name', query_cnvs, 'in')])   
    if mutation_reported_filter == 'true':
        patient_seq_mutations = labkey_client.query_study(query = 'sample_genetrails_sequence_variant', filter_array = [labkey.query.QueryFilter('ParticipantID', id_range[0], 'gte'), labkey.query.QueryFilter('ParticipantID', id_range[1], 'lt'), labkey.query.QueryFilter('reported', 'true', 'eq'), labkey.query.QueryFilter('gene_id/display_name', query_mutations, 'in')])
    else:
        patient_seq_mutations = labkey_client.query_study(query = 'sample_genetrails_sequence_variant', filter_array = [labkey.query.QueryFilter('ParticipantID', id_range[0], 'gte'), labkey.query.QueryFilter('ParticipantID', id_range[1], 'lt'), labkey.query.QueryFilter('gene_id/display_name', query_mutations, 'in')])

    # Make date related column names consistent and make sure there are no ambigious duplicate rows.
    patient_markers['Participant ID'] = patient_markers['Participant ID'].astype('int64')
    patient_markers['Date'] = pd.to_datetime(patient_markers['Date'])
    patient_markers['Status'].replace(positive_marker_values, 'Positive', inplace = True)
    patient_markers['Status'].replace(negative_marker_values, 'Negative', inplace = True)
    patient_markers.drop_duplicates(['Participant ID', 'Date', 'Marker'], inplace = True, keep = 'last')

    patient_cnv['Participant ID'] = patient_cnv['Participant ID'].astype('int64')
    patient_cnv.rename(columns = {'Collection Date': 'Date'}, inplace = True)
    patient_cnv['Date'] = pd.to_datetime(patient_cnv['Date'])
    patient_cnv.drop_duplicates(['Participant ID', 'Date', 'Gene'], inplace = True, keep = 'last')

    patient_seq_mutations['Participant ID'] = patient_seq_mutations['Participant ID'].astype('int64')
    patient_seq_mutations.rename(columns = {'Collection Date': 'Date'}, inplace = True)
    patient_seq_mutations['Date'] = pd.to_datetime(patient_seq_mutations['Date'])
    patient_seq_mutations.drop_duplicates(['Participant ID', 'Date', 'Gene', 'Position Start', 'Position End', 'Reference Base', 'Variant Base'], inplace = True, keep = 'last')

    # Pivot so applying logic across patient biopsy is more straightforward. 
    patient_markers = generate_pivot_table(patient_markers, INDEX, summary_ihc_markers, summary_ihc_columns, sort_asc = True)
    patient_markers['Subtype'] = patient_markers.apply(set_subtype, axis = 1)

    patient_cnv = generate_pivot_table(patient_cnv, INDEX, summary_cnv_genes, summary_cnv_columns, sort_asc = True)

    patient_seq_mutations = generate_pivot_table(patient_seq_mutations, INDEX, summary_mutation_genes, summary_mutation_columns, sort_asc = True)
    
    # Using outer joins so if the dates do match between tables that information can be represented in 1 row while not losing rows from the original tables. 
    patient_summaries = patient_markers.merge(patient_cnv, left_index = True, right_index = True, how = 'outer')
    patient_summaries = patient_summaries.merge(patient_seq_mutations, left_index = True, right_index = True, how = 'outer')
    patient_summaries.sort_index(ascending = True, inplace = True)

    # TODO (hocumj 8/14/2019): Handle more studies than AMTEC. Grab list of ids by matching args.study string to dictionary key.
    if args.study == 'amtec':
        study_patient_summaries = patient_summaries.loc[study_ids, :]
        study_patient_summaries.to_csv(args.out_dir + '/' + str(args.study) + '_summaries.tsv', sep = '\t')
    else:
        patient_summaries.to_csv(args.out_dir + '/patient_summaries.tsv', sep = '\t')
    # Return a table with only the patient of interest.
    if args.id > 0:
        patient_summary = patient_summaries.loc[args.id, :]
        patient_summary.to_csv(args.out_dir + '/' + str(args.id) + '_summary.tsv', sep = '\t')

def set_subtype(row):
    subtype = 'Undetermined'
    if (row['ER'] == 'Positive' or row['PR'] == 'Positive'):
        if (row['HER2'] == 'Negative'):
            subtype = 'Luminal A'
        elif (row['HER2'] == 'Positive'):
            subtype = 'Luminal B'
    elif (row['ER'] == 'Negative' and row['PR'] == 'Negative' and row['HER2'] == 'Negative'):
        subtype = 'Triple Negative'
    # TODO (hocumj 8/14/2019) Provide a warning stating the which value is missing in the config file. 
    else:
        warnings.warn('Not enough information provided to determine subtype. Consider adding more "marker values" in the json config file.')
    return subtype

def generate_pivot_table(df, index, genes_of_interest, columns_of_interest, sort_asc):
    # Grab the columns to use for the index and the columns to summarize the table.
    pt = df.loc[: , index + columns_of_interest]
    # Slice so only the genes of interest are included (should already be filtered by the labkey query). 
    pt = pt.loc[pt[columns_of_interest[0]].isin(genes_of_interest)]
    pt = pt.pivot_table(index = index, columns = columns_of_interest[0], aggfunc = 'last')
    # Collapse the multilevel column index to 1 level.
    # col[::-1] reverses the tuple so the gene/protein name is the beginning of the column name. 
    if len(columns_of_interest) > 2:
        pt.columns = [' '.join(col[::-1]).strip() for col in pt.columns.values]
    # Keep the column name simpler in this case. 
    else:
        pt.columns = [col[1] for col in pt.columns.values]
    pt.sort_index(axis = 1, ascending = sort_asc, inplace = True)
    return pt

if __name__ == '__main__':
    main()