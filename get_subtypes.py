#!/usr/bin/env python

import argparse
import json
import pandas as pd

from labkey_client import SMMARTLabkey

class ReceptorInfoError(Exception):
    pass

def subtypes_groups (receptors):
    if (receptors['HER2'] == "Positive"):
        group = "HER2"
    elif (receptors['ER'] == "Positive" or receptors['PR'] == "Positive"):
        group = "HR Positive"
    elif (receptors['ER'] == "Negative" & receptors['PR'] == "Negative" & receptors['HER2'] == "Negative"):
        group = "Triple Negative"
    else:
        raise ReceptorInfoError("Not enough information provided to determine subtype.")

    return group

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("data_type", help="Data type to filter samples: Exome or Transcriptome")
    parser.add_argument("--sample", help="Sample name to filter on")
    args = parser.parse_args()

    # SMMART Labkey server info
    server = "smmart-research.ohsu.edu"
    project = "SMMART Research"
    contextPath = None
    labkey = SMMARTLabkey(server, project, contextPath)

    # Get Sequencing Metadata
    meta = labkey.query_study(query='sequencing_metadata')

    # Filter by sample name
    if args.sample:
        meta = meta[meta['Sample Code'] == args.sample]

    # Filter by data type
    if args.data_type:
        meta = meta[meta['Protocol[Library_Type]'] == args.data_type]

    # Get Collection Dates from Biolibrary specimen table
    specimen = labkey.query_compbiostudy(query='vbiolibraryspecimens')
    
    # Get Receptor Status
    receptor = labkey.query_study(query='clinical_receptor_status')

    # Merge tables
    df = pd.merge(meta,specimen, how='inner', left_on="Sample[BEMS_ID]", right_on="BEMS ID")
    df = pd.merge(df,receptor, how='inner', left_on=["Collection Date","Participant ID_x"], right_on=["Date","Participant ID"])

    # Write Table
    cols = ['Participant ID', 'Collection Date', 'Library ID', 'BEMS ID','Marker','Status']
    df = df.loc[:,cols]
    df.drop_duplicates(inplace=True)
    df.to_csv("sample_subtypes.tsv", sep="\t", index=False)

    # Determine Subtype
    ## To Do
