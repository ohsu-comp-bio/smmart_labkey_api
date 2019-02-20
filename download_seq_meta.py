#!/usr/bin/env python

import argparse
import json
import pandas as pd

from labkey_client import SMMARTLabkey

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("lookup_list", help="Labkey key for mapping integer values")
    args = parser.parse_args()

    # Load lookup list
    with open(args.lookup_list, "r") as l:
        lookup = json.load(l)

    # SMMART Labkey server info
    server = "bcclabkey.ohsu.edu"
    project = "SMMART Research"
    contextPath = "SMMARTResearch"

    # Query metadata
    labkey = SMMARTLabkey(server, project, contextPath)
    meta = labkey.query_sequencing_meta()
    df = labkey._create_dataframe(meta)

    # Map integer codes to lookup list
    # Better way to do this other than making copies?
    for k,v in lookup.items():
        df[k] = df[k].map('{:,.0f}'.format).astype(str).replace(v)

    # Write metadata table
    df.to_csv("df-meta-query.csv", sep=",", index=False) 
