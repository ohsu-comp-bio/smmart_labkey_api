#!/usr/bin/env python

import argparse
import json
import pandas as pd

from labkey_client import SMMARTLabkey

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("table_name", help="Table name")
    parser.add_argument("output", help="Output file name")
    parser.add_argument("--sample", help="Sample name to filter on")
    parser.add_argument("--format", default="tsv", help="Output file format: TSV or CSV")
    args = parser.parse_args()

    # SMMART Labkey server info
    server = "smmart-research.ohsu.edu"
    project = "SMMART Research"
    contextPath = "SMMARTResearch"

    # Query metadata
    labkey = SMMARTLabkey(server, project, contextPath)
    df = labkey.query_study(query=args.table_name)

    # Filter by sample name
    if args.sample:
        df = df[df['Sample Code'] == args.sample]

    # Write metadata table
    if args.format == "csv":
        df.to_csv(args.output, sep=",", index=False) 
    else:
        df.to_csv(args.output, sep="\t", index=False)
