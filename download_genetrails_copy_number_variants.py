#!/usr/bin/env python

import argparse
import json
import pandas as pd

from labkey_client import SMMARTLabkey

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("chromosomes", help="Labkey Chromosome lookup list")
    parser.add_argument("--sample", help="Sample name to filter on")
    args = parser.parse_args()

    # Load chromosome lookup list
    with open(args.chromosomes, "r") as c:
        chrs = json.load(c)

    # SMMART Labkey server info
    server = "bcclabkey.ohsu.edu"
    project = "SMMART Research"
    contextPath = "SMMARTResearch"

    # Query metadata
    labkey = SMMARTLabkey(server, project, contextPath)
    variants = labkey.query_genetrails_copy_number_variant()
    df = labkey._create_dataframe(variants)

    # Map integer codes to lookup list
    df['Chromosome'] = df['Chromosome'].astype(str).replace(chrs)

    # Filter by sample name
    if args.sample:
        df = df[df['Sample Code'] == args.sample]

    # Write metadata table
    df.to_csv("genetrails-variants.tsv", sep="\t", index=False) 
