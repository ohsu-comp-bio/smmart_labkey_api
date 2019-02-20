#!/usr/bin/env python

# This script targets the client api version 0.4.0 and later

import labkey
import json
import pandas as pd

class SMMARTLabkey():
    def __init__(self, server, project, contextPath):
        self.server_context = labkey.utils.create_server_context(server, project, contextPath, use_ssl=True)

    def query_labkey(self, schema, query, filter_array):
        q = labkey.query.select_rows(
            server_context=self.server_context,
            schema_name=schema,
            query_name=query,
            filter_array=filter_array
        )

        return q

    def query_doc_store(self, filter_array=None):
        schema = "lists"
        query = "documentstore"
        doc = self.query_labkey(schema, query, filter_array)
        return doc

    def query_sequencing_meta(self, filter_array=None):
        schema = "study"
        query = "sequencing_metadata"
        meta = self.query_labkey(schema, query, filter_array)
        return meta

    def query_genetrails_variants(self, filter_array=None):
        schema = "study"
        query = "sample_genetrails_sequence_variant"
        variants = self.query_labkey(schema, query, filter_array)
        return variants

    def query_genetrails_copy_number_variant(self, filter_array=None):
        schema = "study"
        query = "sample_genetrails_copy_number_variant"
         variants = self.query_labkey(schema, query, filter_array)
         return variants

    def _create_dataframe(self, q):

        # Get visible columns and header
        columns = q['columnModel']
        header = {}
        for column in columns:
            if not column['hidden']:
                header[column['dataIndex']] = column['header']

        # Select rows
        rows = q['rows']
        
        # Load rows into pandas dataframe
        df = pd.DataFrame(rows)
        df = df[header.keys()]
        df.rename(header, axis='columns', inplace=True)

        return df
