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

    def query_study(self, query, filter_array=None):
        schema = "study"
        study_table = self.query_labkey(schema, query, filter_array)
        study_table = self._create_dataframe(study_table)
        return study_table

    def query_lookup_lists(self,lookup):
        lists = self.query_labkey(schema=lookup['schema'], query=lookup['queryName'], filter_array=None)
        return lists['rows']

    def replace_lookup_lists(self, fields, df):
        lookups = {}
        key_cols = {}
        for field in fields:
            if 'lookup' in field:
                lists = self.query_lookup_lists(field['lookup'])
                lookups[field['lookup']['keyColumn']] = lists
                key_cols[field['lookup']['keyColumn']] = field['name']

        lookups.pop('ParticipantID')

        for lookup_key,lookup_val in lookups.items():
            if '/' in lookup_key:
                lookup_id = lookup_key.split("/")[1]
            else:
                lookup_id = lookup_key

            ids = {}
            for val in lookup_val:
                if 'display_name' in val.keys():
                    ids[val[lookup_id]] = val['display_name']

            if key_cols[lookup_key] in df.columns:
                df[key_cols[lookup_key]].replace(ids, inplace=True)

        return df
        
    def _create_dataframe(self, q):

        # Get visible columns and header
        columns = q['columnModel']
        fields = q['metaData']['fields']
        rows = q['rows']

        # 
        header = {}
        for column in columns:
            if not column['hidden']:
                header[column['dataIndex']] = column['header']

        # Load rows into pandas dataframe
        df = pd.DataFrame(rows)
        df = df[list(header.keys())]
        df = self.replace_lookup_lists(fields, df)
        df.rename(header, axis='columns', inplace=True)
        
        return df
