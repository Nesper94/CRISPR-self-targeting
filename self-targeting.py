#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read CRISPR Self-Targeting database
df = pd.read_csv('self-target.csv')

# Create auxiliary lists
products = []
protein_ids = []

# Get feature table of protospacer region
for i, start, end in zip(df['Refseq ID'], df['Proto-spacer Start'], df['Proto-spacer End']):
    try:
        result = Entrez.efetch(db='nuccore', id=i, seq_start=start, seq_stop=end, rettype='ft')
        feature_table = result.read().split('\t')
    except:
        products.append('No info')
        protein_ids.append('No info')
        continue

    # Get product name
    try:
        product_idx = feature_table.index('product')
        products.append(feature_table[product_idx + 1].strip())
    except:
        products.append('No protein')

    # Get protein ID
    try:
        protein_id_idx = feature_table.index('protein_id')
        protein_ids.append(feature_table[protein_id_idx + 1].replace('ref|','').strip('\n|'))
    except(ValueError):
        protein_ids.append('No protein')

# Modify original dataframe
df['Product'] = products
df['Protein id'] = protein_ids

# Lists for new columns
GO_function = []
GO_process = []
genes = []
EC_numbers = []

# Download protein info
for i in df['Protein id']:
    if i == 'No info':
        GO_function.append('No NCBI response')
        GO_process.append('No NCBI response')
        genes.append('No NCBI response')
        EC_numbers.append('No NCBI response')
    elif i == 'No protein':
        GO_function.append('No protein')
        GO_process.append('No protein')
        genes.append('No protein')
        EC_numbers.append('No protein')
    else:
        try:
            result = Entrez.efetch(db='protein', id=i, rettype='gp')
        except:
            GO_function.append('No NCBI response')
            GO_process.append('No NCBI response')
            genes.append('No NCBI response')
            EC_numbers.append('No NCBI response')
            continue

        result = result.read()
        res_list = result.split('\n')
        func_id = 'No function id'
        process_id = 'No process id'
        gene = 'No gene name'
        ec_number = 'No EC number'

        for i in res_list:
            if 'GO_function' in i:
                func_id = i.split()[0].replace('/GO_function="','')
            if 'GO_process' in i:
                process_id = i.split()[0].replace('/GO_process="','')
            if 'gene=' in i:
                gene = i.split()[0].split('"')[1]
            if 'EC_number' in i:
                ec_number = i.split()[0].split('"')[1]

        GO_function.append(func_id)
        GO_process.append(process_id)
        genes.append(gene)
        EC_numbers.append(ec_number)

df['GO Function'] = GO_function
df['GO process'] = GO_process
df['Gene'] = genes
df['EC number'] = EC_numbers

# Save table to file
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
