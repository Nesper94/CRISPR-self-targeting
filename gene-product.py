#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read CRISPRminer Self-Targeting database
df = pd.read_csv('self-target.csv')

# Create auxiliary lists
products = []
protein_ids = []

# Get feature table of protospacer region
for i, start, end in zip(df['Refseq ID'], df['Proto-spacer Start'], df['Proto-spacer End']):
    try:
        handle = Entrez.efetch(db='nuccore', id=i, seq_start=start, seq_stop=end, rettype='ft')
        feature_table = handle.read().split('\t')
        handle.close()
    except Exception as e:
        print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end)
        products.append('Error: '+str(e))
        protein_ids.append('Error: '+str(e))
        continue

    # Get product name
    try:
        product_idx = feature_table.index('product')
        products.append(feature_table[product_idx + 1].strip())
    except Exception as e:
        print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end)
        products.append('Error: '+str(e))

    # Get protein ID
    try:
        protein_id_idx = feature_table.index('protein_id')
        protein_ids.append(feature_table[protein_id_idx + 1].replace('ref|','').strip('\n|'))
    except Exception as e:
        print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end)
        protein_ids.append('Error: '+str(e))

# Modify original dataframe
df['Product'] = products
df['Protein id'] = protein_ids

# Save table to file
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
