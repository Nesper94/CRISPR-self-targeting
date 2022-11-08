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
    result = Entrez.efetch(db='nuccore', id=i, seq_start=start, seq_stop=end, rettype='ft')
    feature_table = result.read().split('\t')

    # Get product name
    product_idx = feature_table.index('product')
    products.append(feature_table[product_idx + 1].strip())

    # Get protein ID
    protein_id_idx = feature_table.index('protein_id')
    protein_ids.append(feature_table[protein_id_idx + 1].replace('ref|','').strip('\n|'))

# Modify original dataframe
df['Product'] = products
df['Protein id'] = protein_ids
