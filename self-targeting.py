#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read CRISPR Self-Targeting database
df = pd.read_csv('self-target.csv')

# Get feature table of protospacer region
for i, start, end in zip(df['Refseq ID'], df['Proto-spacer Start'], df['Proto-spacer End']):
    result = Entrez.efetch(db='nuccore', id=i, seq_start=start, seq_stop=end, rettype='ft')
