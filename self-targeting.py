#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read CRISPR Self-Targeting database
df = pd.read_csv('self-target.csv')

# Get feature table of protospacer region
for i in df['Refseq ID']:
    result = Entrez.efetch(db='nuccore', id=i, seq_start=1146152, seq_stop=1146181, rettype='ft')
