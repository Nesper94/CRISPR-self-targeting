#!/usr/bin/env python3

# Download sequences of hypothetical proteins and write them in a FASTA file

import pandas as pd
from Bio import Entrez

Entrez.email = 'juan.arboleda2@udea.edu.co'

cols = ['Product', 'Protein id']
df = pd.read_csv('self-target-proteins.tsv', sep='\t', usecols=cols)

with open('hypothetical.fasta', 'a') as file:

    for i in df[df['Product'].str.contains('hypothetical')]['Protein id']:
        try:
            handle = Entrez.efetch(db='protein', id=i, rettype='fasta')
            seq = handle.read()
            handle.close()
        except:
            continue

        file.write(seq)
