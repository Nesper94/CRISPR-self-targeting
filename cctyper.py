#!/usr/bin/env python3

from Bio import Entrez
import pandas as pd
import subprocess

Entrez.email = 'juan.arboleda2@udea.edu.co'

# Read database
df = pd.read_csv('self-target-proteins.tsv', sep='\t')

# Number of bases adjacent to spacers
upstream = 9000
downstream = 9000

# List of CRISPR Types
crisprt = []

# Get the sequence from NCBI
for id_, start, end in zip(df['Refseq ID'], df['Spacer Start'], df['Spacer End']):
    try:
        handle = Entrez.efetch(db='nuccore', id=id_, seq_start=start-upstream, seq_stop=end+downstream, rettype='fasta')
        result = handle.read()
        handle.close()

    except Exception as e:
        print(e, 'Refseq ID:', id_, 'Start:', start, 'End:', end)
        crisprt.append('Error: '+str(e))
        continue

    # Export sequence to file
    with open('seq.fasta', 'w') as fasta:
        fasta.write(result)

    # Run CCTyper
    process = subprocess.Popen('conda run -n cctyper cctyper --no_plot seq.fasta result'.split())
    process.wait()

    # Import results
    try:
        res = pd.read_csv('result/cas_operons.tab', sep='\t')
    except FileNotFoundError:
        crisprt.append('No CRISPR type prediction')
        continue
    finally:
        # Delete result folder
        process = subprocess.Popen('rm -rf result/'.split())
        process.wait()
    
    # Save to list
    crisprt.extend(res['Prediction'].to_list())

# Create new column from list
df['CRISPR Type'] = crisprt
# Clean result
df['CRISPR Type'] = df['CRISPR Type'].str.replace('\nName: Prediction, dtype: object', '').str.replace('0    ','')

# Save dataframe
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
