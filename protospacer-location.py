#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez
import parse
from sys import stderr

Entrez.email = 'juan.arboleda2@udea.edu.co'

df = pd.read_csv('self-target-proteins.tsv', sep='\t')

def get_cds_info(record):
    '''Returns CDS chain (same or complement), start and end as a tuple.'''

    substring = parse.search('[location={}({}..{})]', record)
    if substring == None:
        substring = parse.search('[location={}..{}]', record)
        substring.fixed = ('same',) + substring.fixed

    return substring.fixed

def ps_loc(protospacer_start, protospacer_end, cds_info):
    '''Returns as percentage the relative position of protospacer start with
    respect to its coding sequence.'''

    cds_loc = cds_info[0]
    cds_start = int(cds_info[1])
    cds_end = int(cds_info[2])
    cds_len = cds_end - cds_start

    if cds_loc == 'complement':
        return (cds_end - protospacer_end)/cds_len*100
    else:
        return (protospacer_start - cds_start)/cds_len*100

# Create list to store locations
protospacer_loc = []

with open('cds.fasta', 'a') as file:

    # Get CDS info
    for i, start, end in zip(df['Refseq ID'], df['Proto-spacer Start'], df['Proto-spacer End']):
        try:
            handle = Entrez.efetch(db='nuccore', id=i, seq_start=start,
                                   seq_stop=end, rettype='fasta_cds_na')
            cds = handle.read()
            handle.close()
        except Exception as e:
            print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end, file=stderr)
            protospacer_loc.append('Error: '+str(e))
            continue

        # Write CDS to fasta
        file.write(cds)

        # Get protospacer relative location
        try:
            cds_info = get_cds_info(cds)
            protospacer_loc.append(ps_loc(int(start), int(end), cds_info))
        except Exception as e:
            print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end, file=stderr)
            protospacer_loc.append('Error: '+str(e))

# Create new column from list
df['Proto-spacer loc'] = protospacer_loc

# Save dataframe
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
