#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read Self-Targeting database
df = pd.read_csv('self-target-proteins.tsv', sep='\t')

# Lists for new columns
GO_function = []
GO_process = []
genes = []
EC_numbers = []

# Download protein info
for i in df['Protein id']:
    if 'Error:' in i:
        GO_function.append(i)
        GO_process.append(i)
        genes.append(i)
        EC_numbers.append(i)
    else:
        try:
            handle = Entrez.efetch(db='protein', id=i, rettype='gp')
            result = handle.read()
            handle.close()
        except Exception as e:
            err = 'Error: '+str(e)
            GO_function.append(err)
            GO_process.append(err)
            genes.append(err)
            EC_numbers.append(err)
            continue

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

# Create new columns in table
df['GO Function'] = GO_function
df['GO process'] = GO_process
df['Gene'] = genes
df['EC number'] = EC_numbers

# Save table to file
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
