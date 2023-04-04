#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez

# Read Self-Targeting database
df = pd.read_csv('self-target-proteins.tsv', sep='\t')

# Get Taxonomy ID
# Create auxiliary TaxIds list
TaxIDs = []

for i in df['Refseq ID']:
    try:
        handle = Entrez.esummary(db='nuccore', id=i)
        result = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(e, 'Refseq ID:', i)
        TaxIDs.append('No NCBI response')
        continue

    try:
        TaxIDs.append(result[0]['TaxId'].real)
    except Exception as e:
        print(e, 'Refseq ID:', i)
        TaxIDs.append('No TaxID')

df['TaxID'] = TaxIDs

# Get species and genus
species_l = []
genus_l = []

for i in df['Organism']:
    split = i.split()
    genus = split[0]
    species = ' '.join(split[0:2])
    species_l.append(species)
    genus_l.append(genus)

df['Species'] = species_l
df['Genus'] = genus_l

# Get additional taxonomic information
superkingdom_l = []
phylum_l = []
class_l = []
order_l = []
family_l = []

for tax_id in df['TaxID']:
    try:
        handle = Entrez.efetch(db='taxonomy', id=tax_id)
        result = Entrez.read(handle)
        handle.close()
    except Exception as e:
        err = 'Error: ' + str(e)
        print(err)
        superkingdom_l.append(err)
        phylum_l.append(err)
        class_l.append(err)
        order_l.append(err)
        family_l.append(err)
        continue

    superkingdom_ = 'Error: No info'
    phylum_ = 'Error: No info'
    class_ = 'Error: No info'
    order_ = 'Error: No info'
    family_ = 'Error: No info'

    for i in result[0]['LineageEx']:
        if i['Rank'] == 'superkingdom':
            superkingdom_ = i['ScientificName']

        elif i['Rank'] == 'phylum':
            phylum_ = i['ScientificName']

        elif i['Rank'] == 'class':
            class_ = i['ScientificName']

        elif i['Rank'] == 'order':
            order_ = i['ScientificName']

        elif i['Rank'] == 'family':
            family_ = i['ScientificName']

    superkingdom_l.append(superkingdom_)
    phylum_l.append(phylum_)
    class_l.append(class_)
    order_l.append(order_)
    family_l.append(family_)

# Add new columns to table
df['Superkingdom'] = superkingdom_l
df['Phylum'] = phylum_l
df['Class'] = class_l
df['Order'] = order_l
df['Family'] = family_l

# Save table to file
df.to_csv('self-target-proteins.tsv', index=False, sep='\t')
