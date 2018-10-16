#!/home/liam/anaconda3/bin/python3.6

import pubchempy as pcp
import rdkit.Chem.AllChem as rdkit
import os

# Import pandas
import pandas as pd

# Assign spreadsheet filename to `file`
file = '/home/liam/software/iupac_to_smiles_converter/DibromoCompoundsOverview.xlsx'

# Load spreadsheet
xl = pd.ExcelFile(file)

# Load a sheet into a DataFrame by name: df1
df1 = xl.parse('Results')

for i in range(670):

    iupac = df1.iloc[i][0]
    try:
        # Generate the smiles code from the PubChem repo
        c = pcp.get_compounds(iupac, 'name')
        a = c[0]
        smiles = a.canonical_smiles
        print('%04d  %s' % (i+2, smiles))

    except IndexError as e:
        print('%04d  ERROR : %s' % (i+2, str(e)))

