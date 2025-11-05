from utils.rdkit_utils import convert_to_mol_from_smiles
from utils.chemical_database import csv_to_df, standardized_df, process_smiles_df

# Handling user input (direct input of smiles, csv file)
def input_smiles(smiles):
    mol = convert_to_mol_from_smiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES input: {smiles}')
    return mol

def input_csv(file_path):
    df = csv_to_df(file_path)
    df = standardized_df(df)
    df = process_smiles_df(df)
    return df