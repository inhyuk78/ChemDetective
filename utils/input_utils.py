from utils.rdkit_utils import convert_to_mol_from_smiles
from utils.chemical_database import csv_to_df, standardized_df, process_smiles_df

# Handling user input (direct input of smiles, csv file)
def input_smiles(smiles):
    mol = convert_to_mol_from_smiles(smiles)
    return mol

#def input_smiles(smiles):
    #while True:
        #try: 
            #mol = convert_to_mol_from_smiles(smiles)
            #if mol != None:
                #return mol
            #else:
                #print(f'{smiles} is not a valid SMILES. Please check and try again.')
        #except Exception as e:
            #print(f'An error occurred: {e}. Please try again.')


def input_csv(file_path):
    print('DEBUG - reading:', file_path)
    df = csv_to_df(file_path)
    df = standardized_df(df)
    df = process_smiles_df(df)
    return df

