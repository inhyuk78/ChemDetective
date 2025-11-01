import pandas as pd
import numpy as np
from rdkit.Chem import PandasTools, Draw
from rdkit_utils import convert_to_mol_from_smiles, check_fg_in_mol,find_mw_in_mol, visualize_mol


def csv_to_df(file_path):
    '''
    Converts CSV file into pandas dataframe
    Parameter:
        file_path : CSV = File (CSV) containing drug names, smiles
    Return:
        df : pd.DataFrame = DataFrame containing CSV file data if file found
        FileNotFoundError = '{file_path} is not recognized. Please try another file.' if file not found
    ''' 
    while True:
        try:
            df = pd.read_csv(file_path)
            return df
        except FileNotFoundError:
             print(f'{file_path} is not recognized. Please try another file.')
             file_path = input('Enter a valid CSV file: ')

def standardized_df(df, fillna_value=None): 
    '''
    Converts dataframe into standardized form for analysis
    Parameter:
        df : pd.DataFrame = Raw user input CSV file converted DataFrame
    Return:
        df : pd.DataFrame = Standardized DataFrame
    '''
    while True:
        if df.shape[1] != 2: 
            print('Incorrect file format. Make sure your file has 2 columns (column 1: Drug Names, column 2: SMILES)')
        else:
            # converts first row into column, removes white spaces & turns into lowercase column names
            df.columns = df.iloc[0]
            df.columns = [col.strip().lower() for col in df.columns] 

            # rename first & second column to 'Drug_Name', 'SMILES' respectively
            # set 'Drug_Name' as index column
            df.columns.values[0] = 'Drug_Name' 
            df.columns.values[1] = 'SMILES'
            df = df.set_index('Drug_Name')

            # fills empty values in df with 'fillna_value'
            if fillna_value is not None:
                df = df.fillna(fillna_values)

            return df

def process_smiles_df(df):
    '''
    Adds new columns (Molecule, functional groups, MW) containing data processed from smiles
    Parameter:
        df : pd.DataFrame = Standardized DataFrame
    Return:
        df : pd.DataFrame = DataFrame containing new columns
    '''
    
    # Add functional groups column
    df['Mol'] = df['SMILES'].apply(convert_to_mol_from_smiles)
    df['Functional Groups'] = df['Mol'].apply(check_fg_in_mol)
    df['Molecule Image'] = df['Mol'].apply(visualize_mol)

    # Add new column with MW
    df['MW'] = df['Mol'].apply(find_mw_in_mol)

    # Remove unncessary RDKit 'Mol' column
    df = df.drop(columns=['Mol'])

    return df
