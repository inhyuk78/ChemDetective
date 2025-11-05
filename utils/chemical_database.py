import pandas as pd
from utils.rdkit_utils import convert_to_mol_from_smiles, check_fg_in_mol,find_mw_in_mol

def csv_to_df(file_path):
    '''
    Converts CSV file into pandas dataframe
    Parameter:
        file_path : CSV = File (CSV) containing drug names, smiles columns
    Return:
        df : pd.DataFrame = DataFrame containing CSV file data if file found
    ''' 
    df = pd.read_csv(file_path)
    return df

def standardized_df(df): 
    '''
    Converts dataframe into standardized form for analysis
    Parameter:
        df : pd.DataFrame = Raw user input CSV file converted DataFrame
    Return:
        df : pd.DataFrame = Standardized DataFrame
    '''
    df.columns = [col.strip().lower() for col in df.columns]
    
    if len(df.columns) != 2:
        raise ValueError('File must contain exactly 2 columns (compound name, SMILES)')

    df = df.rename(columns={
        df.columns[0]: 'Drug_Name',
        df.columns[1]: 'SMILES'
    })

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

    # Add new column with MW
    df['MW'] = df['Mol'].apply(find_mw_in_mol)

    # Remove unncessary RDKit 'Mol' column
    df = df.drop(columns=['Mol'])

    return df