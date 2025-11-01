from utils.rdkit_utils import convert_to_mol_from_smiles, visualize_mol, check_fg_in_mol, find_mw_in_mol
from utils.chemical_database import csv_to_df, standardized_df, process_smiles_df


# Handling user input (smiles, csv file, image)
def input_smiles():
    while True:
        try: 
            smiles = input('Insert your molecule in SMILES format here: ')
            mol = convert_to_mol_from_smiles(smiles)
            if mol != None:
                return mol
            else:
                print(f'{smiles} is not a valid SMILES. Please check and try again.')
        except Exception as e:
            print(f'An error occurred: {e}. Please try again.')


def input_csv():
    return input('Insert your CSV file containing columns (Drug name, SMILES) here: ')

# def input_img():


# Handling file export (csv file, excel file, pdf file)
def export_to_csv(df):
    df.to_csv('output_file.csv', index=True)
    print('Results exported as output_file.csv file in the directory.')

def export_to_excel(df):
    df.to_excel('output_file.xlsx')
    print('Results exported as out_file.xlsx file in the directory.')

# def export_to_pdf(df):
    # with PdfPages('output_file.pdf') as pdf:
        # fig, ax = subplots(figsize=(10,6))
        # ax.axis('off')

        # table = ax.table(cellText  = df.values,
                         # colLabels = df.columns,
                         # loc       = 'center')
        
        # table.auto_set_font_size(False)
        # table.set_fontsize(4.5)

        # table.scale(1.2,1.5)

        # pdf.savefig()


def main():
    # Potential error/exception handling
    while True: 
        try:
            choice = int(input('Choose the format you would like to insert (1: SMILES format, 2: CSV file, 3: Image PNG file): '))

            if choice in [1,2,3]:
                print(f'You chose {choice}')
                break
            else:
                print(f'{choice} is an invalid choice. Please try again.')
        except ValueError:
            print('Invalid entry. Please insert a number between 1-3.')

    # If user input is SMILES format
    if choice == 1:
        mol = input_smiles()
        print(f'\nMolecule contains the following functional groups: {check_fg_in_mol(mol)}, \n\nMolecular weight of the molecule is: {find_mw_in_mol(mol)}, \n\nImage has been saved in the file directory.\n\n{visualize_mol(mol)}')
    
    # If user input is CSV file
    elif choice == 2: 
        file_path = input_csv()
        df = csv_to_df(file_path)
        df = standardized_df(df)
        df = process_smiles_df(df)
        print(f'Here are the results: \n{df}') ### For display directly on website
        
        # User choice of file export (csv, excel)
        while True:
            try:
                file_choice = int(input('Choose the file format you would like to export (1: CSV, 2: Excel, 3: PDF): '))
                if file_choice == 1:
                    export_to_csv(df)
                    break
                elif file_choice == 2:
                    export_to_excel(df)
                    break
                elif file_choice == 3:
                    export_to_pdf(df)
                    break
                else:
                    print('Invalid entry. Please insert a number between 1-3')
            except ValueError:
                print('Invalid entry. Please insert a number between 1-3.')

    # else:



if __name__ == "__main__":
    main()
