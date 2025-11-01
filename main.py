from utils.rdkit_utils import visualize_mol, check_fg_in_mol, find_mw_in_mol

from utils.input_utils import input_smiles, input_csv
from utils.output_utils import export_to_csv, export_to_excel


def main():
    # Potential error/exception handling
    while True: 
        try:
            choice = int(input('Choose the format you would like to insert (1: SMILES format, 2: CSV file): '))

            if choice in [1,2]:
                print(f'You chose {choice}')
                break
            else:
                print(f'{choice} is an invalid choice. Please try again.')
        except ValueError:
            print('Invalid entry. Please insert number 1 or 2.')

    # If user input is SMILES format
    if choice == 1:
        mol = input_smiles()
        print(f'\nMolecule contains the following functional groups: {check_fg_in_mol(mol)}, \n\nMolecular weight of the molecule is: {find_mw_in_mol(mol)}, \n\nImage has been saved in the file directory.\n\n{visualize_mol(mol)}')
    
    # If user input is CSV file
    elif choice == 2: 
        df = input_csv()
        print(f'Here are the results: \n{df}') ### For display directly on website
        
        # User choice of file export (csv, excel)
        while True:
            try:
                file_choice = int(input('Choose the file format you would like to export (1: CSV, 2: Excel): '))
                if file_choice == 1:
                    export_to_csv(df)
                    break
                elif file_choice == 2:
                    export_to_excel(df)
                    break
                else:
                    print('Invalid entry. Please insert a number between 1-3')
            except ValueError:
                print('Invalid entry. Please insert a number between 1-3.')



if __name__ == "__main__":
    main()
