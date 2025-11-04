from utils.rdkit_utils import visualize_mol, check_fg_in_mol, find_mw_in_mol
from utils.input_utils import input_smiles, input_csv
from utils.output_utils import export_to_csv, export_to_excel
from flask import Flask, render_template, request, url_for, redirect, session
import os

app = Flask(__name__)
app.secret_key = 'ChemDetective2025'

@app.route('/', methods=['POST', 'GET'])
def homepage():
    smiles = None
    
    if request.method == 'POST':
        smiles = request.form.get('smiles')

        if smiles:
            session['smiles'] = smiles
            return redirect(url_for('result_smiles'))

    return render_template('index.html', smiles=smiles)


@app.route('/result-smiles', methods=['POST', 'GET'])
def result_smiles():
    mol = None
    matched_fgs = None
    MW = None
    img_path = None
    smiles = None

    smiles = session.get('smiles')
    return_home = request.form.get('return_home')
    
    if not smiles:
        return redirect(url_for('homepage'))
    
    if return_home:
        return redirect(url_for('homepage'))

    mol = input_smiles(smiles)

    if mol:
        matched_fgs = check_fg_in_mol(mol) # Functional groups present 
        MW = find_mw_in_mol(mol) # Molecular weight
        img_path = visualize_mol(mol) # Molecule image
    else:
        return redirect(url_for('homepage'))

    return render_template(
        'result_smiles.html',
        matched_fgs=matched_fgs,
        MW=MW,
        img_path=img_path,
        smiles=smiles)


@app.route('/result-csv', methods=['GET', 'POST'])
def result_csv():
    file_path = session.get('file_path')

    if not file_path or not os.path.exists(file_path):
        return redirect(url_for('homepage'))

    result_df = input_csv(file_path)
    return render_template('result_csv.html', table=result_df.to_html())


@app.route('/upload', methods=['POST'])
def upload():
    file = request.files.get('file')

    if not file:
        return redirect(url_for('homepage'))
    
    uploads_folder = os.path.join(app.root_path, 'uploads') # defines file download location
    os.makedirs(uploads_folder, exist_ok=True) # creates new if no 'uploads' folder
    file_path = os.path.join(uploads_folder, file.filename) # defines filepath of inserted file
    file.save(file_path) # saves file into defined filepath

    session['file_path'] = file_path
    return redirect(url_for('result_csv'))


if __name__ == '__main__':
    app.run(debug=True)




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
