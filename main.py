from utils.rdkit_utils import visualize_mol, check_fg_in_mol, find_mw_in_mol
from utils.input_utils import input_smiles, input_csv
from utils.output_utils import export_to_csv, export_to_excel
from flask import Flask, render_template, request, url_for, redirect, session, flash
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
    
    if request.method == 'POST' and request.form.get('return_home'):
        return redirect(url_for('homepage'))

    if not smiles:
        flash('No SMILES input provided.', 'error')
        return redirect(url_for('homepage'))
    
    try:
        mol = input_smiles(smiles)
        matched_fgs = check_fg_in_mol(mol) # Functional groups present 
        MW = find_mw_in_mol(mol) # Molecular weight
        img_path = visualize_mol(mol) # Molecule image   
    except ValueError as e:
        flash(str(e), 'error')
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

    if request.method == 'POST' and request.form.get('return_home1'):
        return redirect(url_for('homepage'))
    
    if not file_path or not os.path.exists(file_path):
        flash('Uploaded file not found. Please re-upload your CSV file', 'error')
        return redirect(url_for('homepage'))
    
    try:
        result_df = input_csv(file_path)
    except FileNotFoundError:
        flash('File could not be found. Please try uploading again.', 'error')
        return redirect(url_for('homepage'))
    except ValueError as ve:
        flash(str(ve), 'error')
        return redirect(url_for('homepage'))
    
    if request.method == 'POST' and request.form.get('save_as_csv'):
        export_to_csv(result_df)
        flash('File saved successfully in your directory!', 'info')
    
    if request.method == 'POST' and request.form.get('save_as_excel'):
        export_to_excel(result_df)
        flash('File saved successfully in your directory!', 'info')

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