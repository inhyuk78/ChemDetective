from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import os

functional_group_smiles_dict = {
    # Hydrocarbons
    "alkene": Chem.MolFromSmarts("C=C"),
    "alkyne": Chem.MolFromSmarts("C#C"),
    "aromatic ring": Chem.MolFromSmarts("a1aaaaa1"),

    # Oxygen-containing groups
    "alcohol": Chem.MolFromSmarts("[CX4][OX2H]"),                       # R–OH (not part of carbonyl)
    "phenol": Chem.MolFromSmarts("c[OH]"),                              # Ar–OH
    "ether": Chem.MolFromSmarts("[OD2]([#6])[#6]"),                     # R–O–R'
    "aldehyde": Chem.MolFromSmarts("[CX3H1](=O)[#6]"),                  # R–CHO
    "ketone": Chem.MolFromSmarts("[CX3](=O)[#6]"),                      # R–CO–R'
    "carboxylic acid": Chem.MolFromSmarts("C(=O)[OX2H1]"),              # R–COOH
    "ester": Chem.MolFromSmarts("C(=O)O[#6]"),                          # R–COOR'
    "acid anhydride": Chem.MolFromSmarts("C(=O)OC(=O)"),                # (RCO)₂O
    "acid chloride": Chem.MolFromSmarts("C(=O)Cl"),                     # RCOCl

    # Nitrogen-containing groups
    "amine (primary)": Chem.MolFromSmarts("[NX3;H2][CX4]"),             # R–NH₂
    "amine (secondary)": Chem.MolFromSmarts("[NX3;H1]([CX4])[CX4]"),    # R₂NH
    "amine (tertiary)": Chem.MolFromSmarts("[NX3]([CX4])([CX4])[CX4]"), # R₃N
    "amide": Chem.MolFromSmarts("C(=O)N"),                              # R–CONH₂ / R–CONR₂
    "nitrile": Chem.MolFromSmarts("C#N"),                               # R–C≡N
    "nitro": Chem.MolFromSmarts("[NX3](=O)=O"),                         # R–NO₂

    # Sulfur-containing groups
    "thiol": Chem.MolFromSmarts("[#16X2H]"),                            # R–SH
    "thioether": Chem.MolFromSmarts("[#16X2]([#6])[#6]"),               # R–S–R'
    "sulfonic acid": Chem.MolFromSmarts("S(=O)(=O)[OH]"),               # R–SO₃H

    # Halogens
    "alkyl halide": Chem.MolFromSmarts("[CX4][F,Cl,Br,I]"),             # R–X
    "aryl halide": Chem.MolFromSmarts("c[F,Cl,Br,I]"),                  # Ar–X

    # Special cases
    "quaternary ammonium": Chem.MolFromSmarts("[NX4+]"),                # R₄N⁺
}

conflict_rules_dict = {
    'carboxylic acid': ['alcohol', 'ketone', 'ether'],
    'ester': ['ketone', 'ether', 'alcohol'],
    'aldehyde': ['ketone'],
    'aromatic ring': ['alkene'],
    'amide': ['ketone', 'amine (primary)', 'amine (secondary)', 'amine (tertiary)'],
    'acid chloride': ['ketone'],
    'phenol': ['aromatic ring', 'alcohol'],
    'ether': ['alcohol'],
    'thioether': ['thiol']
}

def convert_to_mol_from_smiles(smiles):
    '''
    Converts a SMILES string into an RDKit Mol object
    Parameter:
        smiles : str = SMILES representation of a molecule
    Return:
        mol : Chem.Mol or None = RDKit Mol if successful, None is invalid SMILES
    '''
    mol = Chem.MolFromSmiles(smiles)
    return mol

def visualize_mol(mol, filename='molecule.png'):
    '''
    Converts mol object into image in PNG file
    Parameter:
        mol : Chem.Mol = RDKit mol object
        filename : str = filename for saving (default: molecule.png)
    Return:
        str: relative path to saved image (for flask to render)
    '''
    # Defining static folder path
    static_folder = os.path.join(os.getcwd(), 'static')
    os.makedirs(static_folder, exist_ok=True)
    img_path = os.path.join(static_folder, filename)

    img = Draw.MolToImage(mol)
    img.save(img_path, format='PNG')
    return f'{filename}'

def check_fg_in_mol(mol):
    '''
    Looks up functional_group_smiles_dict and conflict_rules_dict to identify FGs present
    Parameter:
        mol : Chem.Mol = RDKit mol object
    Return:
        if match: 'Molecule {smiles} contains the following functional groups: {', '.join(matched_fgs)}'
        if no match: 'Molecule {smiles} contains no listed functional groups'
    '''
    matched_fgs = []

    for name, pattern in functional_group_smiles_dict.items():     # iterates through each fg_mol, fg_name in functional_groups list
        if mol.HasSubstructMatch(pattern):             # If match, adds fg_name into empty matched_fgs list
            matched_fgs.append(name)

    for dominant_group, filtered_group in conflict_rules_dict.items():
        if dominant_group in matched_fgs:
            for group in filtered_group:
                if group in matched_fgs:
                    matched_fgs.remove(group)

    if matched_fgs:
        return f"{', '.join(matched_fgs)}"
    else:
        return f"Molecule contains no listed functional groups"

def find_mw_in_mol(mol):
    '''
    Calculate molecular weight of a molecule
    Parameter:
        mol : Chem.Mol = RDKit Mol object
    Return:
        mw : Float = Molecular weight in grams per mole (g/mol)
    '''
    mw = Descriptors.MolWt(mol)
    return f'{round(mw, 2)} g/mol'
