from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit import DataStructs
from PIL import Image

functional_group_smiles_dict = {
    # Hydrocarbons
    "Alkene": Chem.MolFromSmarts("C=C"),
    "Alkyne": Chem.MolFromSmarts("C#C"),
    "Aromatic ring": Chem.MolFromSmarts("a1aaaaa1"),

    # Oxygen-containing groups
    "Alcohol": Chem.MolFromSmarts("[CX4][OX2H]"),                       # R–OH (not part of carbonyl)
    "Phenol": Chem.MolFromSmarts("c[OH]"),                              # Ar–OH
    "Ether": Chem.MolFromSmarts("[OD2]([#6])[#6]"),                     # R–O–R'
    "Aldehyde": Chem.MolFromSmarts("[CX3H1](=O)[#6]"),                  # R–CHO
    "Ketone": Chem.MolFromSmarts("[CX3](=O)[#6]"),                      # R–CO–R'
    "Carboxylic acid": Chem.MolFromSmarts("C(=O)[OX2H1]"),              # R–COOH
    "Ester": Chem.MolFromSmarts("C(=O)O[#6]"),                          # R–COOR'
    "Acid anhydride": Chem.MolFromSmarts("C(=O)OC(=O)"),                # (RCO)₂O
    "Acid chloride": Chem.MolFromSmarts("C(=O)Cl"),                     # RCOCl

    # Nitrogen-containing groups
    "Amine (primary)": Chem.MolFromSmarts("[NX3;H2][CX4]"),             # R–NH₂
    "Amine (secondary)": Chem.MolFromSmarts("[NX3;H1]([CX4])[CX4]"),    # R₂NH
    "Amine (tertiary)": Chem.MolFromSmarts("[NX3]([CX4])([CX4])[CX4]"), # R₃N
    "Amide": Chem.MolFromSmarts("C(=O)N"),                              # R–CONH₂ / R–CONR₂
    "Nitrile": Chem.MolFromSmarts("C#N"),                               # R–C≡N
    "Nitro": Chem.MolFromSmarts("[NX3](=O)=O"),                         # R–NO₂

    # Sulfur-containing groups
    "Thiol": Chem.MolFromSmarts("[#16X2H]"),                            # R–SH
    "Thioether": Chem.MolFromSmarts("[#16X2]([#6])[#6]"),               # R–S–R'
    "Sulfonic acid": Chem.MolFromSmarts("S(=O)(=O)[OH]"),               # R–SO₃H

    # Halogens
    "Alkyl halide": Chem.MolFromSmarts("[CX4][F,Cl,Br,I]"),             # R–X
    "Aryl halide": Chem.MolFromSmarts("c[F,Cl,Br,I]"),                  # Ar–X

    # Special cases
    "Quaternary ammonium": Chem.MolFromSmarts("[NX4+]"),                # R₄N⁺
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

def visualize_mol(mol):
    '''
    Converts mol object into PIL image object
    Parameter:
        mol : Chem.Mol = RDKit mol object
    Return:
        PIL Image : Png = Visual image of molecule
    '''
    img = Draw.MolToImage(mol)
    img.show()
    img.save('molecule.png', format='PNG')
    return img

def check_fg_in_mol(mol):
    '''
    Checks whether functional groups listed in functional_group_smiles_dict present in SMILES
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
