# ChemDetective — README

**Version:** 1.0

ChemDetective is a compact, Python-based tool for quickly analysing small organic molecules. It accepts either a single SMILES string or a CSV file of molecules, and returns identified functional groups, molecular weight, and a rendered image of the molecule. The purpose of this project was primarily to apply libraries: Pandas, RDKit, Flask. 


## Features

* Accepts input as:

  * A single molecule in **SMILES** format, or
  * A CSV file with two columns: **Molecule Name**, **SMILES**.
* <ins>For a single SMILES input</ins>: returns identified **Functional Groups**, **Molecular Weight**, and an **Image of the Molecule**.
* <ins>For a CSV input</ins>: returns the original table extended with two columns — **Functional Groups** and **Molecular Weight** — and allows downloading the results as CSV or Excel.


## Installation

1) Download/Clone all files in this repo.
2) Install all python libraries listed in requirements.txt
3) Run code.
4) Open your browser and navigate to the address (eg. `http://127.0.0.1:5000`) shown in the terminal to use the ChemDetective UI.


## Design notes — functional group detection

ChemDetective currently detects functional groups by matching RDKit **SMARTS** patterns stored in a dictionary. The detection pipeline works like this:

1. Each functional group is defined by one SMARTS pattern.
2. RDKit is used to check whether each SMARTS matches the molecule.
3. A `conflict_rules_dict` is applied to avoid reporting less-specific groups when a more specific (dominant) group is present. For example, incorrectly detecting alcohol (-OH) group in carboxylic acid (-COOH) group. 

### Known limitation

* The present approach can produce false negatives when a `dominant` and a `filtered` functional group appear in *different parts* of the same molecule. The current `conflict_rules_dict` removes the filtered group entirely when the dominant group is present anywhere in the molecule — even if they are on distinct substructures — causing incorrect suppression of true groups.

## Demonstration

### Input organic molecule in SMILES format:

![input_smiles](input_smiles.gif)

### Input CSV file:

![input_csv](input_csv.gif)
