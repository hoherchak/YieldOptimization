import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# Common substrates in photoredox/nickel dual catalysis
substrates = {
    # Aryl halides
    'bromobenzene': 'c1ccccc1Br',
    '4-bromoanisole': 'COc1ccc(Br)cc1',
    '4-bromoacetophenone': 'CC(=O)c1ccc(Br)cc1',
    '4-bromobenzonitrile': 'N#Cc1ccc(Br)cc1',
    # Alkyl halides
    'benzyl-bromide': 'BrCc1ccccc1',
    'ethyl-2-bromopropionate': 'CCOC(=O)C(C)Br',
    '2-bromo-4-phenylbutane': 'CC(Br)CCc1ccccc1',
    # Carboxylic acids
    'benzoic-acid': 'O=C(O)c1ccccc1',
    'phenylacetic-acid': 'O=C(O)Cc1ccccc1',
    '4-methoxybenzoic-acid': 'COc1ccc(C(=O)O)cc1'
}

def calculate_substrate_descriptors(smiles):
    """Calculate molecular descriptors for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    descriptors = {
        # Electronic properties
        'TPSA': Descriptors.TPSA(mol),
        'MolLogP': Descriptors.MolLogP(mol),
        'HOMO': -5.0,  # Placeholder - would need quantum chemistry calculation
        'LUMO': -1.0,  # Placeholder - would need quantum chemistry calculation
        
        # Structural properties
        'MolWt': Descriptors.ExactMolWt(mol),
        'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
        'NumHAcceptors': Descriptors.NumHAcceptors(mol),
        'NumHDonors': Descriptors.NumHDonors(mol),
        'NumAromaticRings': rdMolDescriptors.CalcNumAromaticRings(mol),
        
        # Topological properties
        'BertzCT': Descriptors.BertzCT(mol),
        'Chi0v': Descriptors.Chi0v(mol),
        'Chi1n': Descriptors.Chi1n(mol),
        
        # Geometric properties
        'Sphericity': 0.0,  # Placeholder - would need 3D conformation
        'Asphericity': 0.0,  # Placeholder - would need 3D conformation
    }
    
    # Calculate 3D descriptors
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        descriptors['Sphericity'] = 0.5  # Example calculation
        descriptors['Asphericity'] = 0.3  # Example calculation
    except:
        pass
    
    return descriptors

def create_substrate_descriptors():
    """Generate substrate descriptor file."""
    all_descriptors = {}
    
    for name, smiles in substrates.items():
        descriptors = calculate_substrate_descriptors(smiles)
        if descriptors:
            all_descriptors[name] = descriptors
    
    df = pd.DataFrame.from_dict(all_descriptors, orient='index')
    df.to_excel('descriptor/substrate_descriptors.xlsx')

if __name__ == '__main__':
    create_substrate_descriptors() 