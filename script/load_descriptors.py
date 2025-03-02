import pandas as pd 
import numpy as np
import csv
import sys 
from os import sep
from script.load_dataset import chemical_space

def f_des_std(des_array):
    """Standardize descriptors by scaling to [0,1] range."""
    # Find columns with variation
    varying_cols = des_array.max(axis=0) != des_array.min(axis=0)
    react_feat_all = des_array[:, varying_cols]
    
    # Scale to [0,1] range
    react_feat_all = (react_feat_all - react_feat_all.min(axis=0)) / \
                    (react_feat_all.max(axis=0) - react_feat_all.min(axis=0))
    return react_feat_all

def get_reaction_descriptors(df):
    """Convert reaction conditions to descriptors."""
    # One-hot encode categorical variables
    electrode_df = pd.get_dummies(df['Anode/Cathode'], prefix='electrode')
    solvent_df = pd.get_dummies(df['Solvent'], prefix='solvent')
    electrolyte_df = pd.get_dummies(df['Electrolyte'], prefix='electrolyte')
    
    # Convert current/potential to numeric
    def extract_value(x):
        return float(x.split()[0])
    
    current_potential = df['Current/Potential'].apply(extract_value).values.reshape(-1, 1)
    
    # Combine all descriptors
    reaction_desc = np.hstack([
        electrode_df.values,
        solvent_df.values,
        electrolyte_df.values,
        current_potential
    ])
    
    return reaction_desc

def get_descriptors():
    """Load and combine all descriptors."""
    # Load reaction data
    df = pd.read_csv('./dataset/all_input_data.csv')
    
    # Get reaction condition descriptors
    reaction_desc = get_reaction_descriptors(df)
    
    # Load chemical descriptors
    nickel_df = pd.read_excel('descriptor/nickel_catalyst_descriptors.xlsx')
    photoredox_df = pd.read_excel('descriptor/photoredox_descriptors.xlsx')
    micellar_df = pd.read_excel('descriptor/micellar_descriptors.xlsx')
    substrate_df = pd.read_excel('descriptor/substrate_descriptors.xlsx')
    
    # Convert to numeric, dropping first column (names)
    nickel_desc = nickel_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
    photoredox_desc = photoredox_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
    micellar_desc = micellar_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
    substrate_desc = substrate_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').fillna(0).values
    
    # Repeat chemical descriptors for each reaction condition
    n_reactions = len(df)
    nickel_desc_rep = np.tile(nickel_desc[0], (n_reactions, 1))
    photoredox_desc_rep = np.tile(photoredox_desc[0], (n_reactions, 1))
    micellar_desc_rep = np.tile(micellar_desc[0], (n_reactions, 1))
    substrate_desc_rep = np.tile(substrate_desc[0], (n_reactions, 1))
    
    # Combine all descriptors
    all_desc = np.hstack([
        reaction_desc,
        nickel_desc_rep,
        photoredox_desc_rep,
        micellar_desc_rep,
        substrate_desc_rep
    ])
    
    # Standardize descriptors
    desc_std = f_des_std(all_desc)
    
    print(f"Descriptor shape: {desc_std.shape}")
    print(f"Number of samples: {n_reactions}")
    
    return desc_std