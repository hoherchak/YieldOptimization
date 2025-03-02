import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# Common Nickel catalysts and their properties
nickel_catalysts = {
    'Ni(COD)2': {
        'oxidation_state': 0,
        'electron_config': 'd10',
        'redox_potential': -1.2,  # V vs SCE
        'ligand_cone_angle': 180,  # degrees
        'buried_volume': 32.5,  # %Vbur
        'coordination_geometry': 'square planar',
        'molecular_weight': 275.91,
    },
    'NiCl2•glyme': {
        'oxidation_state': 2,
        'electron_config': 'd8',
        'redox_potential': -0.9,
        'ligand_cone_angle': 160,
        'buried_volume': 28.4,
        'coordination_geometry': 'octahedral',
        'molecular_weight': 230.59,
    },
    'Ni(dtbbpy)Br2': {
        'oxidation_state': 2,
        'electron_config': 'd8',
        'redox_potential': -1.1,
        'ligand_cone_angle': 175,
        'buried_volume': 35.2,
        'coordination_geometry': 'square planar',
        'molecular_weight': 492.08,
    }
}

# Common photoredox catalysts and their properties
photoredox_catalysts = {
    'Ir(ppy)3': {
        'absorption_max': 375,  # nm
        'emission_max': 518,  # nm
        'excited_state_lifetime': 1900,  # ns
        'quantum_yield': 0.40,
        'ground_state_potential': -2.19,  # V vs SCE
        'excited_state_potential': 0.31,  # V vs SCE
        'molecular_weight': 654.79,
    },
    'Ru(bpy)3Cl2': {
        'absorption_max': 452,
        'emission_max': 615,
        'excited_state_lifetime': 1100,
        'quantum_yield': 0.095,
        'ground_state_potential': -1.33,
        'excited_state_potential': 0.77,
        'molecular_weight': 748.63,
    },
    '4CzIPN': {
        'absorption_max': 435,
        'emission_max': 535,
        'excited_state_lifetime': 5100,
        'quantum_yield': 0.94,
        'ground_state_potential': -1.21,
        'excited_state_potential': 1.35,
        'molecular_weight': 788.93,
    }
}

# Common surfactants for micellar media
micellar_media = {
    'TPGS-750-M': {
        'cmc': 0.12,  # mM
        'hlb': 13.2,
        'aggregation_number': 45,
        'surface_tension': 41.2,  # mN/m
        'micelle_size': 52,  # nm
        'zeta_potential': -15.3,  # mV
        'polydispersity': 0.24,
    },
    'PS-750-M': {
        'cmc': 0.15,
        'hlb': 12.8,
        'aggregation_number': 42,
        'surface_tension': 39.8,
        'micelle_size': 48,
        'zeta_potential': -12.8,
        'polydispersity': 0.21,
    },
    'Cremophor EL': {
        'cmc': 0.09,
        'hlb': 14.1,
        'aggregation_number': 50,
        'surface_tension': 43.5,
        'micelle_size': 55,
        'zeta_potential': -18.2,
        'polydispersity': 0.28,
    }
}

# Generate descriptor files
def create_nickel_catalyst_descriptors():
    try:
        print("Creating nickel catalyst DataFrame...")
        df = pd.DataFrame.from_dict(nickel_catalysts, orient='index')
        print("Saving nickel catalyst descriptors...")
        df.to_excel('descriptor/nickel_catalyst_descriptors.xlsx')
        print("Nickel catalyst descriptors saved successfully!")
    except Exception as e:
        print(f"Error creating nickel catalyst descriptors: {str(e)}")

def create_photoredox_descriptors():
    try:
        print("Creating photoredox catalyst DataFrame...")
        df = pd.DataFrame.from_dict(photoredox_catalysts, orient='index')
        print("Saving photoredox descriptors...")
        df.to_excel('descriptor/photoredox_descriptors.xlsx')
        print("Photoredox descriptors saved successfully!")
    except Exception as e:
        print(f"Error creating photoredox descriptors: {str(e)}")

def create_micellar_descriptors():
    try:
        print("Creating micellar media DataFrame...")
        df = pd.DataFrame.from_dict(micellar_media, orient='index')
        print("Saving micellar descriptors...")
        df.to_excel('descriptor/micellar_descriptors.xlsx')
        print("Micellar descriptors saved successfully!")
    except Exception as e:
        print(f"Error creating micellar descriptors: {str(e)}")

def create_reaction_condition_descriptors():
    try:
        print("Creating reaction conditions DataFrame...")
        conditions = {
            'condition_set_1': {
                'temperature': 25,  # °C
                'light_intensity': 40,  # W/m2
                'wavelength': 450,  # nm
                'stirring_rate': 800,  # rpm
                'ni_loading': 2.0,  # mol%
                'pc_loading': 1.0,  # mol%
                'surfactant_conc': 2.0,  # wt%
                'water_content': 98,  # %
            },
            'condition_set_2': {
                'temperature': 30,
                'light_intensity': 50,
                'wavelength': 455,
                'stirring_rate': 1000,
                'ni_loading': 3.0,
                'pc_loading': 1.5,
                'surfactant_conc': 2.5,
                'water_content': 97,
            },
            'condition_set_3': {
                'temperature': 35,
                'light_intensity': 60,
                'wavelength': 460,
                'stirring_rate': 1200,
                'ni_loading': 4.0,
                'pc_loading': 2.0,
                'surfactant_conc': 3.0,
                'water_content': 96,
            }
        }
        df = pd.DataFrame.from_dict(conditions, orient='index')
        print("Saving reaction condition descriptors...")
        df.to_excel('descriptor/reaction_condition_descriptors.xlsx')
        print("Reaction condition descriptors saved successfully!")
    except Exception as e:
        print(f"Error creating reaction condition descriptors: {str(e)}")

if __name__ == '__main__':
    create_nickel_catalyst_descriptors()
    create_photoredox_descriptors()
    create_micellar_descriptors()
    create_reaction_condition_descriptors() 