import sys
from pathlib import Path
import os

# Add the script directory to Python path
script_dir = Path(__file__).parent
sys.path.append(str(script_dir))

# Ensure descriptor directory exists
descriptor_dir = Path(script_dir).parent / 'descriptor'
descriptor_dir.mkdir(exist_ok=True)

from generate_descriptors import (
    create_nickel_catalyst_descriptors,
    create_photoredox_descriptors,
    create_micellar_descriptors,
    create_reaction_condition_descriptors
)
from substrate_descriptors import create_substrate_descriptors

def generate_all_descriptors():
    """Generate all descriptor files for the dual nickel/photoredox catalysis system."""
    print("Generating descriptor files...")
    
    # Create nickel catalyst descriptors
    print("Creating nickel catalyst descriptors...")
    create_nickel_catalyst_descriptors()
    
    # Create photoredox catalyst descriptors
    print("Creating photoredox catalyst descriptors...")
    create_photoredox_descriptors()
    
    # Create micellar media descriptors
    print("Creating micellar media descriptors...")
    create_micellar_descriptors()
    
    # Create reaction condition descriptors
    print("Creating reaction condition descriptors...")
    create_reaction_condition_descriptors()
    
    # Create substrate descriptors
    print("Creating substrate descriptors...")
    create_substrate_descriptors()
    
    print("All descriptor files have been generated successfully!")
    
    # List generated files
    print("\nGenerated descriptor files:")
    for file in descriptor_dir.glob('*.xlsx'):
        print(f"- {file.name}")

if __name__ == '__main__':
    generate_all_descriptors() 