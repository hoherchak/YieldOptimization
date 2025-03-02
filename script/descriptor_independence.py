import pandas as pd
import numpy as np
from scipy.stats import pearsonr

def select_uncorrelated_features(df, threshold=0.9):
    """Select uncorrelated features based on correlation threshold."""
    # Calculate correlation matrix
    corr_matrix = df.corr().abs()
    
    # Create upper triangle mask
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    
    # Find features to drop
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]
    
    # Keep uncorrelated features
    df_uncorr = df.drop(columns=to_drop)
    
    return df_uncorr

def analyze_descriptor_independence():
    """Analyze descriptor independence and save results."""
    # Load descriptor files
    nickel_df = pd.read_excel('descriptor/nickel_catalyst_descriptors.xlsx')
    photoredox_df = pd.read_excel('descriptor/photoredox_descriptors.xlsx')
    micellar_df = pd.read_excel('descriptor/micellar_descriptors.xlsx')
    reaction_df = pd.read_excel('descriptor/reaction_condition_descriptors.xlsx')
    substrate_df = pd.read_excel('descriptor/substrate_descriptors.xlsx')
    
    # Store names
    names = {
        'nickel': nickel_df.iloc[:, 0],
        'photoredox': photoredox_df.iloc[:, 0],
        'micellar': micellar_df.iloc[:, 0],
        'reaction': reaction_df.iloc[:, 0],
        'substrate': substrate_df.iloc[:, 0]
    }
    
    # Convert to numeric, dropping first column (names)
    dfs = {
        'nickel': nickel_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce'),
        'photoredox': photoredox_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce'),
        'micellar': micellar_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce'),
        'reaction': reaction_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce'),
        'substrate': substrate_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    }
    
    # Print original dimensions
    for name, df in dfs.items():
        print(f"{name} descriptors: {df.shape[1]}")
    
    # Select uncorrelated features for each category
    uncorr_dfs = {}
    for name, df in dfs.items():
        # Drop columns with all NaN
        df = df.dropna(axis=1, how='all')
        # Fill remaining NaN with 0
        df = df.fillna(0)
        # Select uncorrelated features
        df_uncorr = select_uncorrelated_features(df)
        uncorr_dfs[name] = df_uncorr
        print(f"{name} uncorrelated descriptors: {df_uncorr.shape[1]}")
    
    # Combine all uncorrelated descriptors
    all_descriptors = pd.concat([df for df in uncorr_dfs.values()], axis=1)
    print(f"Total uncorrelated descriptors: {all_descriptors.shape[1]}")
    
    # Save uncorrelated descriptors
    all_descriptors.to_excel('descriptor/independent_descriptors.xlsx', index=False)
    
    # Save descriptor names
    names_df = pd.DataFrame({
        'Category': np.repeat(list(names.keys()), [len(v) for v in names.values()]),
        'Name': pd.concat(names.values())
    })
    names_df.to_excel('descriptor/descriptor_names.xlsx', index=False)
    
    return all_descriptors

if __name__ == '__main__':
    analyze_descriptor_independence()