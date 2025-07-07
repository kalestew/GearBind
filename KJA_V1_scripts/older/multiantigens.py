#!/usr/bin/env python3
"""
process_multiple_antigens.py
Run the same CR3022 mutations against different antigens (BA.4, BA.1.1)
"""

import os
import shutil
import pandas as pd

ANTIGENS = {
    'WT': '6xc3_wt.pdb',
    'BA4': '6xc3_ba4.pdb',
    'BA11': '6xc3_ba11.pdb'
}

def process_antigen(antigen_name, antigen_pdb, mutations_csv):
    """Process mutations for a specific antigen variant."""
    
    print(f"\nProcessing {antigen_name} antigen...")
    
    # Create directories
    output_dir = f"cr3022_mutants_{antigen_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Update the first script to use this antigen
    # (You would modify the generate_cr3022_mutations.py script)
    
    # For demonstration, assume we have the mutant structures
    # and create a new data.csv for this antigen
    
    df = pd.read_csv(mutations_csv)
    df['pdb_id'] = f'cr3022_{antigen_name}'
    
    # Save new CSV
    new_csv = f"cr3022_mutations_{antigen_name}.csv"
    df.to_csv(new_csv, index=False)
    
    # Run GearBind predictions
    os.system(f"python run_gearbind_predictions.py --data {new_csv} --output predictions_{antigen_name}")
    
    return f"predictions_{antigen_name}/results.csv"

def compare_antigens():
    """Compare predictions across different antigens."""
    
    results = {}
    
    for antigen_name in ANTIGENS:
        results_file = f"predictions_{antigen_name}/results.csv"
        if os.path.exists(results_file):
            results[antigen_name] = pd.read_csv(results_file)
    
    # Merge results
    merged = None
    for antigen_name, df in results.items():
        df = df[['mutation', 'prediction']].rename(columns={'prediction': f'ddG_{antigen_name}'})
        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on='mutation')
    
    # Calculate average modified z-score across antigens
    for col in [c for c in merged.columns if c.startswith('ddG_')]:
        median_val = merged[col].median()
        mad = (merged[col] - median_val).abs().median()
        z_col = col.replace('ddG_', 'z_')
        merged[z_col] = (merged[col] - median_val) / (1.4826 * mad)
    
    # Average z-score
    z_cols = [c for c in merged.columns if c.startswith('z_')]
    merged['avg_z_score'] = merged[z_cols].mean(axis=1)
    
    # Rank by average z-score
    merged = merged.sort_values('avg_z_score')
    
    # Save results
    merged.to_csv('cr3022_mutations_all_antigens.csv', index=False)
    
    print("\nTop mutations across all antigens:")
    print(merged.head(20))

if __name__ == "__main__":
    # First run mutations against each antigen
    for antigen_name, antigen_pdb in ANTIGENS.items():
        process_antigen(antigen_name, antigen_pdb, "cr3022_mutations.csv")
    
    # Compare results
    compare_antigens()