#!/usr/bin/env python3
"""
run_gearbind_predictions.py
Run GearBind predictions on all CR3022 mutants
"""

import os
import yaml
import pandas as pd
from pathlib import Path

# Paths
DATA_CSV = "cr3022_mutations.csv"
OUTPUT_DIR = "predictions"
CONFIG_TEMPLATE = "config/predict/CR3022_GearBindP.yaml"

def create_custom_dataset_class():
    """Create a custom dataset class for CR3022 mutations."""
    
    dataset_code = '''
# Add this to gearbind/dataset.py or create a new file

class CR3022Mutations(Dataset):
    """Dataset for CR3022 saturation mutagenesis."""
    
    def __init__(self, path, **kwargs):
        self.path = path
        self.data_dir = os.path.join(os.path.dirname(path), "cr3022_mutants")
        self.data = pd.read_csv(path)
        
        # Process each mutation
        data_list = []
        for _, row in self.data.iterrows():
            wt_path = os.path.join(self.data_dir, row["wt_protein"])
            mt_path = os.path.join(self.data_dir, row["mt_protein"])
            
            # Load structures
            wt_complex = data.Protein.from_pdb(wt_path)
            mt_complex = data.Protein.from_pdb(mt_path)
            
            # Extract chains
            chain_a = row["chain_a"]  # "HL"
            chain_b = row["chain_b"]  # "C"
            
            # Process wild-type
            wt_chains = wt_complex.split()
            wt_ab = sum([wt_chains[c] for c in chain_a])
            wt_ag = wt_chains[chain_b]
            
            # Process mutant
            mt_chains = mt_complex.split()
            mt_ab = sum([mt_chains[c] for c in chain_a])
            mt_ag = mt_chains[chain_b]
            
            item = {
                "pdb_id": row["pdb_id"],
                "mutation": row["mutation"],
                "wt_antibody": wt_ab,
                "wt_antigen": wt_ag,
                "mt_antibody": mt_ab,
                "mt_antigen": mt_ag,
            }
            data_list.append(item)
        
        self.data_list = data_list
        self.transform = kwargs.get("transform", None)
    
    def __len__(self):
        return len(self.data_list)
    
    def __getitem__(self, index):
        return self.data_list[index]
'''
    
    # Save to file
    with open('cr3022_mutations_dataset.py', 'w') as f:
        f.write(dataset_code)
    
    print("Created custom dataset class in cr3022_mutations_dataset.py")
    print("Add this to gearbind/dataset.py before running predictions")

def create_config_file():
    """Create configuration file for predictions."""
    
    config = {
        'output_dir': OUTPUT_DIR,
        'dataset': {
            'class': 'CR3022Mutations',
            'path': DATA_CSV,
            'transform': {
                'class': 'ProteinProteinInterfaceTransform',
                'max_radius': 15,
                'num_nearest': 128
            }
        },
        'task': {
            'class': 'PropertyPrediction',
            'model': {
                'class': 'GearBindModel',
                'path': 'checkpoints/GearBindP_downstream_30.pth'
            },
            'num_mlp_layer': 2,
            'graph_batch_size': 1
        },
        'optimizer': {
            'class': 'Adam',
            'lr': 0.0005,
            'weight_decay': 0
        },
        'engine': {
            'gpus': [0],
            'batch_size': 16,
            'log_interval': 50
        },
        'metric': ['mae', 'rmse', 'pearsonr', 'spearmanr']
    }
    
    # Save config
    with open('cr3022_predict_config.yaml', 'w') as f:
        yaml.dump(config, f)
    
    return 'cr3022_predict_config.yaml'

def run_predictions():
    """Run GearBind predictions."""
    
    config_file = create_config_file()
    
    # Run prediction script
    cmd = f"python script/predict.py -c {config_file}"
    print(f"Running: {cmd}")
    os.system(cmd)

def analyze_results():
    """Analyze prediction results and rank mutations."""
    
    # Read results
    results_file = f"{OUTPUT_DIR}/GearBindP_CR3022Mutations_test.csv"
    df = pd.read_csv(results_file)
    
    # Read original mutations for details
    mutations_df = pd.read_csv('cr3022_mutations_detailed.csv')
    
    # Merge results
    df = df.merge(mutations_df, on='mutation')
    
    # Calculate modified z-score
    median_ddg = df['prediction'].median()
    mad = (df['prediction'] - median_ddg).abs().median()
    df['modified_z_score'] = (df['prediction'] - median_ddg) / (1.4826 * mad)
    
    # Rank mutations (lower ddG is better for affinity)
    df = df.sort_values('prediction')
    
    # Save ranked results
    df.to_csv('cr3022_ranked_mutations.csv', index=False)
    
    # Print top mutations for each CDR
    print("\nTop affinity-improving mutations by CDR:")
    
    for chain in ['H', 'L']:
        print(f"\n{chain} chain CDRs:")
        chain_df = df[df['chain'] == chain]
        
        # Group by approximate CDR regions
        if chain == 'H':
            cdr_ranges = [(26, 35, 'CDR-H1'), (50, 66, 'CDR-H2'), (99, 108, 'CDR-H3')]
        else:
            cdr_ranges = [(24, 40, 'CDR-L1'), (56, 62, 'CDR-L2'), (95, 103, 'CDR-L3')]
        
        for start, end, cdr_name in cdr_ranges:
            cdr_df = chain_df[(chain_df['position'] >= start) & (chain_df['position'] <= end)]
            if not cdr_df.empty:
                print(f"\n  {cdr_name}:")
                for _, row in cdr_df.head(5).iterrows():
                    if not row['is_wt']:  # Skip wild-type
                        print(f"    {row['mutation']}: ddG = {row['prediction']:.3f}, z-score = {row['modified_z_score']:.3f}")

def main():
    # Step 1: Create custom dataset class
    create_custom_dataset_class()
    
    print("\n" + "="*50)
    print("Before proceeding:")
    print("1. Add the CR3022Mutations class from cr3022_mutations_dataset.py to gearbind/dataset.py")
    print("2. Make sure you have run the mutation generation script first")
    print("3. Ensure GearBind checkpoints are downloaded")
    print("="*50 + "\n")
    
    input("Press Enter when ready to continue...")
    
    # Step 2: Run predictions
    run_predictions()
    
    # Step 3: Analyze results
    analyze_results()

if __name__ == "__main__":
    main()