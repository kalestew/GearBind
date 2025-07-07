#!/usr/bin/env python3
"""
generate_cr3022_mutations.py
Generate all single-point mutations for CR3022 CDRs and create mutant structures using FoldX
"""

import os
import pandas as pd
import subprocess
from pathlib import Path
from tqdm import tqdm
import multiprocessing as mp
from functools import partial

# Configuration
FOLDX_PATH = "/path/to/foldx"  # Update this to your FoldX path
PDB_FILE = "6xc3_wt.pdb"  # Wild-type CR3022-RBD complex
OUTPUT_DIR = "cr3022_mutants"
DATA_CSV = "cr3022_mutations.csv"

# CDR regions from the paper
CDR_REGIONS = {
    'H': [  # Heavy chain
        list(range(26, 36)),   # 26-35
        list(range(50, 67)),   # 50-66  
        list(range(99, 109))   # 99-108
    ],
    'L': [  # Light chain
        list(range(24, 41)),   # 24-40
        list(range(56, 63)),   # 56-62
        list(range(95, 104))   # 95-103
    ]
}

# Standard amino acids (excluding the wild-type will be handled per position)
AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def get_residue_at_position(pdb_file, chain, position):
    """Extract the wild-type residue at a specific position."""
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[21] == chain:
                res_num = int(line[22:26].strip())
                if res_num == position:
                    res_name = line[17:20].strip()
                    # Convert 3-letter to 1-letter code
                    three_to_one = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                    }
                    return three_to_one.get(res_name, 'X')
    return None

def generate_mutations_list():
    """Generate all single-point mutations for CDR regions."""
    mutations = []
    
    for chain, regions in CDR_REGIONS.items():
        for region in regions:
            for position in region:
                # Get wild-type residue
                wt_residue = get_residue_at_position(PDB_FILE, chain, position)
                if wt_residue:
                    # Generate all 20 amino acid mutations (including self-mutation)
                    for mut_residue in AMINO_ACIDS:
                        mutation = f"{wt_residue}{chain}{position}{mut_residue}"
                        mutations.append({
                            'mutation': mutation,
                            'chain': chain,
                            'position': position,
                            'wt_residue': wt_residue,
                            'mut_residue': mut_residue,
                            'is_wt': wt_residue == mut_residue
                        })
    
    print(f"Generated {len(mutations)} mutations")
    return mutations

def prepare_foldx_files(pdb_file):
    """Prepare PDB file for FoldX (RepairPDB)."""
    repaired_pdb = pdb_file.replace('.pdb', '_repaired.pdb')
    
    # Create FoldX repair command
    cmd = [
        FOLDX_PATH,
        "--command=RepairPDB",
        f"--pdb={pdb_file}",
        "--output-dir=.",
        "--clean-mode=3"
    ]
    
    print("Running RepairPDB...")
    subprocess.run(cmd, check=True)
    
    # FoldX outputs with specific naming convention
    expected_repair = pdb_file.replace('.pdb', '_Repair.pdb')
    if os.path.exists(expected_repair):
        os.rename(expected_repair, repaired_pdb)
    
    return repaired_pdb

def run_foldx_buildmodel(repaired_pdb, mutation, output_name):
    """Run FoldX BuildModel for a single mutation."""
    # Create individual mutation file
    mut_file = f"individual_list_{mutation}.txt"
    with open(mut_file, 'w') as f:
        f.write(f"{mutation};\n")
    
    # Run BuildModel
    cmd = [
        FOLDX_PATH,
        "--command=BuildModel",
        f"--pdb={repaired_pdb}",
        f"--mutant-file={mut_file}",
        "--output-dir=.",
        "--numberOfRuns=1",
        "--clean-mode=3"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Find and rename output file
        base_name = repaired_pdb.replace('.pdb', '')
        mutant_file = f"{base_name}_1.pdb"
        
        if os.path.exists(mutant_file):
            final_path = os.path.join(OUTPUT_DIR, output_name)
            os.rename(mutant_file, final_path)
            success = True
        else:
            success = False
            
    except subprocess.CalledProcessError:
        success = False
    
    # Cleanup
    if os.path.exists(mut_file):
        os.remove(mut_file)
    
    return success

def process_mutation(args):
    """Process a single mutation (for parallel processing)."""
    mutation_data, repaired_pdb = args
    mutation = mutation_data['mutation']
    
    # Create output filename
    output_name = f"cr3022_{mutation}.pdb"
    output_path = os.path.join(OUTPUT_DIR, output_name)
    
    # Skip if already exists
    if os.path.exists(output_path):
        return mutation_data, output_name, True
    
    # For wild-type mutations (self-mutations), just copy the file
    if mutation_data['is_wt']:
        import shutil
        shutil.copy(repaired_pdb, output_path)
        return mutation_data, output_name, True
    
    # Run FoldX
    success = run_foldx_buildmodel(repaired_pdb, mutation, output_name)
    
    return mutation_data, output_name, success

def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Generate mutations list
    print("Generating mutations list...")
    mutations = generate_mutations_list()
    
    # Prepare PDB file
    print("Preparing PDB file with FoldX RepairPDB...")
    repaired_pdb = prepare_foldx_files(PDB_FILE)
    
    # Process mutations in parallel
    print(f"Processing {len(mutations)} mutations...")
    
    # Prepare arguments for parallel processing
    args = [(mut, repaired_pdb) for mut in mutations]
    
    # Use multiprocessing
    with mp.Pool(mp.cpu_count()) as pool:
        results = list(tqdm(
            pool.imap(process_mutation, args),
            total=len(mutations),
            desc="Generating mutant structures"
        ))
    
    # Create data.csv file for GearBind
    data_rows = []
    successful = 0
    
    for mutation_data, output_name, success in results:
        if success:
            successful += 1
            data_rows.append({
                'pdb_id': 'cr3022',
                'mutation': mutation_data['mutation'],
                'chain_a': 'HL',  # Antibody chains
                'chain_b': 'C',   # Antigen chain
                'wt_protein': os.path.basename(repaired_pdb),
                'mt_protein': output_name
            })
    
    # Save data.csv
    df = pd.DataFrame(data_rows)
    df.to_csv(DATA_CSV, index=False)
    
    print(f"\nCompleted! Successfully generated {successful}/{len(mutations)} mutant structures")
    print(f"Data saved to {DATA_CSV}")
    
    # Also save detailed mutation info
    detailed_df = pd.DataFrame(mutations)
    detailed_df.to_csv('cr3022_mutations_detailed.csv', index=False)

if __name__ == "__main__":
    main()