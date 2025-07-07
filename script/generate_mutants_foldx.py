#!/usr/bin/env python3
"""
script/generate_mutants_foldx.py
--------------------------------
Generate single-point mutant structures for an antibody–antigen complex
using FoldX.

Features
========
1. Single command to perform full saturation mutagenesis on user-defined
   residue positions (e.g. CDRs) of the antibody chains.
2. Runs FoldX RepairPDB once, then FoldX BuildModel for every mutation.
3. Generates a GearBind-compatible data.csv describing WT / mutant pairs.
4. Multi-processing support for faster mutant generation.
5. Robust restart logic (skips mutants whose structure file already exists).

Typical usage
=============

    python script/generate_mutants_foldx.py \\
        --pdb 6xc3_wt.pdb \\
        --chain-a HL \\
        --chain-b C \\
        --regions-file config/cr3022_cdr_regions.yaml \\
        --foldx-path /opt/FoldX/foldx \\
        --output-dir data/cr3022_mutants

The *regions-file* is a YAML mapping `{ H: [[26, 35], [50, 66], ...], L: [...] }`.
If *regions-file* is omitted you can pass `--regions` directly, e.g.

    --regions "H:26-35,50-66,99-108;L:24-40,56-62,95-103"

The script produces

    <output-dir>/repaired.pdb           # FoldX-repaired WT complex
    <output-dir>/mutants/<mutation>.pdb # Each mutant structure
    <output-dir>/data.csv               # GearBind input table

where <mutation> follows FoldX 1-letter convention, e.g. *SH26A*.

Author: GearBind helper script
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd  # type: ignore
import yaml

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AMINO_ACIDS = list(THREE_TO_ONE.values())


def parse_region_string(s: str) -> Dict[str, List[range]]:
    """Parse a CLI region string like ``"H:26-35,50-66;L:24-40"``."""
    chains = {}
    chain_segs = [seg for seg in re.split(r";|\s", s) if seg]
    for seg in chain_segs:
        chain, ranges = seg.split(":")
        parsed = []
        for r in ranges.split(","):
            if not r:
                continue
            start, end = (int(x) for x in r.split("-"))
            parsed.append(range(start, end + 1))
        chains[chain] = parsed
    return chains


def extract_wt_residue(pdb_path: Path, chain: str, res_id: int) -> str:
    """Return the WT 1-letter residue code for *chain* and *res_id*."""
    with pdb_path.open() as fh:
        for line in fh:
            if line.startswith("ATOM") and line[21] == chain:
                rid = int(line[22:26])
                if rid == res_id:
                    three = line[17:20].strip()
                    return THREE_TO_ONE.get(three, "X")
    return "X"  # Unknown / non-standard


# ---------------------------------------------------------------------------
# FoldX helpers
# ---------------------------------------------------------------------------

def run(cmd: Sequence[str], **kwargs):
    """Run *cmd* via subprocess.run."""
    print(" [36m$", " ".join(cmd), "[0m")
    return subprocess.run(cmd, **kwargs)


def foldx_repair(pdb: Path, foldx: Path, workdir: Path) -> Path:
    repaired = workdir / (pdb.stem + "_repaired.pdb")
    if repaired.exists():
        return repaired
    
    # Copy PDB to workdir for FoldX
    local_pdb = workdir / pdb.name
    shutil.copy(pdb, local_pdb)
    
    # Run FoldX in the workdir
    result = run([str(foldx), "--command=RepairPDB", f"--pdb={pdb.name}", "--clean-mode=3"], 
                 cwd=str(workdir), capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"\n[ERROR] FoldX RepairPDB failed with exit code {result.returncode}")
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        raise RuntimeError(f"FoldX RepairPDB failed for {pdb.name}")
    
    produced = workdir / (pdb.stem + "_Repair.pdb")
    if not produced.exists():
        print(f"\n[ERROR] Expected repaired PDB file not found: {produced}")
        print(f"Available files in {workdir}: {list(workdir.glob('*.pdb'))}")
        raise RuntimeError(f"FoldX RepairPDB did not produce expected output file")
    
    produced.rename(repaired)
    
    # Clean up the copied PDB
    local_pdb.unlink(missing_ok=True)
    
    return repaired


def build_mutant(repaired: Path, mutation: str, foldx: Path, out_pdb: Path):
    """Build a mutant structure using FoldX BuildModel."""
    if out_pdb.exists():
        return True
    
    workdir = repaired.parent
    mut_file = workdir / f"individual_list_{mutation}.txt"
    mut_file.write_text(f"{mutation};\n")
    
    try:
        # Run FoldX BuildModel
        result = subprocess.run([str(foldx), "--command=BuildModel", 
             f"--pdb={repaired.name}", 
             f"--mutant-file={mut_file.name}", 
             "--numberOfRuns=1", 
             "--clean-mode=3"], 
            capture_output=True, text=True, cwd=str(workdir))
        
        # Check if FoldX failed
        if result.returncode != 0:
            print(f"\n[ERROR] FoldX failed for mutation {mutation} with exit code {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")
            return False
        
        # FoldX creates output files with different patterns depending on version
        # Common patterns include:
        # - <basename>_1.pdb (in the main directory)
        # - <basename>_1_0.pdb
        # - <basename>/<basename>_1.pdb (in a subdirectory)
        
        # Try different possible output locations
        possible_files = [
            workdir / f"{repaired.stem}_1.pdb",
            workdir / f"{repaired.stem}_1_0.pdb",
            workdir / repaired.stem / f"{repaired.stem}_1.pdb",
            workdir / repaired.stem / f"{repaired.stem}_1_0.pdb",
        ]
        
        # Also check for WT_ prefix (FoldX sometimes adds this)
        possible_files.extend([
            workdir / f"WT_{repaired.stem}_1.pdb",
            workdir / f"WT_{repaired.stem}_1_0.pdb",
        ])
        
        # Find the first existing file
        for candidate in possible_files:
            if candidate.exists():
                # Copy to the desired output location
                shutil.copy(candidate, out_pdb)
                # Clean up the original
                candidate.unlink()
                return True
        
        # If we couldn't find the output file, log available files for debugging
        print(f"\n[Warning] Could not find output for mutation {mutation}")
        print(f"Expected one of: {[str(p) for p in possible_files]}")
        print(f"Available PDB files in {workdir}:")
        for pdb_file in sorted(workdir.glob('*.pdb')):
            print(f"  - {pdb_file.name}")
        # Also check subdirectories
        if (workdir / repaired.stem).exists():
            print(f"Available PDB files in {workdir / repaired.stem}:")
            for pdb_file in sorted((workdir / repaired.stem).glob('*.pdb')):
                print(f"  - {pdb_file.name}")
        
        return False
        
    except Exception as e:
        print(f"\n[ERROR] Exception while building mutation {mutation}: {type(e).__name__}: {e}")
        return False
        
    finally:
        # Clean up
        mut_file.unlink(missing_ok=True)
        # Clean up FoldX temporary directory if it exists
        temp_dir = workdir / repaired.stem
        if temp_dir.exists() and temp_dir.is_dir():
            shutil.rmtree(temp_dir, ignore_errors=True)
        # Clean up other FoldX output files
        for pattern in ["*.fxout", "WT_*.pdb", f"{repaired.stem}_*.pdb", "individual_list_*.txt"]:
            for f in workdir.glob(pattern):
                if f != repaired and f != out_pdb:  # Don't delete the repaired structure or successful output
                    f.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def build_mutations_table(pdb: Path, regions: Dict[str, List[range]], skip_self: bool = True) -> List[Dict]:
    """Return a list of mutation dicts (wt_res, mut_res, chain, pos, label)."""
    table = []
    for chain, segs in regions.items():
        for seg in segs:
            for pos in seg:
                wt = extract_wt_residue(pdb, chain, pos)
                for mut in AMINO_ACIDS:
                    if skip_self and wt == mut:
                        continue  # Skip self-mutations
                    label = f"{wt}{chain}{pos}{mut}"
                    table.append({"chain": chain, "position": pos, "wt": wt, "mut": mut, "mutation": label})
    return table


def worker(args):
    mutation, repaired, foldx, out_dir = args
    try:
        ok = build_mutant(repaired, mutation["mutation"], foldx, out_dir / f"{mutation['mutation']}.pdb")
        return mutation, ok, None
    except Exception as e:
        print(f"\n[ERROR] Worker failed for mutation {mutation['mutation']}: {type(e).__name__}: {e}")
        return mutation, False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Generate saturation mutagenesis mutant structures with FoldX for GearBind")
    parser.add_argument("--pdb", required=True, type=Path, help="Path to WT PDB file (antibody-antigen complex)")
    parser.add_argument("--chain-a", required=True, help="Interacting antibody chains, e.g. HL")
    parser.add_argument("--chain-b", required=True, help="Antigen chain, e.g. C")
    
    # Default FoldX path for this repository
    default_foldx = Path(__file__).parent.parent / "foldx" / "foldx_51binaryMac"
    if not default_foldx.exists():
        default_foldx = "foldx"  # Fallback to PATH
    
    parser.add_argument("--foldx-path", default=str(default_foldx), type=Path, 
                       help=f"Path to FoldX binary (default: {default_foldx})")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--regions-file", type=Path, help="YAML file defining regions to mutate {H: [[26,35], ...], L: [...]}")
    group.add_argument("--regions", type=str, help="Region string, e.g. 'H:26-35,50-66;L:24-40,56-62' ")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to store mutants and CSV")
    parser.add_argument("--workers", type=int, default=max(mp.cpu_count() - 1, 1), help="Parallel workers (default: all cores minus 1)")
    parser.add_argument("--include-self", action="store_true", help="Include self-mutations (e.g., I->I)")
    parser.add_argument("--sequential", action="store_true", help="Run mutations sequentially instead of in parallel")

    args = parser.parse_args()

    # Check FoldX executable
    if not args.foldx_path.exists():
        print(f"\n[ERROR] FoldX binary not found at: {args.foldx_path}")
        print("Please provide the correct path using --foldx-path")
        return
    
    if not os.access(args.foldx_path, os.X_OK):
        print(f"\n[ERROR] FoldX binary at {args.foldx_path} is not executable!")
        print(f"Please make it executable with: chmod +x {args.foldx_path}")
        if os.name == 'posix' and os.uname().sysname == 'Darwin':  # macOS
            print("\nOn macOS, you may also need to allow it in System Preferences:")
            print("1. Try running it once: ./foldx_51binaryMac")
            print("2. Go to System Preferences > Security & Privacy")
            print("3. Click 'Allow Anyway' for the blocked app")
            print("4. Try running it again and click 'Open' when prompted")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mutants_dir = args.output_dir / "mutants"
    mutants_dir.mkdir(exist_ok=True)

    # Load regions
    if args.regions_file:
        regions = yaml.safe_load(args.regions_file.read_text())
    else:
        regions = parse_region_string(args.regions)

    # Repair WT structure
    repaired_pdb = foldx_repair(args.pdb, args.foldx_path, args.output_dir)

    # Enumerate mutations
    mutation_table = build_mutations_table(args.pdb, regions, skip_self=not args.include_self)
    print(f"Total mutations to generate: {len(mutation_table)}")
    
    if args.sequential:
        # Sequential processing
        print("Running mutations sequentially...")
        results = []
        for i, m in enumerate(mutation_table):
            print(f"Processing mutation {i+1}/{len(mutation_table)}: {m['mutation']}")
            result = worker((m, repaired_pdb, args.foldx_path, mutants_dir))
            results.append(result)
    else:
        # Parallel processing with fallback
        print(f"Running mutations in parallel with {args.workers} workers...")
        pool_args = [ (m, repaired_pdb, args.foldx_path, mutants_dir) for m in mutation_table ]
        
        try:
            with mp.Pool(args.workers) as pool:
                results = pool.map(worker, pool_args)
        except Exception as e:
            print(f"\n[WARNING] Parallel processing failed: {type(e).__name__}: {e}")
            print("Falling back to sequential processing...")
            results = []
            for i, args_tuple in enumerate(pool_args):
                print(f"Processing mutation {i+1}/{len(pool_args)}: {args_tuple[0]['mutation']}")
                result = worker(args_tuple)
                results.append(result)

    # Collect successes
    rows = []
    success = 0
    failed = []
    for (m, ok, err) in results:
        if ok:
            success += 1
            rows.append({
                "pdb_id": args.pdb.stem,
                "mutation": m["mutation"],
                "chain_a": args.chain_a,
                "chain_b": args.chain_b,
                "wt_protein": repaired_pdb.name,
                "mt_protein": f"mutants/{m['mutation']}.pdb"
            })
        else:
            failed.append(m["mutation"])
    
    print(f"\nSuccessfully built {success}/{len(mutation_table)} mutant structures")
    if failed:
        print(f"Failed mutations ({len(failed)}): {', '.join(failed[:10])}" + 
              (" ..." if len(failed) > 10 else ""))
    
    if not rows:
        print("\n[ERROR] No mutant structures were successfully generated!")
        print("Check the FoldX error messages above for details.")
        return

    pd.DataFrame(rows).to_csv(args.output_dir / "data.csv", index=False)
    print("GearBind data.csv written to", args.output_dir / "data.csv")


if __name__ == "__main__":
    main() 