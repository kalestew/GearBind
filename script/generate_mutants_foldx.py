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
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd
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
    chain_segs = [seg for seg in re.split(";|\s", s) if seg]
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
    """Run *cmd* via subprocess.run raising on error."""
    print(" [36m$", " ".join(cmd), "[0m")
    subprocess.run(cmd, check=True, **kwargs)


def foldx_repair(pdb: Path, foldx: Path, workdir: Path) -> Path:
    repaired = workdir / (pdb.stem + "_repaired.pdb")
    if repaired.exists():
        return repaired
    run([str(foldx), "--command=RepairPDB", f"--pdb={pdb}", "--output-dir", str(workdir), "--clean-mode=3"])
    produced = workdir / (pdb.stem + "_Repair.pdb")
    produced.rename(repaired)
    return repaired


def build_mutant(repaired: Path, mutation: str, foldx: Path, out_pdb: Path):
    if out_pdb.exists():
        return True
    mut_file = repaired.parent / f"individual_list_{mutation}.txt"
    mut_file.write_text(f"{mutation};\n")
    try:
        run([str(foldx), "--command=BuildModel", f"--pdb={repaired}", f"--mutant-file={mut_file}", "--output-dir", str(repaired.parent), "--numberOfRuns", "1", "--clean-mode", "3"], capture_output=True, text=True)
        candidate = repaired.with_suffix("").with_suffix("_1.pdb")
        if candidate.exists():
            candidate.rename(out_pdb)
            return True
        return False
    finally:
        mut_file.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def build_mutations_table(pdb: Path, regions: Dict[str, List[range]]) -> List[Dict]:
    """Return a list of mutation dicts (wt_res, mut_res, chain, pos, label)."""
    table = []
    for chain, segs in regions.items():
        for seg in segs:
            for pos in seg:
                wt = extract_wt_residue(pdb, chain, pos)
                for mut in AMINO_ACIDS:
                    label = f"{wt}{chain}{pos}{mut}"
                    table.append({"chain": chain, "position": pos, "wt": wt, "mut": mut, "mutation": label})
    return table


def worker(args):
    mutation, repaired, foldx, out_dir = args
    ok = build_mutant(repaired, mutation["mutation"], foldx, out_dir / f"{mutation['mutation']}.pdb")
    return mutation, ok


def main():
    parser = argparse.ArgumentParser(description="Generate saturation mutagenesis mutant structures with FoldX for GearBind")
    parser.add_argument("--pdb", required=True, type=Path, help="Path to WT PDB file (antibody-antigen complex)")
    parser.add_argument("--chain-a", required=True, help="Interacting antibody chains, e.g. HL")
    parser.add_argument("--chain-b", required=True, help="Antigen chain, e.g. C")
    parser.add_argument("--foldx-path", default="foldx", type=Path, help="Path to FoldX binary (default: foldx on PATH)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--regions-file", type=Path, help="YAML file defining regions to mutate {H: [[26,35], ...], L: [...]}")
    group.add_argument("--regions", type=str, help="Region string, e.g. 'H:26-35,50-66;L:24-40,56-62' ")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to store mutants and CSV")
    parser.add_argument("--workers", type=int, default=max(mp.cpu_count() - 1, 1), help="Parallel workers (default: all cores minus 1)")

    args = parser.parse_args()

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
    mutation_table = build_mutations_table(args.pdb, regions)
    print(f"Total mutations (including self): {len(mutation_table)}")

    # Parallel BuildModel
    pool_args = [ (m, repaired_pdb, args.foldx_path, mutants_dir) for m in mutation_table ]
    with mp.Pool(args.workers) as pool:
        results = pool.map(worker, pool_args)

    # Collect successes
    rows = []
    success = 0
    for (m, ok) in results:
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
    print(f"Successfully built {success}/{len(mutation_table)} mutant structures")

    pd.DataFrame(rows).to_csv(args.output_dir / "data.csv", index=False)
    print("GearBind data.csv written to", args.output_dir / "data.csv")


if __name__ == "__main__":
    main() 