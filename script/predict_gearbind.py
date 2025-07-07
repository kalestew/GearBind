#!/usr/bin/env python3
"""
script/predict_gearbind.py
--------------------------
Run GearBind predictions on mutant structures prepared by
`generate_mutants_foldx.py`.

This script dynamically builds a YAML configuration that points GearBind
(`script/predict.py`) to the dataset produced in *output_dir* and runs the
prediction. Results are written back to *output_dir*.

Example
=======

    python script/predict_gearbind.py \
        --mutants-dir data/cr3022_mutants \
        --checkpoint checkpoints/GearBindP_downstream_30.pth

The script will create `<mutants-dir>/gearbind_config.yaml`, launch
GearBind, and save predictions to `<mutants-dir>/GearBindP_<dataset>.csv`.
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import yaml

CONFIG_NAME = "gearbind_predict.yaml"


def write_config(mutants_dir: Path, checkpoint: Path, batch_size: int = 16):
    """Create a GearBind YAML config inside *mutants_dir*."""
    config = {
        "output_dir": str(mutants_dir),
        "dataset": {
            "class": "CSVAntibodyMutants",
            "path": str(mutants_dir / "data.csv"),
            "transform": {
                "class": "ProteinProteinInterfaceTransform",
                "max_radius": 15,
                "num_nearest": 128,
            },
        },
        "task": {
            "class": "PropertyPrediction",
            "model": {
                "class": "GearBindModel",
                "path": str(checkpoint),
            },
            "num_mlp_layer": 2,
            "graph_batch_size": 1,
        },
        "engine": {
            "gpus": [0],
            "batch_size": batch_size,
            "log_interval": 50,
        },
        "metric": ["mae", "rmse", "pearsonr", "spearmanr"],
    }
    cfg_path = mutants_dir / CONFIG_NAME
    with cfg_path.open("w") as fh:
        yaml.safe_dump(config, fh)
    return cfg_path


def ensure_dataset_class():
    """Patch gearbind.dataset with a lightweight CSV dataset class if missing."""
    import importlib
    import types

    ds_mod = importlib.import_module("gearbind.dataset")
    if hasattr(ds_mod, "CSVAntibodyMutants"):
        return  # Already present

    from torchdrug import data  # type: ignore

    class CSVAntibodyMutants(data.Dataset):  # pylint: disable=too-few-public-methods
        """Minimal dataset reading protein structures from a GearBind csv."""

        def __init__(self, path: str, transform=None, **kwargs):  # noqa: D401
            super().__init__()
            import pandas as pd  # Local import to avoid mandatory dependency for CLI
            self.csv = pd.read_csv(path)
            self.transform = transform
            self.proteins = []
            for _, row in self.csv.iterrows():
                wt = data.Protein.from_pdb(Path(path).with_name(row["wt_protein"]))
                mt = data.Protein.from_pdb(Path(path).with_name(row["mt_protein"]))
                self.proteins.append({"wt": wt, "mt": mt})

        def __len__(self):
            return len(self.proteins)

        def __getitem__(self, idx):
            return self.proteins[idx]

    # Inject into module namespace so GearBind can look it up via reflection
    setattr(ds_mod, "CSVAntibodyMutants", CSVAntibodyMutants)



def main():
    parser = argparse.ArgumentParser(description="Run GearBind predictions on mutant structures")
    parser.add_argument("--mutants-dir", required=True, type=Path, help="Directory output by generate_mutants_foldx.py")
    parser.add_argument("--checkpoint", required=True, type=Path, help="Path to GearBind checkpoint (.pth)")
    parser.add_argument("--batch-size", type=int, default=16, help="Graph batch size for GearBind engine")
    args = parser.parse_args()

    ensure_dataset_class()

    cfg_path = write_config(args.mutants_dir, args.checkpoint, batch_size=args.batch_size)

    cmd = ["python", "script/predict.py", "-c", str(cfg_path)]
    print("Running GearBind:", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main() 