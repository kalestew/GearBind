# GearBind Quick Start Guide

This guide walks you through antibody saturation mutagenesis using GearBind.

## Prerequisites

1. **FoldX**: Download from https://foldxsuite.crg.eu/ and place at `./foldx/foldx_51binaryMac`
2. **Python Environment**: Set up using the provided installation script
3. **GPU**: CUDA 11.1 compatible GPU for predictions

## Step 1: Install Dependencies

```bash
# Quick installation with UV
./install_with_uv.sh
source .venv/bin/activate
```

## Step 2: Download Model Checkpoints

```bash
cd checkpoints
gdown 1nFEjbjdlRWFwYz7LUNv_D6oLnEsZ5beJ
unzip new-gearbind-model-weights.zip
mv new-gearbind-model-weights/*.pth ./
rm -rf new-gearbind-model-weights new-gearbind-model-weights.zip
cd ..
```

## Step 3: Define CDR Regions

Create `config/my_antibody_cdrs.yaml`:
```yaml
H:  # Heavy chain CDRs
  - [26, 35]    # CDR-H1
  - [50, 66]    # CDR-H2
  - [99, 108]   # CDR-H3
L:  # Light chain CDRs
  - [24, 40]    # CDR-L1
  - [56, 62]    # CDR-L2
  - [95, 103]   # CDR-L3
```

## Step 4: Generate Mutant Structures

```bash
# Using the CR3022 example structure
python script/generate_mutants_foldx.py \
    --pdb data/6xc3_wt.pdb \
    --chain-a HL \
    --chain-b C \
    --regions-file config/my_antibody_cdrs.yaml \
    --output-dir data/cr3022_saturation \
    --workers 8
```

This will:
- Use FoldX from `./foldx/foldx_51binaryMac` (no need to specify path!)
- Repair the wild-type structure
- Generate ~1,400 mutant structures (70 positions × 20 amino acids)
- Create `data/cr3022_saturation/data.csv` for GearBind

## Step 5: Run GearBind Predictions

```bash
# Run predictions with GearBind-P model
python script/predict_gearbind.py \
    --mutants-dir data/cr3022_saturation \
    --checkpoint checkpoints/GearBindP_downstream_30.pth \
    --batch-size 16
```

Results will be saved to: `data/cr3022_saturation/GearBindP_CSVAntibodyMutants.csv`

## Step 6: Ensemble Predictions (Optional)

Run multiple models for more robust predictions:
```bash
# GearBind model
python script/predict_gearbind.py \
    --mutants-dir data/cr3022_saturation \
    --checkpoint checkpoints/GearBind_downstream_10.pth

# BindDDG model  
python script/predict_gearbind.py \
    --mutants-dir data/cr3022_saturation \
    --checkpoint checkpoints/BindDDG_downstream_4.pth
```

## Understanding the Output

The output CSV contains:
- `pdb_id`: Identifier for the complex
- `mutation`: Mutation in FoldX format (e.g., SH26A)
- `chain_a`, `chain_b`: Interacting chains
- `ddG_pred`: Predicted ΔΔGbind (more negative = better binding)

## Tips

1. **CDR Regions**: Adjust the regions in the YAML file for your antibody
2. **Memory**: For large runs, reduce `--workers` if you run out of memory
3. **Batch Size**: Adjust `--batch-size` based on your GPU memory
4. **Custom Antibody**: Replace `data/6xc3_wt.pdb` with your antibody-antigen complex

## Complete Example Script

For a fully automated workflow, run:
```bash
./example_workflow.sh
```

This script demonstrates the entire pipeline from structure to predictions. 