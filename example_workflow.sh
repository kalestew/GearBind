#!/usr/bin/env bash
# Example workflow for antibody saturation mutagenesis with GearBind

# Step 1: Define your antibody CDR regions
# Create a YAML file defining which regions to mutate
cat > config/my_antibody_cdr_regions.yaml << EOF
H:  # Heavy chain CDRs
  - [26, 35]    # CDR-H1
  - [50, 66]    # CDR-H2  
  - [99, 108]   # CDR-H3
L:  # Light chain CDRs
  - [24, 40]    # CDR-L1
  - [56, 62]    # CDR-L2
  - [95, 103]   # CDR-L3
EOF

# Step 2: Generate all mutant structures with FoldX
# Note: FoldX binary is expected at ./foldx/foldx_51binaryMac
# Download FoldX from https://foldxsuite.crg.eu/ and place it in the foldx/ directory
echo "🧬 Generating mutant structures..."
python script/generate_mutants_foldx.py \
    --pdb data/my_antibody_complex.pdb \
    --chain-a HL \
    --chain-b C \
    --regions-file config/my_antibody_cdr_regions.yaml \
    --output-dir data/my_antibody_mutants \
    --workers 8
# Optional: Override FoldX path with --foldx-path /custom/path/to/foldx

# This creates:
# - data/my_antibody_mutants/repaired.pdb (FoldX-repaired WT)
# - data/my_antibody_mutants/mutants/*.pdb (all mutants)
# - data/my_antibody_mutants/data.csv (mutation list)

# Step 3: Download GearBind model checkpoints (if not already done)
if [ ! -f checkpoints/GearBindP_downstream_30.pth ]; then
    echo "📥 Downloading GearBind checkpoints..."
    cd checkpoints
    gdown 1nFEjbjdlRWFwYz7LUNv_D6oLnEsZ5beJ
    unzip -q new-gearbind-model-weights.zip
    mv new-gearbind-model-weights/*.pth ./
    rm -rf new-gearbind-model-weights new-gearbind-model-weights.zip
    cd ..
fi

# Step 4: Run GearBind predictions
echo "🔮 Running GearBind predictions..."
python script/predict_gearbind.py \
    --mutants-dir data/my_antibody_mutants \
    --checkpoint checkpoints/GearBindP_downstream_30.pth \
    --batch-size 16

# Results will be saved to:
# data/my_antibody_mutants/GearBindP_CSVAntibodyMutants.csv

# Step 5: Optional - Run with different models for ensemble predictions
echo "🎯 Running ensemble predictions..."
for model in GearBind_downstream_10.pth BindDDG_downstream_4.pth; do
    python script/predict_gearbind.py \
        --mutants-dir data/my_antibody_mutants \
        --checkpoint checkpoints/$model \
        --batch-size 16
done

echo "✅ Workflow complete! Check results in data/my_antibody_mutants/" 