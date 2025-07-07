# GearBind Installation Guide

This guide provides multiple installation options for GearBind dependencies.

## Prerequisites

- Python 3.8 or 3.9 (PyTorch 1.8.0 compatibility)
- CUDA 11.1 (for GPU support)
- FoldX (for mutant structure generation)

## Installation Options

### Option 1: Using UV (Recommended)

[UV](https://github.com/astral-sh/uv) is a fast Python package installer and project manager.

1. **Install UV** (if not already installed):
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. **Run the installation script**:
   ```bash
   ./install_with_uv.sh
   ```

3. **Activate the environment**:
   ```bash
   source .venv/bin/activate
   ```

### Option 2: Using UV with pyproject.toml

```bash
# Create virtual environment
uv venv --python 3.9

# Activate environment
source .venv/bin/activate

# Install dependencies
uv pip install -e .

# Run the PyTorch Geometric installation script
uv run install_uv_simple.py
```

### Option 3: Using pip with requirements.txt

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install PyTorch with CUDA 11.1
pip install torch==1.8.0+cu111 torchvision==0.9.0+cu111 -f https://download.pytorch.org/whl/torch_stable.html

# Install other requirements
pip install -r requirements.txt

# Install PyTorch Geometric (requires special handling)
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-1.8.0+cu111.html
pip install torch-geometric
```

### Option 4: Using Conda (Original method from README)

```bash
conda install pyg pytorch=1.8.0 cudatoolkit=11.1 torchdrug -c pyg -c pytorch -c conda-forge
conda install rdkit easydict pyyaml biopython gdown -c conda-forge
```

## Verifying Installation

After installation, verify everything is working:

```python
python -c "import torch, torch_geometric, torchdrug, rdkit; print('All packages imported successfully!')"
```

## Troubleshooting

1. **CUDA Issues**: Ensure your NVIDIA drivers support CUDA 11.1
2. **PyTorch Geometric**: If installation fails, try installing components one by one
3. **RDKit**: On some systems, RDKit may require conda installation instead of pip

## Next Steps

After installation, you can use the GearBind scripts:

```bash
# Generate mutant structures
python script/generate_mutants_foldx.py --help

# Run GearBind predictions
python script/predict_gearbind.py --help
``` 