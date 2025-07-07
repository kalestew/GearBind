#!/usr/bin/env bash
# install_with_uv.sh - Install GearBind dependencies using UV
# Requires UV to be installed: https://github.com/astral-sh/uv

set -e  # Exit on error

echo "🔧 Setting up GearBind environment with UV..."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ UV is not installed. Please install it first:"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Create virtual environment with Python 3.9 (compatible with PyTorch 1.8)
echo "📦 Creating virtual environment..."
uv venv --python 3.9

# Activate the virtual environment
echo "🔄 Activating virtual environment..."
source .venv/bin/activate

# Install PyTorch with CUDA 11.1 support (matching README specs)
echo "🔥 Installing PyTorch with CUDA 11.1 support..."
uv pip install torch==1.8.0+cu111 torchvision==0.9.0+cu111 -f https://download.pytorch.org/whl/torch_stable.html

# Install PyTorch Geometric and dependencies
echo "📊 Installing PyTorch Geometric..."
TORCH_VERSION=$(python -c "import torch; print(torch.__version__.split('+')[0])")
CUDA_VERSION=$(python -c "import torch; print(torch.version.cuda.replace('.', ''))")

# Install torch-scatter, torch-sparse, torch-cluster, torch-spline-conv
uv pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-${TORCH_VERSION}+cu${CUDA_VERSION}.html
uv pip install torch-geometric

# Install TorchDrug
echo "💊 Installing TorchDrug..."
uv pip install torchdrug

# Install other dependencies
echo "📚 Installing remaining dependencies..."
uv pip install pandas numpy pyyaml easydict biopython gdown

# RDKit requires special handling
echo "🧪 Installing RDKit..."
uv pip install rdkit

echo "✅ Installation complete!"
echo ""
echo "To activate the environment in future sessions, run:"
echo "   source .venv/bin/activate"
echo ""
echo "To test the installation, run:"
echo "   python -c 'import torch, torch_geometric, torchdrug; print(\"All packages imported successfully!\")'" 