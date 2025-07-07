#!/usr/bin/env bash
# install_with_uv_macos.sh - Install GearBind dependencies using UV for macOS
# Adapted for macOS (including Apple Silicon)

set -e  # Exit on error

echo "🔧 Setting up GearBind environment with UV (macOS)..."

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "❌ UV is not installed. Please install it first:"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# Create virtual environment with Python 3.9
echo "📦 Creating virtual environment..."
uv venv --python 3.9

# Activate the virtual environment
echo "🔄 Activating virtual environment..."
source .venv/bin/activate

# Install PyTorch for macOS (CPU version)
echo "🔥 Installing PyTorch for macOS..."
uv pip install torch==1.8.0 torchvision==0.9.0

# Install PyTorch Geometric and dependencies
echo "📊 Installing PyTorch Geometric..."
TORCH_VERSION="1.8.0"

# Install torch-scatter, torch-sparse, torch-cluster, torch-spline-conv
# For macOS, we need to use CPU versions
uv pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-${TORCH_VERSION}+cpu.html
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
echo "Note: This is a CPU-only installation suitable for macOS."
echo "For GPU support, use a Linux machine with CUDA."
echo ""
echo "To activate the environment in future sessions, run:"
echo "   source .venv/bin/activate"
echo ""
echo "To test the installation, run:"
echo "   python -c 'import torch, torch_geometric, torchdrug; print(\"All packages imported successfully!\")'" 