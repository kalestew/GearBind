#!/usr/bin/env bash
# install_macos_conda.sh - Install GearBind on macOS using conda
# Works with Apple Silicon (ARM64) and Intel Macs

set -e

echo "🔧 Setting up GearBind environment with conda (macOS)..."

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda is not installed. Please install Miniconda or Anaconda first."
    echo "   Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Create conda environment
echo "📦 Creating conda environment 'gearbind'..."
conda create -n gearbind python=3.9 -y

# Activate environment
echo "🔄 Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate gearbind

# Install PyTorch and related packages
echo "🔥 Installing PyTorch for macOS..."
conda install pytorch==1.13.1 torchvision torchaudio -c pytorch -y

# Install PyTorch Geometric
echo "📊 Installing PyTorch Geometric..."
conda install pyg -c pyg -y

# Install TorchDrug and other dependencies
echo "💊 Installing remaining packages..."
conda install -c conda-forge torchdrug rdkit easydict pyyaml biopython gdown pandas numpy -y

echo "✅ Installation complete!"
echo ""
echo "Note: Using PyTorch 1.13.1 for macOS compatibility."
echo "The models were trained on PyTorch 1.8.0, but should work with minor warnings."
echo ""
echo "To activate the environment in future sessions, run:"
echo "   conda activate gearbind"
echo ""
echo "To test the installation, run:"
echo "   python -c 'import torch, torch_geometric, torchdrug; print(\"All packages imported successfully!\")'" 