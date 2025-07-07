#!/usr/bin/env python3
"""
Simple UV-based installation script for GearBind.
Usage: uv run install_uv_simple.py
"""
import subprocess
import sys


def run_command(cmd):
    """Execute a shell command and handle errors."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        sys.exit(1)
    return result.stdout


def main():
    print("🚀 Installing GearBind with UV...")
    
    # Install base dependencies
    print("\n📦 Installing base dependencies...")
    run_command("uv pip install -r requirements.txt")
    
    # Special handling for PyTorch Geometric
    print("\n📊 Installing PyTorch Geometric components...")
    pyg_packages = [
        "torch-scatter",
        "torch-sparse", 
        "torch-cluster",
        "torch-spline-conv",
        "torch-geometric"
    ]
    
    # Get PyTorch and CUDA versions
    torch_version = "1.8.0"
    cuda_version = "cu111"
    
    # Install PyG packages with proper wheel URLs
    for package in pyg_packages[:-1]:  # All except torch-geometric
        url = f"https://data.pyg.org/whl/torch-{torch_version}+{cuda_version}.html"
        run_command(f"uv pip install {package} -f {url}")
    
    # Install torch-geometric last
    run_command("uv pip install torch-geometric")
    
    print("\n✅ Installation complete!")
    print("\nTo test the installation:")
    print("  python -c 'import torch, torch_geometric, torchdrug; print(\"Success!\")'")


if __name__ == "__main__":
    main() 