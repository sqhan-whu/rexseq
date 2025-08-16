#!/bin/bash

# REXseq Installation Script
# This script installs REXseq and its dependencies

set -e  # Exit on any error

echo "=========================================="
echo "REXseq: RNA Editing Detection Pipeline"
echo "Installation Script"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda is not installed or not in PATH"
    echo "Please install Miniconda or Anaconda first:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Create conda environment
echo "Creating conda environment 'rexseq'..."
conda create -n rexseq python=3.8 -y

# Activate environment
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rexseq

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
conda install -c bioconda -y \
    snakemake \
    trim-galore \
    hisat2 \
    samtools \
    bamtools \
    bedtools \
    seqtk

# Install Python packages
echo "Installing Python packages..."
pip install -r requirements.txt

# Create necessary directories
echo "Creating project directories..."
mkdir -p data/fq
mkdir -p results
mkdir -p logs

# Make scripts executable
echo "Making scripts executable..."
chmod +x run_rexseq.py
chmod +x main_snakefile.smk

# Copy example config if it doesn't exist
if [ ! -f config.yaml ]; then
    echo "Creating config.yaml from example..."
    cp config_example.yaml config.yaml
    echo "WARNING: Please edit config.yaml with your specific paths and settings"
fi

echo "=========================================="
echo "Installation completed successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Activate the environment: conda activate rexseq"rexseq"
echo "2. Edit config.yaml with your paths and settings"
echo "3. Place your fastq files in data/fq/"
echo "4. Run the pipeline: python run_rexseq.py"
echo ""
echo "For more information, see README.md and QUICKSTART.md"
echo "=========================================="