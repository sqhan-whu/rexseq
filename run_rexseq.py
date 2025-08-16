#!/usr/bin/env python3
"""
REXseq: RNA Editing Detection Pipeline Runner
A comprehensive pipeline for detecting RNA editing sites from REX-seq data.

Usage:
    python run_rexseq.py [options]

Author: Shaoqing  Han
"""

import os
import sys
import subprocess
import argparse
import yaml
from pathlib import Path

def check_dependencies():
    """Check if all required dependencies are installed."""
    required_tools = [
        'snakemake', 'trim_galore', 'hisat2', 'samtools', 
        'bamtools', 'bedtools', 'seqtk', 'python', 'perl'
    ]
    
    missing_tools = []
    for tool in required_tools:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"ERROR: Missing required tools: {', '.join(missing_tools)}")
        print("Please install the missing tools before running REXseq.")
        sys.exit(1)
    
    print("✓ All required dependencies are installed.")

def validate_config(config_file):
    """Validate the configuration file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        required_keys = [
            'samples', 'threads', 'path', 'py_script', 'BED',
            'hg38_fa', 'hg38_size', 'hg38_hisat2_index'
        ]
        
        missing_keys = []
        for key in required_keys:
            if key not in config:
                missing_keys.append(key)
        
        if missing_keys:
            print(f"ERROR: Missing required configuration keys: {', '.join(missing_keys)}")
            sys.exit(1)
        
        # Check if input directory exists
        if not os.path.exists(config['path']):
            print(f"ERROR: Input directory {config['path']} does not exist.")
            sys.exit(1)
        
        # Check if input files exist
        for sample in config['samples']:
            r1_file = os.path.join(config['path'], 'fq', f'{sample}_R1.fq.gz')
            r2_file = os.path.join(config['path'], 'fq', f'{sample}_R2.fq.gz')
            
            if not os.path.exists(r1_file):
                print(f"ERROR: Input file {r1_file} does not exist.")
                sys.exit(1)
            if not os.path.exists(r2_file):
                print(f"ERROR: Input file {r2_file} does not exist.")
                sys.exit(1)
        
        print("✓ Configuration file is valid.")
        return config
        
    except Exception as e:
        print(f"ERROR: Failed to validate configuration file: {e}")
        sys.exit(1)

def run_snakemake(snakefile, config_file, cores, dry_run=False):
    """Run the Snakemake pipeline."""
    cmd = [
        'snakemake',
        '--snakefile', snakefile,
        '--configfile', config_file,
        '--cores', str(cores),
        '--rerun-incomplete',
        '--keep-going'
    ]
    
    if dry_run:
        cmd.append('--dry-run')
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("✓ Pipeline completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Pipeline failed with exit code {e.returncode}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='REXseq: RNA Editing Detection Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_rexseq.py                    # Run with default settings
  python run_rexseq.py --cores 8          # Run with 8 cores
  python run_rexseq.py --dry-run          # Dry run to check workflow
  python run_rexseq.py --config my_config.yaml  # Use custom config
        """
    )
    
    parser.add_argument(
        '--config', 
        default='config.yaml',
        help='Configuration file (default: config.yaml)'
    )
    
    parser.add_argument(
        '--cores', 
        type=int, 
        default=1,
        help='Number of cores to use (default: 1)'
    )
    
    parser.add_argument(
        '--dry-run', 
        action='store_true',
        help='Perform a dry run without executing commands'
    )
    
    parser.add_argument(
        '--snakefile',
        default='main_snakefile.smk',
        help='Snakemake file to use (default: main_snakefile.smk)'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("REXseq: RNA Editing Detection Pipeline")
    print("=" * 60)
    
    # Check dependencies
    print("\n1. Checking dependencies...")
    check_dependencies()
    
    # Validate configuration
    print("\n2. Validating configuration...")
    config = validate_config(args.config)
    
    # Create output directories
    print("\n3. Creating output directories...")
    os.makedirs('results', exist_ok=True)
    os.makedirs('results/clean_fq', exist_ok=True)
    os.makedirs('results/bam', exist_ok=True)
    os.makedirs('results/bam/hit', exist_ok=True)
    os.makedirs('results/bam/hit/seperate', exist_ok=True)
    os.makedirs('results/final_sites', exist_ok=True)
    os.makedirs('results/annover', exist_ok=True)
    print("✓ Output directories created.")
    
    # Run pipeline
    print(f"\n4. Running REXseq pipeline...")
    if args.dry_run:
        print("DRY RUN MODE - No commands will be executed")
    
    success = run_snakemake(
        args.snakefile, 
        args.config, 
        args.cores, 
        args.dry_run
    )
    
    if success:
        print("\n" + "=" * 60)
        print("REXseq pipeline completed successfully!")
        print("=" * 60)
        print("\nOutput files:")
        print("- RNA editing sites: results/final_sites/*.site.txt")
        print("- Annotated sites: results/final_sites/*.annotation.txt")
        print("- Variant annotations: results/annover/*.variant_function")
        print("\nFor detailed results, check the 'results' directory.")
    else:
        print("\n" + "=" * 60)
        print("REXseq pipeline failed!")
        print("=" * 60)
        print("Please check the error messages above and fix any issues.")
        sys.exit(1)

if __name__ == '__main__':
    main()