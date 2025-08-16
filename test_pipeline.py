#!/usr/bin/env python3
"""
REXseq Pipeline Test Script
Tests the pipeline configuration and dependencies
"""

import os
import sys
import subprocess
import yaml
from pathlib import Path

def test_dependencies():
    """Test if all required dependencies are available."""
    print("Testing dependencies...")
    
    required_tools = [
        'snakemake', 'trim_galore', 'hisat2', 'samtools', 
        'bamtools', 'bedtools', 'seqtk', 'python', 'perl'
    ]
    
    missing_tools = []
    for tool in required_tools:
        try:
            if tool == 'seqtk':
                subprocess.run([tool], capture_output=True)
            else:
                subprocess.run([tool, '--version'], capture_output=True)
            print(f"✓ {tool}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"✗ {tool} - NOT FOUND")
            missing_tools.append(tool)
    
    if missing_tools:
        print(f"\nERROR: Missing tools: {', '.join(missing_tools)}")
        return False
    
    print("✓ All dependencies are available")
    return True

def test_config():
    """Test configuration file."""
    print("\nTesting configuration...")
    
    if not os.path.exists('config.yaml'):
        print("✗ config.yaml not found")
        print("Please copy config_example.yaml to config.yaml and edit it")
        return False
    
    try:
        with open('config.yaml', 'r') as f:
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
            print(f"✗ Missing configuration keys: {', '.join(missing_keys)}")
            return False
        
        print("✓ Configuration file is valid")
        return True
        
    except Exception as e:
        print(f"✗ Error reading config.yaml: {e}")
        return False

def test_input_data():
    """Test if input data exists."""
    print("\nTesting input data...")
    
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        data_path = config['path']
        if not os.path.exists(data_path):
            print(f"✗ Data directory not found: {data_path}")
            return False
        
        fq_dir = os.path.join(data_path, 'fq')
        if not os.path.exists(fq_dir):
            print(f"✗ Fastq directory not found: {fq_dir}")
            return False
        
        # Check for sample files
        samples = config['samples']
        missing_files = []
        for sample in samples:
            r1_file = os.path.join(fq_dir, f'{sample}_R1.fq.gz')
            r2_file = os.path.join(fq_dir, f'{sample}_R2.fq.gz')
            
            if not os.path.exists(r1_file):
                missing_files.append(r1_file)
            if not os.path.exists(r2_file):
                missing_files.append(r2_file)
        
        if missing_files:
            print(f"✗ Missing input files: {', '.join(missing_files)}")
            return False
        
        print(f"✓ Input data found for {len(samples)} samples")
        return True
        
    except Exception as e:
        print(f"✗ Error checking input data: {e}")
        return False

def test_reference_files():
    """Test if reference files exist."""
    print("\nTesting reference files...")
    
    try:
        with open('config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        
        reference_files = [
            config['hg38_hisat2_index'],
            config['hg38_fa'],
            config['hg38_size'],
            config['BED'],
            config['trans_id']
        ]
        
        missing_files = []
        for file_path in reference_files:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
        
        if missing_files:
            print(f"✗ Missing reference files: {', '.join(missing_files)}")
            return False
        
        print("✓ All reference files found")
        return True
        
    except Exception as e:
        print(f"✗ Error checking reference files: {e}")
        return False

def test_snakemake():
    """Test Snakemake workflow."""
    print("\nTesting Snakemake workflow...")
    
    try:
        # Run a dry run
        cmd = [
            'snakemake',
            '--snakefile', 'main_snakefile.smk',
            '--configfile', 'config.yaml',
            '--cores', '1',
            '--dry-run'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print("✓ Snakemake workflow is valid")
            return True
        else:
            print(f"✗ Snakemake workflow error:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
            return False
            
    except Exception as e:
        print(f"✗ Error testing Snakemake: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 50)
    print("REXseq Pipeline Test Suite")
    print("=" * 50)
    
    tests = [
        test_dependencies,
        test_config,
        test_input_data,
        test_reference_files,
        test_snakemake
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test failed with exception: {e}")
    
    print("\n" + "=" * 50)
    print(f"Test Results: {passed}/{total} tests passed")
    print("=" * 50)
    
    if passed == total:
        print("✓ All tests passed! The pipeline is ready to run.")
        print("\nTo run the pipeline:")
        print("  python run_rexseq.py")
    else:
        print("✗ Some tests failed. Please fix the issues before running the pipeline.")
        sys.exit(1)

if __name__ == '__main__':
    main()