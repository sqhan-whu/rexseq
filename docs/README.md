# REXseq: RNA Editing Detection Pipeline

REXseq is a comprehensive Snakemake-based pipeline for detecting RNA editing sites from REXseq data. The pipeline specifically focuses on detecting A→G and T→C RNA editing events, which are the most common types of RNA editing in mammals.

## Overview

RNA editing is a post-transcriptional modification process that alters the nucleotide sequence of RNA molecules. REXseq provides a robust and reproducible workflow for identifying RNA editing sites from paired-end REXseq data.

### Key Features

- **Comprehensive Quality Control**: Includes read trimming, deduplication, and quality filtering
- **Strand-specific Analysis**: Separates positive and negative strand reads for accurate editing detection
- **Multiple Editing Types**: Detects both A→G (positive strand) and T→C (negative strand) editing events
- **Quality Filtering**: Implements stringent quality filters to reduce false positives
- **Annotation Integration**: Provides gene annotation and functional analysis of detected sites
- **Reproducible Workflow**: Built with Snakemake for reproducibility and scalability

## Pipeline Workflow

The REXseq pipeline consists of five main steps:

### Step 1: Read Mapping
- **trim_galore**: Quality trimming and adapter removal
- **fastuniq**: Read deduplication
- **seqtk**: Additional trimming for optimal mapping
- **hisat2**: REXseq read alignment to reference genome
- **samtools/bamtools**: BAM file processing and filtering

### Step 2: Gene Hitting
- **bamtools**: Separate mate1 and mate2 reads
- **bedtools**: Intersect reads with gene annotations
- **bamtools**: Merge filtered reads

### Step 3: Strand Separation
- **samtools**: Separate positive and negative strand reads
- **bamtools**: Sort and index strand-specific BAM files

### Step 4: Site Detection
- **AG.py**: Detect A→G editing on positive strand
- **TC.py**: Detect T→C editing on negative strand
- **trans_neg2pos.py**: Transform negative strand coordinates
- **find_error2.py**: Quality filtering and site classification
- **error2sites.py**: Final site selection

### Step 5: Annotation
- **pre_annover.py**: Prepare sites for annotation
- **ANNOVAR**: Gene annotation
- **sites_table2.py**: Comprehensive annotation and reporting

## Installation

### Prerequisites

REXseq requires the following software tools:

```bash
# Core bioinformatics tools
conda install -c bioconda snakemake
conda install -c bioconda trim-galore
conda install -c bioconda hisat2
conda install -c bioconda samtools
conda install -c bioconda bamtools
conda install -c bioconda bedtools
conda install -c bioconda seqtk

# Python packages
pip install pysam numpy pyyaml

# ANNOVAR (for annotation)
# Download from: http://annovar.openbioinformatics.org/
```

### Installation

1. Clone the REXseq repository:
```bash
git clone https://github.com/your-repo/rexseq.git
cd rexseq
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

## Configuration

### Configuration File

Create a `config.yaml` file with the following parameters:

```yaml
# Sample information
samples:
  - sample1
  - sample2
  - sample3

# Computational resources
threads: 8

# Input data path
path: /path/to/your/data

# Script directory
py_script: /path/to/scripts

# Reference genome files
hg38_hisat2_index: /path/to/hg38/hisat2_index
hg38_fa: /path/to/hg38.fa
hg38_size: /path/to/hg38.fa.size

# Gene annotation
BED: /path/to/gencode.v31.annotation.bed
trans_id: /path/to/anno.gene.txt

# Annotation databases
ANNOVER: /path/to/annovar/annotate_variation.pl
ANNOVER_hg38_DATABASE: /path/to/annovar/hg38db
Alu_hg38_bed: /path/to/repeat.hg38.bed
```

### Input Data Structure

Organize your input data as follows:

```
data/
├── fq/
│   ├── sample1_R1.fq.gz
│   ├── sample1_R2.fq.gz
│   ├── sample2_R1.fq.gz
│   ├── sample2_R2.fq.gz
│   └── ...
```

## Usage

### Basic Usage

Run the pipeline with default settings:

```bash
python run_rexseq.py
```

### Advanced Usage

```bash
# Run with multiple cores
python run_rexseq.py --cores 16

# Dry run to check workflow
python run_rexseq.py --dry-run

# Use custom configuration
python run_rexseq.py --config my_config.yaml

# Use custom Snakemake file
python run_rexseq.py --snakefile my_snakefile.smk
```

### Direct Snakemake Usage

You can also run the pipeline directly with Snakemake:

```bash
# Run the complete pipeline
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 8

# Run specific steps
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 8 mapping
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 8 sites
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 8 annotation
```

## Output Files

The pipeline generates the following output files:

### Main Results
- `results/final_sites/*.site.txt`: Final RNA editing sites
- `results/final_sites/*.annotation.txt`: Annotated editing sites
- `results/annover/*.variant_function`: Gene annotations

### Intermediate Files
- `results/clean_fq/`: Quality-controlled fastq files
- `results/bam/`: Aligned BAM files
- `results/bam/hit/`: Gene-overlapping reads
- `results/bam/hit/seperate/`: Strand-separated BAM files

### Output Format

The main output file (`*.annotation.txt`) contains the following columns:

1. **ENCODE_ID**: Gene identifier
2. **Chromosome**: Chromosome name
3. **Position**: Genomic position
4. **Strand**: DNA strand (+/-)
5. **Coverage**: Read coverage at the site
6. **RT**: Number of reads supporting the edit
7. **Mutation**: Mutation type (e.g., A>G, T>C)
8. **Gene**: Gene name
9. **Region**: Genomic region type
10. **Alu**: Alu repeat annotation
11. **Upstream_seq**: Upstream sequence context
12. **Downstream_seq**: Downstream sequence context

## Quality Control

The pipeline implements several quality control measures:

### Read Quality
- Minimum quality score: 20
- Minimum read length: 30 bp
- Stringency: 3 (for adapter trimming)

### Mapping Quality
- Minimum mapping quality: 5
- Primary alignments only
- Proper pairs only
- Mapped reads only

### Editing Site Quality
- Minimum coverage: 3 reads
- Minimum editing ratio: 5% for low-confidence, 10% for high-confidence
- Maximum coverage: 5000 reads (to avoid repetitive regions)
- Quality score threshold: 26

## Troubleshooting

### Common Issues

1. **Missing dependencies**: Ensure all required tools are installed and in your PATH
2. **Memory issues**: Reduce the number of threads or increase system memory
3. **Disk space**: Ensure sufficient disk space for intermediate files
4. **File permissions**: Check that you have write permissions in the output directory

### Debug Mode

Run with verbose output for debugging:

```bash
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 1 --verbose
```

### Clean Restart

To restart the pipeline from scratch:

```bash
snakemake --snakefile main_snakefile.smk --configfile config.yaml --cores 8 --clean
```

## Performance

### Resource Requirements

- **CPU**: 8-16 cores recommended
- **Memory**: 16-32 GB RAM recommended
- **Storage**: 50-100 GB free space (depending on dataset size)
- **Time**: 2-8 hours per sample (depending on data size and system)

## Citation

If you use REXseq in your research, please cite:

```
REXseq: A comprehensive pipeline for RNA editing detection from REXseq data.
[Shaoqing Han]
[Deciphering the regulatory code of RNA inosine through enzymatic precision mapping and explainable deep learning model] [2025]
```

## License

This project is licensed under the [Non-Commercial License](LICENSE).  
Commercial use requires permission.
