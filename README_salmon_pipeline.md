# Salmon RNA-seq Quantification Pipeline

A comprehensive Python script for processing bulk RNA-seq data from mouse models using Salmon on M1 Mac.

## Features

✅ Automatic transcriptome reference download (GENCODE mouse)  
✅ Salmon index building  
✅ Batch processing of paired-end FASTQ files  
✅ GC and sequence bias correction  
✅ Automatic count and TPM matrix generation  
✅ Progress logging and error handling  
✅ Resume capability (skips already-processed samples)

## Prerequisites

### 1. Install Conda (if not already installed)

```bash
# Download Miniforge for M1 Mac
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

### 2. Create RNA-seq Environment

```bash
# Create environment with Salmon
conda create -n rnaseq -c bioconda salmon pandas python=3.10
conda activate rnaseq
```

### 3. Verify Installation

```bash
salmon --version
python --version
```

## Quick Start

### Basic Usage

```bash
python salmon_rnaseq_pipeline.py \
  --raw-data /path/to/your/fastq/files \
  --output /path/to/output/directory \
  --threads 8
```

### Advanced Options

```bash
python salmon_rnaseq_pipeline.py \
  --raw-data /Volumes/Data/HCC_PDX_RNAseq/raw_data \
  --output /Volumes/Data/HCC_PDX_RNAseq/salmon_results \
  --reference-dir ~/rnaseq_refs/mouse \
  --salmon-index ~/rnaseq_refs/mouse/mouse_salmon_index \
  --threads 10
```

## Input File Requirements

Your FASTQ files should follow this naming convention:
- `{sample_name}_R1_001.fastq.gz` (forward reads)
- `{sample_name}_R2_001.fastq.gz` (reverse reads)

Examples:
- `HCC001_tumor_R1_001.fastq.gz` and `HCC001_tumor_R2_001.fastq.gz`
- `Control_liver_R1_001.fastq.gz` and `Control_liver_R2_001.fastq.gz`

## Output Files

```
output_directory/
├── pipeline.log                      # Complete execution log
├── salmon_quants/                    # Per-sample quantification
│   ├── sample1/
│   │   ├── quant.sf                 # Transcript quantifications
│   │   ├── aux_info/
│   │   └── logs/
│   ├── sample2/
│   └── ...
├── salmon_counts_matrix.csv          # Count matrix (for DESeq2)
├── salmon_tpm_matrix.csv             # TPM matrix (normalized)
└── quantification_summary.txt        # Summary statistics
```

## Command-Line Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--raw-data` | Yes | - | Directory with FASTQ files |
| `--output` | Yes | - | Output directory |
| `--reference-dir` | No | `~/rnaseq_refs/mouse` | Reference files location |
| `--salmon-index` | No | `{reference-dir}/mouse_salmon_index` | Salmon index path |
| `--threads` | No | `8` | Number of CPU threads |
| `--skip-download` | No | `False` | Skip transcriptome download |

## Pipeline Steps

1. **Check Salmon Installation** - Verifies Salmon is available
2. **Download Reference** - Gets GENCODE mouse transcriptome (if needed)
3. **Build Index** - Creates Salmon index (if needed, ~10-15 min)
4. **Find FASTQ Pairs** - Identifies all paired-end samples
5. **Quantify Samples** - Runs Salmon on each sample (~5-10 min per sample)
6. **Aggregate Results** - Creates count and TPM matrices

## Example Workflow

### First Time Setup

```bash
# Activate environment
conda activate rnaseq

# Run pipeline (will download reference and build index)
python salmon_rnaseq_pipeline.py \
  --raw-data ~/Data/RNAseq/fastq \
  --output ~/Data/RNAseq/results \
  --threads 8
```

### Adding New Samples

```bash
# Just point to the same output directory
# Already-processed samples will be skipped automatically
python salmon_rnaseq_pipeline.py \
  --raw-data ~/Data/RNAseq/fastq_batch2 \
  --output ~/Data/RNAseq/results \
  --threads 8
```

## Downstream Analysis

### Option 1: Use in R with tximport (Recommended)

```r
# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures"))

# Load libraries
library(tximport)
library(DESeq2)

# Import Salmon results
samples <- list.files("results/salmon_quants", full.names = TRUE)
files <- file.path(samples, "quant.sf")
names(files) <- basename(samples)

# Create tx2gene mapping (transcript to gene)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("gencode.vM33.annotation.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# Import with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition)
```

### Option 2: Use Python Matrices Directly

```python
import pandas as pd

# Load count matrix
counts = pd.read_csv("results/salmon_counts_matrix.csv", index_col=0)

# Load TPM matrix
tpm = pd.read_csv("results/salmon_tpm_matrix.csv", index_col=0)

# Your analysis here...
```

## Troubleshooting

### "Salmon not found"
```bash
conda activate rnaseq
which salmon  # Should show path in conda environment
```

### Low mapping rates (<70%)
- Check if you're using the correct reference (mouse vs human)
- Verify FASTQ file quality with FastQC
- Check library type in `logs/salmon_quant.log`

### Out of memory errors
- Reduce `--threads` (try 4 or 6)
- Close other applications
- M1 base model: use 4 threads, Pro/Max: 8-10 threads

### "No R2 pair found"
- Check your file naming convention
- Make sure R1 and R2 files have identical prefixes
- Expected format: `{sample}_R1_*.fastq.gz` and `{sample}_R2_*.fastq.gz`

## Performance Tips

**For M1 Mac:**
- M1 Base (8GB): Use `--threads 4`
- M1 Pro/Max (16-32GB): Use `--threads 8-10`
- Each sample takes ~5-10 minutes depending on file size

**Resume Processing:**
The pipeline automatically skips samples that have already been processed. To reprocess a specific sample, delete its directory from `salmon_quants/`.

## Citation

If you use this pipeline, please cite:

**Salmon:**
Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417-419.

## Support

For issues or questions:
1. Check the `pipeline.log` file for detailed error messages
2. Review the Salmon documentation: https://salmon.readthedocs.io/
3. Check sample-specific logs in `salmon_quants/{sample}/logs/`

## License

This pipeline script is provided as-is for academic and research use.
