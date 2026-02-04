# If you don't have conda, install miniforge (ARM-native for M1)
# Download from: https://github.com/conda-forge/miniforge

# Create a new environment for RNA-seq analysis
conda create -n rnaseq salmon python=3.10
conda activate rnaseq

# Verify installation:
salmon --version

### Step 2: Download Mouse Transcriptome Reference###

# Create a reference directory
mkdir -p ~/rnaseq_refs/mouse
cd ~/rnaseq_refs/mouse

# Download mouse transcriptome (GENCODE - recommended)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz

# Or from Ensembl (alternative)
# wget http://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

### Step 3: Build Salmon Index ###
# Build the index (this takes ~10-15 minutes on M1)
salmon index \
  -t gencode.vM33.transcripts.fa.gz \
  -i mouse_salmon_index \
  -k 31 \
  --gencode

###Step 4: Python Script to Process All Samples###

import os
import subprocess
from pathlib import Path
import glob

# Configuration
RAW_DATA_DIR = "/Users/weichengli/wl_RNA2_CCA_mouse/fastq1"  # Change this!
OUTPUT_DIR = "/Users/weichengli/wl_RNA2_CCA_mouse/fastq1/salmon_quants"  # Change this!
SALMON_INDEX = "/Users/weichengli/wl_RNA2_CCA_mouse/fastq1/rnaseq_refs/mouse/mouse_salmon_index"  # Change if different
THREADS = 8  # Adjust based on your M1 (Pro/Max have more cores)

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def find_fastq_pairs(data_dir):
    """
    Find paired-end FASTQ files.
    Assumes naming: {sample}_R1_001.fastq.gz and {sample}_R2_001.fastq.gz
    """
    # Get all R1 files
    r1_files = sorted(glob.glob(os.path.join(data_dir, "*_R1_*.fastq.gz")))
    
    pairs = []
    for r1_file in r1_files:
        # Derive R2 filename
        r2_file = r1_file.replace("_R1_", "_R2_")
        
        if os.path.exists(r2_file):
            # Extract sample name
            sample_name = os.path.basename(r1_file).split("_R1_")[0]
            pairs.append((sample_name, r1_file, r2_file))
        else:
            print(f"Warning: No R2 pair found for {r1_file}")
    
    return pairs

def run_salmon(sample_name, r1_file, r2_file, output_dir, salmon_index, threads):
    """
    Run Salmon quantification for a single sample
    """
    sample_output = os.path.join(output_dir, sample_name)
    
    cmd = [
        "salmon", "quant",
        "-i", salmon_index,
        "-l", "A",  # Auto-detect library type
        "-1", r1_file,
        "-2", r2_file,
        "-o", sample_output,
        "--validateMappings",
        "--gcBias",  # Correct for GC bias
        "--seqBias",  # Correct for sequence-specific bias
        "-p", str(threads)
    ]
    
    print(f"Processing {sample_name}...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
        print(f"✓ {sample_name} completed successfully\n")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error processing {sample_name}: {e}\n")

def main():
    # Find all FASTQ pairs
    print("Searching for FASTQ files...")
    pairs = find_fastq_pairs(RAW_DATA_DIR)
    
    print(f"Found {len(pairs)} sample pairs:\n")
    for sample_name, r1, r2 in pairs:
        print(f"  - {sample_name}")
    print()
    
    # Process each sample
    for sample_name, r1_file, r2_file in pairs:
        run_salmon(sample_name, r1_file, r2_file, OUTPUT_DIR, SALMON_INDEX, THREADS)
    
    print("All samples processed!")
    print(f"Results are in: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()

###Step 5: Aggregate Results with tximport (for DESeq2)###
# After Salmon completes, you'll have a quant.sf file in each sample directory. To use with DESeq2, create another Python script:

import pandas as pd
import os
import glob

SALMON_OUTPUT_DIR = "/Users/weichengli/wl_RNA2_CCA_mouse/fastq1/salmon_quants"
OUTPUT_FILE = "salmon_counts_matrix.csv"

def aggregate_salmon_counts(salmon_dir, output_file):
    """
    Aggregate Salmon quant.sf files into a count matrix
    """
    # Find all quant.sf files
    quant_files = glob.glob(os.path.join(salmon_dir, "*/quant.sf"))
    
    if not quant_files:
        print("No quant.sf files found!")
        return
    
    # Read first file to get transcript IDs
    first_df = pd.read_csv(quant_files[0], sep="\t")
    count_matrix = pd.DataFrame({"transcript_id": first_df["Name"]})
    
    # Collect counts from all samples
    for quant_file in sorted(quant_files):
        sample_name = os.path.basename(os.path.dirname(quant_file))
        df = pd.read_csv(quant_file, sep="\t")
        
        # Use NumReads (estimated counts) for count-based analysis
        count_matrix[sample_name] = df["NumReads"].values
    
    # Save
    count_matrix.to_csv(output_file, index=False)
    print(f"Count matrix saved to {output_file}")
    print(f"Shape: {count_matrix.shape}")
    print(f"Samples: {len(count_matrix.columns) - 1}")

if __name__ == "__main__":
    aggregate_salmon_counts(SALMON_OUTPUT_DIR, OUTPUT_FILE)
