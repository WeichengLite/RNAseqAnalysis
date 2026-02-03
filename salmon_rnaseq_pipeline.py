#!/usr/bin/env python3
"""
Salmon RNA-seq Quantification Pipeline
For bulk RNA-seq analysis of mouse samples on M1 Mac

Author: Oliver's RNA-seq Pipeline
Date: 2026-02-03

This script automates:
1. Checking/downloading mouse transcriptome reference
2. Building Salmon index
3. Running Salmon quantification on paired-end FASTQ files
4. Aggregating results into count matrices
"""

import os
import sys
import subprocess
import glob
import argparse
from pathlib import Path
import pandas as pd
from datetime import datetime


class SalmonPipeline:
    """Complete pipeline for Salmon RNA-seq quantification"""
    
    def __init__(self, config):
        self.config = config
        self.log_file = os.path.join(config['output_dir'], 'pipeline.log')
        
    def log(self, message):
        """Log message to both console and file"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_message = f"[{timestamp}] {message}"
        print(log_message)
        with open(self.log_file, 'a') as f:
            f.write(log_message + "\n")
    
    def check_salmon_installed(self):
        """Check if Salmon is installed"""
        try:
            result = subprocess.run(['salmon', '--version'], 
                                  capture_output=True, text=True, check=True)
            version = result.stdout.strip()
            self.log(f"✓ Salmon is installed: {version}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.log("✗ Salmon not found!")
            self.log("Please install Salmon using:")
            self.log("  conda install -c bioconda salmon")
            return False
    
    def download_reference(self):
        """Download mouse transcriptome reference if not exists"""
        ref_dir = self.config['reference_dir']
        os.makedirs(ref_dir, exist_ok=True)
        
        transcriptome_file = os.path.join(ref_dir, "gencode.vM33.transcripts.fa.gz")
        
        if os.path.exists(transcriptome_file):
            self.log(f"✓ Transcriptome reference already exists: {transcriptome_file}")
            return transcriptome_file
        
        self.log("Downloading mouse transcriptome from GENCODE...")
        url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.transcripts.fa.gz"
        
        try:
            subprocess.run(['wget', '-P', ref_dir, url], check=True)
            self.log(f"✓ Downloaded transcriptome to {transcriptome_file}")
            return transcriptome_file
        except subprocess.CalledProcessError as e:
            self.log(f"✗ Failed to download transcriptome: {e}")
            self.log("You can manually download from:")
            self.log(f"  {url}")
            sys.exit(1)
    
    def build_salmon_index(self, transcriptome_file):
        """Build Salmon index from transcriptome"""
        index_dir = self.config['salmon_index']
        
        if os.path.exists(index_dir) and os.path.isdir(index_dir):
            # Check if index is complete
            required_files = ['hash.bin', 'sa.bin', 'txpInfo.bin']
            if all(os.path.exists(os.path.join(index_dir, f)) for f in required_files):
                self.log(f"✓ Salmon index already exists: {index_dir}")
                return True
        
        self.log("Building Salmon index (this may take 10-15 minutes)...")
        
        cmd = [
            'salmon', 'index',
            '-t', transcriptome_file,
            '-i', index_dir,
            '-k', '31',
            '--gencode'
        ]
        
        self.log(f"Command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.log(f"✓ Salmon index built successfully: {index_dir}")
            return True
        except subprocess.CalledProcessError as e:
            self.log(f"✗ Failed to build Salmon index: {e}")
            return False
    
    def find_fastq_pairs(self):
        """
        Find paired-end FASTQ files.
        Assumes naming: {sample}_R1_*.fastq.gz and {sample}_R2_*.fastq.gz
        """
        data_dir = self.config['raw_data_dir']
        self.log(f"Searching for FASTQ files in: {data_dir}")
        
        # Get all R1 files
        r1_pattern = os.path.join(data_dir, "*_R1_*.fastq.gz")
        r1_files = sorted(glob.glob(r1_pattern))
        
        if not r1_files:
            self.log(f"✗ No R1 FASTQ files found matching pattern: {r1_pattern}")
            return []
        
        pairs = []
        for r1_file in r1_files:
            # Derive R2 filename (handle various naming conventions)
            r2_file = r1_file.replace("_R1_", "_R2_")
            
            if os.path.exists(r2_file):
                # Extract sample name
                sample_name = os.path.basename(r1_file).split("_R1_")[0]
                pairs.append((sample_name, r1_file, r2_file))
            else:
                self.log(f"⚠ Warning: No R2 pair found for {os.path.basename(r1_file)}")
        
        return pairs
    
    def run_salmon_quant(self, sample_name, r1_file, r2_file):
        """Run Salmon quantification for a single sample"""
        output_dir = self.config['output_dir']
        salmon_index = self.config['salmon_index']
        threads = self.config['threads']
        
        sample_output = os.path.join(output_dir, 'salmon_quants', sample_name)
        
        # Skip if already processed
        quant_file = os.path.join(sample_output, 'quant.sf')
        if os.path.exists(quant_file):
            self.log(f"⊙ {sample_name} already processed, skipping...")
            return True
        
        cmd = [
            'salmon', 'quant',
            '-i', salmon_index,
            '-l', 'A',  # Auto-detect library type
            '-1', r1_file,
            '-2', r2_file,
            '-o', sample_output,
            '--validateMappings',
            '--gcBias',  # Correct for GC bias
            '--seqBias',  # Correct for sequence-specific bias
            '-p', str(threads)
        ]
        
        self.log(f"Processing {sample_name}...")
        self.log(f"  R1: {os.path.basename(r1_file)}")
        self.log(f"  R2: {os.path.basename(r2_file)}")
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            self.log(f"✓ {sample_name} completed successfully")
            
            # Extract and log mapping rate
            log_file = os.path.join(sample_output, 'logs', 'salmon_quant.log')
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'Mapping rate' in line:
                            self.log(f"  {line.strip()}")
                            break
            return True
        except subprocess.CalledProcessError as e:
            self.log(f"✗ Error processing {sample_name}")
            self.log(f"  Error: {e.stderr.decode() if e.stderr else 'Unknown error'}")
            return False
    
    def aggregate_counts(self):
        """Aggregate Salmon quant.sf files into count matrices"""
        salmon_dir = os.path.join(self.config['output_dir'], 'salmon_quants')
        output_dir = self.config['output_dir']
        
        self.log("Aggregating count data...")
        
        # Find all quant.sf files
        quant_files = sorted(glob.glob(os.path.join(salmon_dir, "*/quant.sf")))
        
        if not quant_files:
            self.log("✗ No quant.sf files found for aggregation!")
            return False
        
        self.log(f"Found {len(quant_files)} quantification files")
        
        # Read first file to get transcript IDs
        first_df = pd.read_csv(quant_files[0], sep="\t")
        
        # Create count matrix (NumReads - estimated counts)
        count_matrix = pd.DataFrame({"transcript_id": first_df["Name"]})
        
        # Create TPM matrix
        tpm_matrix = pd.DataFrame({"transcript_id": first_df["Name"]})
        
        # Collect data from all samples
        for quant_file in quant_files:
            sample_name = os.path.basename(os.path.dirname(quant_file))
            df = pd.read_csv(quant_file, sep="\t")
            
            # NumReads for count-based analysis (DESeq2)
            count_matrix[sample_name] = df["NumReads"].round().astype(int)
            
            # TPM for normalized expression
            tpm_matrix[sample_name] = df["TPM"]
        
        # Save matrices
        count_file = os.path.join(output_dir, "salmon_counts_matrix.csv")
        tpm_file = os.path.join(output_dir, "salmon_tpm_matrix.csv")
        
        count_matrix.to_csv(count_file, index=False)
        tpm_matrix.to_csv(tpm_file, index=False)
        
        self.log(f"✓ Count matrix saved: {count_file}")
        self.log(f"  Shape: {count_matrix.shape}")
        self.log(f"  Samples: {len(count_matrix.columns) - 1}")
        self.log(f"✓ TPM matrix saved: {tpm_file}")
        
        # Create a summary statistics file
        summary_file = os.path.join(output_dir, "quantification_summary.txt")
        with open(summary_file, 'w') as f:
            f.write("Salmon Quantification Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total samples: {len(count_matrix.columns) - 1}\n")
            f.write(f"Total transcripts: {len(count_matrix)}\n\n")
            
            # Per-sample total counts
            f.write("Per-sample read counts:\n")
            for col in count_matrix.columns[1:]:
                total_counts = count_matrix[col].sum()
                f.write(f"  {col}: {total_counts:,.0f} reads\n")
        
        self.log(f"✓ Summary statistics saved: {summary_file}")
        
        return True
    
    def run_pipeline(self):
        """Execute the complete pipeline"""
        self.log("=" * 70)
        self.log("Starting Salmon RNA-seq Quantification Pipeline")
        self.log("=" * 70)
        
        # Step 1: Check Salmon installation
        self.log("\n[STEP 1] Checking Salmon installation...")
        if not self.check_salmon_installed():
            return False
        
        # Step 2: Download reference
        self.log("\n[STEP 2] Preparing transcriptome reference...")
        transcriptome_file = self.download_reference()
        
        # Step 3: Build Salmon index
        self.log("\n[STEP 3] Building Salmon index...")
        if not self.build_salmon_index(transcriptome_file):
            return False
        
        # Step 4: Find FASTQ pairs
        self.log("\n[STEP 4] Finding FASTQ files...")
        pairs = self.find_fastq_pairs()
        
        if not pairs:
            self.log("✗ No valid FASTQ pairs found!")
            return False
        
        self.log(f"✓ Found {len(pairs)} sample pairs:")
        for sample_name, _, _ in pairs:
            self.log(f"  - {sample_name}")
        
        # Step 5: Run Salmon quantification
        self.log("\n[STEP 5] Running Salmon quantification...")
        success_count = 0
        for sample_name, r1_file, r2_file in pairs:
            if self.run_salmon_quant(sample_name, r1_file, r2_file):
                success_count += 1
        
        self.log(f"\n✓ Successfully processed {success_count}/{len(pairs)} samples")
        
        # Step 6: Aggregate results
        self.log("\n[STEP 6] Aggregating results...")
        if not self.aggregate_counts():
            return False
        
        self.log("\n" + "=" * 70)
        self.log("Pipeline completed successfully!")
        self.log("=" * 70)
        self.log(f"\nResults are in: {self.config['output_dir']}")
        self.log(f"Log file: {self.log_file}")
        
        return True


def main():
    parser = argparse.ArgumentParser(
        description='Salmon RNA-seq Quantification Pipeline for Mouse Samples',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python salmon_rnaseq_pipeline.py \\
    --raw-data /path/to/fastq/files \\
    --output /path/to/output \\
    --threads 8

The script will:
  1. Download mouse transcriptome (if needed)
  2. Build Salmon index (if needed)
  3. Quantify all paired-end FASTQ files
  4. Generate count and TPM matrices
        """
    )
    
    parser.add_argument(
        '--raw-data',
        required=True,
        help='Directory containing raw FASTQ files (*_R1_*.fastq.gz, *_R2_*.fastq.gz)'
    )
    
    parser.add_argument(
        '--output',
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '--reference-dir',
        default=os.path.expanduser('~/rnaseq_refs/mouse'),
        help='Directory for reference files (default: ~/rnaseq_refs/mouse)'
    )
    
    parser.add_argument(
        '--salmon-index',
        default=None,
        help='Path to Salmon index (default: reference-dir/mouse_salmon_index)'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads for Salmon (default: 8)'
    )
    
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip transcriptome download (use existing reference)'
    )
    
    args = parser.parse_args()
    
    # Set up configuration
    config = {
        'raw_data_dir': os.path.abspath(args.raw_data),
        'output_dir': os.path.abspath(args.output),
        'reference_dir': os.path.abspath(args.reference_dir),
        'salmon_index': args.salmon_index or os.path.join(
            os.path.abspath(args.reference_dir), 
            'mouse_salmon_index'
        ),
        'threads': args.threads,
        'skip_download': args.skip_download
    }
    
    # Create output directory
    os.makedirs(config['output_dir'], exist_ok=True)
    
    # Run pipeline
    pipeline = SalmonPipeline(config)
    success = pipeline.run_pipeline()
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
