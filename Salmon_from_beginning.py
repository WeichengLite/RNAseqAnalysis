# If you don't have conda, install miniforge (ARM-native for M1)
# Download from: https://github.com/conda-forge/miniforge

# Create a new environment for RNA-seq analysis
conda create -n rnaseq salmon python=3.10
conda activate rnaseq

# Verify installation:
salmon --version

