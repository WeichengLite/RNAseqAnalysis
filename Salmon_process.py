# For MacBook users, highly recommend having its own environment with conda for each different work
# Ask Gemini/Chat/Claude, to build a salmon work environment


#To check if you're in the correct conda environment:
conda env list
#The active environment will have an asterisk (*) next to it. You want to be in the salmon environment.
#If you're not in it, activate it first:
conda activate salmon


#check: this will show you the path to salmon if it's installed, or return nothing if it's not found.
which salmon
#or check: If salmon is installed, this will show the version number. If not, you'll get a "command not found" error.
salmon --version


#To integrate Salmon into Jupeter Lab notebook
# 1. Make sure you are in the environment
conda activate salmon

# 2. Install the tool that connects Conda to Jupyter
conda install python ipykernel

# 3. Create the "bridge" (the kernel)
python -m ipykernel install --user --name salmon --display-name "Bio (salmon)"

