#!/usr/bin/bash

# Exit immediately if a command exits with a non-zero status,
# and treat unset variables as errors
set -euo pipefail

echo "Installing the hsf5 library:"

# Activate conda (assuming conda is installed and configured)
# This assumes conda is in your PATH
source "$(conda info --base)/etc/profile.d/conda.sh"

#creating a new env with conda
conda create -n qndm2 python=3.9.13
conda activate qndm2

#Networkx package
conda install -c anaconda networkx

#mat2qubit dependencies
conda install -c conda-forge json5

#Mat2qubit package
git clone https://github.com/IntelLabs/mat2qubit.git

pip install -r requirements.txt

cd
