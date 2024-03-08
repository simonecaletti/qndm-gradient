#!/usr/bin/bash

echo "Installing the hsf5 library:"

#creating a new env with conda
conda create -n qndm python=3.9.13
conda activate qndm

#Networkx package
conda install -c anaconda networkx

#mat2qubit dependencies
conda install -c conda-forge json5

#Mat2qubit package
git clone https://github.com/IntelLabs/mat2qubit.git
cd

#
