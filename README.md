# Quantum Non Demolition Measurement (QNDM) algorithm for gradient evaluation

Official repository for the QNDM algorithm described in [arXiv:2301.07128](https://arxiv.org/abs/2301.07128).
The QNDM algorithm allows to compute the gradient (and higher derivatives) of a function embedded in a quantum circuit. It is particularly useful in the context of Variational Quantum Algorithm (VQA).
The standard method to compute derivative is the so-called Direct Measurement algorithm [arXiv:2008.06517](https://arxiv.org/abs/2008.06517) which is the default options in the most popular QC libraries.

## Installation

To install the QNDM library just clone the repository using the command
```
git clone https://github.com/simonecaletti/qndm-gradient.git 
```
All the function for the QNDM algorithm are defined in the **qndm/** folder. The **test_scripts/** folder instead contains a set of scripts with some simple tests, like QNDM and DM derivatives evaluation and optimization tasks.

To use the interface with the hamlib library [arXiv:2306.13126](https://arxiv.org/abs/2306.13126) you need to install the **mat2qubit** package [arXiv:2205.09776](https://arxiv.org/abs/2205.09776). The instruction are contained in the **hdf5-install.sh** script, so just run 
```
chmod +x hdf5-install.sh     
bash hdf5-install.sh
```

To test the installation run a test script in the folder **test_scripts**, for example 
```
python3 test_qndm.py 
```
If the installation is working correctly a **QNDM_der.csv** and a **RunCard_Der.txt** file have been created in folder **output_test**. The first one contains information about the gradient computation using the QNDM algorithm, while the second is an automatically generated runcard containing the detail of the run.

