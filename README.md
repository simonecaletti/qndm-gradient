# Quantum Non Demolition Measurement (QNDM) algorithm for gradient evaluation

his is the official repository for the Quantum Non-Demolition Measurement (QNDM) algorithm, as described in  [arXiv:2301.07128](https://arxiv.org/abs/2301.07128).
The QNDM algorithm enables the computation of gradients (and higher-order derivatives) of functions embedded in quantum circuits. It is particularly useful in the context of Variational Quantum Algorithms (VQAs).
QNDM is an alternative approach to the standard Direct Measurement (DM) algorithm for computing derivatives, detailed in [arXiv:2008.06517](https://arxiv.org/abs/2008.06517). The DM algorithm is the default method in most popular quantum computing libraries.

## Repository Structure

- qndm/: Contains the core functions for the QNDM algorithm.
- scripts/: Includes scripts and examples for various use cases:
- - examples/: Contains example code for using QNDM and DM derivative methods.
- - error_analysis/: Provides the code used for error analysis, as presented in the paper.
- - optimization/:  Includes example scripts for finding the minimum energy level of physical Hamiltonians using a VQA based on QNDM.

## Installation

To install the QNDM library just clone the repository using the command
```
git clone https://github.com/simonecaletti/qndm-gradient.git 
```
To use the interface with the hamiltonian library, hamlib [arXiv:2306.13126](https://arxiv.org/abs/2306.13126), supported in this repository you need to install the **mat2qubit** package [arXiv:2205.09776](https://arxiv.org/abs/2205.09776). 
The instruction are contained in the **hdf5-install.sh** script. Simply run:
```
bash hdf5-install.sh
```
To test the installation run a example script, for example 
```
python scripts/examples/example_qndm.py 
```
If the installation is working correctly a **QNDM_der_0.csv** and a **RunCard_Der.txt** file have been created. The first one contains information about the gradient computation using the QNDM algorithm, while the second is an automatically generated runcard containing the detail of the run.

## Citation

If you find out work useful, please cite our paper at:

```
@misc{minuto2024novelapproachreducederivative,
      title={A Novel Approach to Reduce Derivative Costs in Variational Quantum Algorithms}, 
      author={Giovanni Minuto and Simone Caletti and Paolo Solinas},
      year={2024},
      eprint={2404.02245},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2404.02245}, 
}

```

We also recommend referring to the theoretical development of **QNDM**, detailed in the following article:

```
@article{Solinas_2023,
   title={Quantum gradient evaluation through quantum non-demolition measurements},
   volume={77},
   ISSN={1434-6079},
   url={http://dx.doi.org/10.1140/epjd/s10053-023-00648-y},
   DOI={10.1140/epjd/s10053-023-00648-y},
   number={5},
   journal={The European Physical Journal D},
   publisher={Springer Science and Business Media LLC},
   author={Solinas, Paolo and Caletti, Simone and Minuto, Giovanni},
   year={2023},
   month=may }

```