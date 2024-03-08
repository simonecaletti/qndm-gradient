#!/usr/bin/python3


import numpy as np
from qiskit.quantum_info import SparsePauliOp

# get_coeffs normalized the Hamiltonian coefficient 
def get_coeffs_norm(spop):
    imp_num = len(spop.coeffs)
    sum_coeff = 0 
    for i in range(1,imp_num+1):
           
            sum_coeff += abs(np.sort(spop.coeffs)[len(spop.coeffs)-i]) #to check
            
    av_bigges_coeff = (sum_coeff/imp_num) #average bigest coefficient
    spop.coeffs = spop.coeffs/av_bigges_coeff #coefficient normalisation 
    return av_bigges_coeff
