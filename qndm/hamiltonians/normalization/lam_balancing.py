#!/usr/bin/python3


import numpy as np
from qiskit.quantum_info import SparsePauliOp

# get_lambda_balancing tunes the lambda term in function of hamiltonian
def get_lambda_balancing(spop,lambda1):

        #sort the coefficients
        sorted_indices = np.argsort(np.abs(spop.coeffs))
        sorted_coeffs = spop.coeffs[sorted_indices]

        n_coeff = len(spop.coeffs[-2:])

        #norm of coefficients
        norm_coeff = np.real(sum(sorted_coeffs[-2:]))/n_coeff

        #new lambda
        lambda1 = lambda1/norm_coeff

        return lambda1
