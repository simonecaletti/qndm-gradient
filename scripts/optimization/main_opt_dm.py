#!/usr/bin/python3

#==========================================================#
#
# QISKIT code to optimize with the DM method
#==========================================================#

# 
# packages importation
#----------------------------------------------------------#

import sys, os
import chime
import numpy as np
from datetime import datetime
from tqdm.auto import tqdm
import inspect
import time
import pandas as pd
from math import pi, sqrt, sin
from multiprocessing import Pool
from dotenv import load_dotenv

# Carica le variabili d'ambiente dal file .env
load_dotenv()

# Legge la variabile e la aggiunge a sys.path
project_path = os.getenv("QNDM_PATH")
if project_path and project_path not in sys.path:
    sys.path.append(project_path)


# 
# QNDM packages importation
#----------------------------------------------------------#

from qndm.core import *


# 
# utils packages importation
#----------------------------------------------------------#

from qndm.utils import *


#==========================================================#
#
# MAIN CODE
#==========================================================#

# backend
backend = AerSimulator()

# 
# parameters read from .sh
#----------------------------------------------------------#

n_qubits   = int(sys.argv[1])    # number of qubits
n_layers   = int(sys.argv[2])    # number of layers
shots      = int(sys.argv[3])    # numver of shots
lay_u      = int(sys.argv[4])    # inside a layer: number of rotational layers
ent_gate   = int(sys.argv[5])    # entanglement layer (if 0 ---> CNOT, otherwise if 1 ---> SWAP)
shift      = float(sys.argv[6])  # shift for the parameter-shift rule 
noise      = sys.argv[7]         # if False ---> no noise, otherwise Fakesimulator      
lr         = float(sys.argv[8])  # learning rate
n_iter     = int(sys.argv[9])    # number of optimization iterations
read_pars  = sys.argv[10]        # flag that if it's True reads the parameters from the file, otherwise are randomly generated
read_gates = sys.argv[11]        # flag that if it's True reads the gates_var  from the file, otherwise are randomly generated
output_dir = sys.argv[12]        # path to the output directory
method     = 'dm'                # either 'qndm' or 'dm'
 


#==========================================================#
#
# Hamiltonian H
#==========================================================#

hamlib_ = True

if hamlib_ == True:
    # Input to use with the qndm.hamiltonians.hydrogen package
    shape = 'linear'  # linear, pyramid, ring, sheet
    r = 0.6  # between 0.5 and 2.0 (step 0.1)
    key = f"/ham_BK-{n_qubits}/"
    
    # get the model (we need for first because sets nqubit)
    model = get_model_lithium(key)

    PS = model["PS"]
    cps = model["cps"]
    pauli_string = len(cps)
    n_qubits = model["nqubit"]
  
elif hamlib_ == False:
    
    # select the pauli string number
    pauli_string = 2

    # hamiltonians H
    PS, cps = get_hamiltonian(n_qubits, pauli_string)

else:
    print('Error')



spop = get_SparsePauliOp(PS, cps)  # spop = sparse pauli operator
print(spop)


#==========================================================#
#
# Quantum circuit parameters
#==========================================================#


# total number of parameters 
n_pars = lay_u * n_layers * n_qubits 


val_g = read_gates_var(output_dir, read_pars, n_qubits, n_layers)
val_g = np.reshape(val_g, n_layers * n_qubits)
val_g = np.array(val_g)

parameters = read_parameters(output_dir, read_gates, n_qubits, n_layers)
parameters = np.reshape(parameters, n_layers * n_qubits)


#==========================================================#
#
# RUN CARD printing
#==========================================================#

print(f'printing runcard into {output_dir}')
print_run_card(output_dir, n_qubits, n_layers, parameters, val_g, shots, spop,  ent_gate, lr, method, 0, shift)




#==========================================================#
#
# Optimization part
#==========================================================#

cost_list = np.zeros(n_iter)   # List of the cost function values for each iteration
time_list = np.zeros(n_iter)   # List of the computational time for each iteration


for i in tqdm(range(n_iter), desc = 'Iteration : ', leave=True):

    start_time = time.time()
    
    cost_list[i] = cost_function(parameters, n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop, shots)
    
    print(f'{i + 1}/{n_iter}, {cost_list[i]}')

    # qndm gradient 
    grad_dm = dm_gradient(parameters, spop, n_qubits, n_layers, lay_u, ent_gate, shift, shots, val_g)
    
    # gradient descent optimization step 
    parameters = parameters - lr * grad_dm

    end_time = time.time()

    time_list[i] = float(end_time - start_time)

    if i % 10 == 0:
        print_opt_file(output_dir, i, cost_list[i], parameters, method, loop = True)



print_file(output_dir, n_iter, cost_list, time_list, method)

print_opt_file(output_dir, n_iter, cost_list[-1], parameters, method, loop = False)

print('done')


