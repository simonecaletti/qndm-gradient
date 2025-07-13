#!/usr/bin/python3

#==========================================================#
#
# QISKIT code to verify compatibility between the QNDM 
# gradient calculated using parallel and sequential methods
#==========================================================#

# 
# generic packages importation
#----------------------------------------------------------#

import sys, os
import chime
from datetime import datetime
from tqdm.auto import tqdm
import inspect
import time
import pandas as pd
from math import pi, sqrt, sin
from multiprocessing import Pool, cpu_count
from concurrent.futures import ThreadPoolExecutor

# Carica le variabili d'ambiente dal file .env
load_dotenv()

# Legge la variabile e la aggiunge a sys.path
project_path = os.getenv("QNDM_PATH")
if project_path and project_path not in sys.path:
    sys.path.append(project_path)


# 
# QISKIT packages importation
#----------------------------------------------------------#

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import Estimator
from qiskit import transpile
from qiskit_aer import Aer, AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import *
from qiskit.circuit import Parameter, ParameterVector, QuantumCircuit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from qiskit_ibm_runtime import SamplerV2


# 
# QNDM packages importation
#----------------------------------------------------------#

from qndm.core import *
from qndm.hamiltonians.examples import *
from qndm.utils.error import get_qndm_error
from qndm.layers.unitaries_gradient import *
from qndm.utils.error import get_dm_error
from qndm.layers.unitaries_gradient import *


# 
# utils packages importation
#----------------------------------------------------------#

from utils.tools import *


#==========================================================#
#
# MAIN CODE
#==========================================================#

# 
# backend 
#----------------------------------------------------------#

backend = AerSimulator()  # Imposta un seed

# 
# parameters read from .sh
#----------------------------------------------------------#

n_qubits    = int(sys.argv[1])    # number of qubits
n_layers    = int(sys.argv[2])    # number of layers
shots       = int(sys.argv[3])    # numver of shots
lay_u       = int(sys.argv[4])    # inside a layer: number of rotational layers
ent_gate    = int(sys.argv[5])    # entanglement layer (if 0 ---> CNOT, otherwise if 1 ---> SWAP)
shift       = float(sys.argv[6])  # shift for the parameter-shift rule 
noise       = sys.argv[7]         # if False ---> no noise, otherwise Fakesimulator      
lambda1     = float(sys.argv[8])  # lambda parameter for QNDM
n_gradients = int(sys.argv[9])   # number of optimization iterations
n_pauli_strings = int(sys.argv[10])   # number of Pauli strings in the hamiltonian
read_pars   = sys.argv[11]        # flag that if it's True reads the parameters from the file, otherwise are randomly generated
read_gates  = sys.argv[12]        # flag that if it's True reads the gates_var  from the file, otherwise are randomly generated
output_dir  = sys.argv[13]        # path to the output directory
method      = 'qndm'              # either 'qndm' or 'dm'
 



#==========================================================#
#
# Hamiltonian H
#==========================================================#

spop    = load_hamiltonian(output_dir, f"hamiltonian_{n_pauli_strings}_pauli_strings_dm.pkl")
newspop = load_hamiltonian(output_dir, f"hamiltonian_{n_pauli_strings}_pauli_strings_qndm.pkl")

print(spop)


#==========================================================#
#
# Quantum circuit parameters
#==========================================================#

n_pars = lay_u * n_layers * n_qubits  # total number of parameters 

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
print_run_card(output_dir, n_qubits, n_layers, parameters, val_g, shots, spop, ent_gate, 0, method, lambda1, shift)




#==========================================================#
#
# Cost function with statevector 
#==========================================================#

def expectation_value_statevector(cas, n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop):
    # Initialize the circuit
    circuit = QuantumCircuit(n_qubits)
    params = ParameterVector("theta", length=n_qubits * n_layers * lay_u)

    # Unitary transformation
    l_d = U1(val_g, params, n_qubits, n_layers, shift, 10000, ent_gate)
    qubits = list(range(n_qubits))
    circuit.compose(l_d, qubits=qubits, inplace=True)

    # Bind the provided parameters to the circuit
    param_bindings = {param: value for param, value in zip(params, cas)}
    circuit = circuit.assign_parameters(param_bindings)

    # Get the statevector of the circuit
    state = Statevector(circuit)

    # Calculate the expectation value of the observable
    expectation_value = state.expectation_value(spop).real

    return expectation_value


#==========================================================#
#
# DM derivative calculation with statevector 
#==========================================================#

def dm_derivative_statevector(args):

    i, cas, shift, n_qubits, n_layers, lay_u, ent_gate, spop, val_g = args
    
    cas_plus  = np.copy(cas)
    cas_minus = np.copy(cas)

    # Apply shift to calculate parameter shift gradient
    cas_plus[i]  += shift
    cas_minus[i] -= shift

    # Calculate expectation values for shifted parameters
    mean_plus  = expectation_value_statevector(cas_plus,  n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop)
    mean_minus = expectation_value_statevector(cas_minus, n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop)

    # Calculate gradient for the i-th parameter
    gradient_component = (mean_plus - mean_minus) / (2 * np.sin(shift))
    return gradient_component




#==========================================================#
#
# DM gradient calculation with statevector 
#==========================================================#

def dm_gradient_statevector(cas, spop, n_qubits, n_layers, lay_u, ent_gate, shift, val_g):
    """
    Calculates the gradient using the parameter shift rule.
    """

    gradient = np.zeros_like(cas)

    # Prepariamo i parametri per il multiprocessing
    args = [(i, cas, shift, n_qubits, n_layers, lay_u, ent_gate, spop, val_g) for i in range(len(cas))]
    
    # Utilizziamo Pool per parallelizzare
    with Pool(processes = 4) as pool:
        results = pool.map(dm_derivative_statevector, args)
        
    # Inseriamo i risultati nei gradienti
    for i, res in enumerate(results):
        gradient[i] = res

    return gradient





#==========================================================#
#
# Main program if we want to compare DM statevector with 
# DM shots and ALSO DM statevector with QNDM shots
#==========================================================#

# calculation and saving of the DM statevector gradient
dm_grad_statevector = dm_gradient_statevector(parameters, spop, n_qubits, n_layers, lay_u, ent_gate, shift, val_g)
dm_gradient_statevector_file = os.path.join(output_dir, "gradient_dm_statevector.txt")
np.savetxt(dm_gradient_statevector_file, dm_grad_statevector, header="DM gradient statevector", comments='')

print("dm gradient with shots saved")


# Calculation and saving of the DM/QNDM shots gradients
for i in range(n_gradients):

    dm_grad_shots = dm_gradient(parameters, spop, n_qubits, n_layers, lay_u, ent_gate, shift, shots, val_g)

    dm_gradient_shots_file = os.path.join(output_dir, f"gradient_dm_shots_{i + 1}.txt")

    np.savetxt(dm_gradient_shots_file, dm_grad_shots, header=f"DM gradient iter = {i + 1}", comments='')

    print(f"dm gradient with shots at iteration {i} saved")

    qndm_grad_shots = qndm_gradient(lambda1, parameters, newspop, n_qubits, n_layers, ent_gate, shift, shots, val_g)

    qndm_gradient_shots_file = os.path.join(output_dir, f"gradient_qndm_shots_{i + 1}.txt")

    np.savetxt(qndm_gradient_shots_file, qndm_grad_shots, header=f"QNDM gradient iter = {i + 1}", comments='')

    print(f"qndm gradient with shots at iteration {i} saved")

    print('')





