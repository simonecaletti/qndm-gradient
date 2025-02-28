
#!/usr/bin/python3

# 
# generic packages importation
#----------------------------------------------------------#

from math import sin, asin
import numpy as np
from multiprocessing import Pool, cpu_count


# 
# QISKIT packages importation
#----------------------------------------------------------#

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit.compiler import transpile
from qiskit_ibm_runtime.fake_provider import FakeManilaV2
from qiskit.circuit import Parameter, ParameterVector, QuantumCircuit


# 
# QNDM packages importation
#----------------------------------------------------------#

from qndm.derivatives.gradient.qndm import qndm_gradient_circuit
from qndm.derivatives.gradient.dm import dm_gradient_circuit
from qndm.derivatives.hessian.qndm import qndm_hessian_circuit
from qndm.derivatives.hessian.dm import dm_hessian_circuit
from qndm.hamiltonians.normalization.lam_balancing import get_lambda_balancing


# 
# QISKIT packages importation
#----------------------------------------------------------#

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp
from qiskit.primitives import Estimator
from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import *
from qiskit.circuit import Parameter, ParameterVector, QuantumCircuit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector


# 
# QNDM packages importation
#----------------------------------------------------------#

from qndm.utils import *
from qndm.hamiltonians.examples import *
from qndm.layers.unitaries_gradient import *
from qndm.hamiltonians.lithium import get_model_lithium, get_model_OH
from qndm.layers.unitaries_gradient import *

# simulator
simulator = AerSimulator()




#==========================================================#
#
# Cost function 
#==========================================================#

def cost_function(parameters, n_qubits : int, n_layers : int, lay_u : int, val_g, shift : float, ent_gate : int, spop, shots : int):
    # circuit initialization
    circuit = QuantumCircuit(n_qubits)
    params  = ParameterVector("theta", length=n_qubits * n_layers * lay_u)

    # unitary transformation
    l_d = U1(val_g, params, n_qubits, n_layers, shift, 10000, ent_gate)
    qubits = list(range(n_qubits))
    circuit.compose(l_d, qubits=qubits, inplace=True)

    # evaluation
    estimator = Estimator()  
    expectation_value = estimator.run(circuit, spop, parameters, shots=shots).result().values.real

    return expectation_value





#==========================================================#
#
# QNDM gradient calculation
#==========================================================#

def qndm_derivative(args):
    shift_position, lambda1, pars, num_qub, num_l, val_g, shift, ent_gate, newspop, shots = args

    # Setup qubit register
    q_reg_size   = num_qub + 1  # numbers of qubit (sys + det)
    q_reg_size_c = 1            # numbers of classic bit
    detect_index = num_qub      # index number of the detector qubit

    # quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="QNDM")

    # quantum circuit: "QNDM for gradient"
    qndm_gradient_circuit(bc, shift_position, newspop, num_qub, num_l, val_g, detect_index, shift, ent_gate)

    #parameters vector
    initial_values = [lambda1/2,lambda1/2]
    for i in range(len(pars)):
        initial_values.append(pars[i])

    param_dict = dict(zip(bc.parameters, initial_values))

    # Prepare the circuit with parameters
    circ = bc.assign_parameters(param_dict)

    # measure the detector qubit
    circ.measure(detect_index, 0)

    # transpile the circuit for optimization
    transpiled_circ = transpile(circ, simulator)

    # run quantum circuit with correct parameter binding
    sim_result = simulator.run(transpiled_circ, parameter_binds=[param_dict], shots=shots).result()
    data = sim_result.get_counts(transpiled_circ)

    p0, p1 = 0, 0
    # extract counts
    for l in data.keys():
        if l == '0':
            p0 += data[l] / shots  # probability of |0> in the detector state
        elif l == '1':
            p1 += data[l] / shots  # probability of |1> in the detector state

    # derivative in the direction e_(shift_position)
    gradient_component = asin(2 * p1 - 1) / (2 * lambda1)

    #print(f'gradient[{shift_position}] = {gradient_component}')
    #print(shift_position)

    return gradient_component


#---------------------------------------------------------------------------------------------------------------------------------------------------
# Gradient calculation QNDM

def qndm_gradient(lambda1, pars, newspop, num_qub, num_l, ent_gate, shift, shots, val_g):
    """Calculate first order derivatives with QNDM method. \n

    # Parameters: \n
    lambda1 -- Value of coupling parameter
    pars -- Parameters of the rotational gates \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shift of 'Parameter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    """

    gradient = np.zeros_like(pars)

    # Prepariamo i parametri per il multiprocessing
    args = [(i, lambda1, pars, num_qub, num_l, val_g, shift, ent_gate, newspop, shots) for i in range(len(pars))]

    # Utilizziamo Pool per parallelizzare
    with Pool(processes = 4) as pool:
        results = pool.map(qndm_derivative, args)

    # Inseriamo i risultati nei gradienti
    for i, res in enumerate(results):
        gradient[i] = res

    return gradient




#==========================================================#
#
# DM gradient calculation
#==========================================================#

def dm_derivative(args):
    
    i, cas, shift, n_qubits, n_layers, lay_u, ent_gate, shots, spop, val_g = args

    cas_plus  = np.copy(cas)
    cas_minus = np.copy(cas)

    # Apply shift
    cas_plus[i]  += shift
    cas_minus[i] -= shift

    # Calculate expectation values for shifted parameters
    mean_plus  = cost_function(cas_plus,  n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop, shots)
    mean_minus = cost_function(cas_minus, n_qubits, n_layers, lay_u, val_g, shift, ent_gate, spop, shots)

    # Calculate gradient for the i-th parameter
    gradient_component = (mean_plus.item() - mean_minus.item()) / (2 * np.sin(shift))
    return gradient_component



def dm_gradient(pars, spop, n_qubits, n_layers, lay_u, ent_gate, shift, shots, val_g):
    """
    Calculates the gradient using the parameter shift rule.
    """

    gradient = np.zeros_like(pars)

    # Prepariamo i parametri per il multiprocessing
    args = [(i, pars, shift, n_qubits, n_layers, lay_u, ent_gate, shots, spop, val_g) for i in range(len(pars))]

    # Utilizziamo Pool per parallelizzare
    with Pool(processes = 4) as pool:
        results = pool.map(dm_derivative, args)

    # Inseriamo i risultati nei gradienti
    for i, res in enumerate(results):
        gradient[i] = res

    return gradient



#---------------------------------------------------------------------------------------------
#####################################
#                                   #
#        Hessian QNDM               #      
#                                   #
#####################################

def qndm_derivative_hessian(lambda1, initial_parameter, sh1, sh2, newspop, num_qub, num_l, ent_gate, shift, G_real_qndm, H_real_qndm, shots, val_g, gradient_calc): 
    # Setup qubit register
    q_reg_size = num_qub + 1  # numbers of qubit (sys + det)
    q_reg_size_c = 1  # numbers of classic bit
    detect_index = num_qub  # index number of the detector qubit

    # quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")

    # balancing of lambda in function of Hamiltonian
    lambda1 = get_lambda_balancing(newspop, lambda1) 

    # hessian
    bc_hess = QuantumCircuit(q_reg, c_reg, name="QNDM_hess")

    # quantum circuit: "QNDM for hessian"
    qndm_hessian_circuit(bc_hess, sh1, sh2, newspop, num_qub, num_l, detect_index, shift, val_g, ent_gate)
    
    # parameters vector
    initial_values_hess = [lambda1 / 2] * 4 + list(initial_parameter)
    param_dict = dict(zip(bc_hess.parameters, initial_values_hess))
    circ_hess = bc_hess.bind_parameters(param_dict)

    # measure the detector qubit
    circ_hess.measure(detect_index, 0)

    # simulator
    simulator = AerSimulator()

    # run quantum circuit
    sim_result = simulator.run(circ_hess, shots=shots).result()
    data = sim_result.get_counts(circ_hess)

    p0, p1 = 0, 0
    # extract counts
    for l in data.keys():
        if l == '0':
            p0 += data[l] / shots  # probability of |0> in the detector state
        elif l == '1':
            p1 += data[l] / shots  # probability of |1> in the detector state

    # Hessian in the direction e_(sh1, sh2)
    if gradient_calc:
        G_real_qndm[sh1] += 2 * (2 * p1 - 1) / (4 * lambda1)
    H_real_qndm[sh1][sh2] += 2 * (2 * p1 - 1) / (4 * lambda1)

    return None 

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation QNDM 
def qndm_hessian(lambda1, pars, newspop, num_qub, num_l, ent_gate, shift, shots, val_g):
    """Calculate second order derivatives with QNDM method. \n

    # Parameters: \n
    lambda1 -- Value of coupling parameter \n
    pars -- Parameters of the rotational gates \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shift of 'Parameter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    """
    
    G_real_qndm = np.zeros(len(val_g))
    H_real_qndm = np.zeros((len(val_g), len(val_g)))

    for sh1 in range(len(val_g)):
        for sh2 in range(len(val_g)):
            qndm_derivative_hessian(lambda1, pars, sh1, sh2, newspop, num_qub, num_l, ent_gate, shift, G_real_qndm, H_real_qndm, shots, val_g, sh1 == sh2)
    
    return None

#---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################
#                                   #
#            Hessian DM            #      
#                                   #
#####################################

def dm_derivative_hessian(initial_parameter, sh1, sh2, num_qub, num_l, ent_gate, shift, kk, val_g, shots, H_real_dm, cps, k): 
    # Setup qubit register
    q_reg_size = num_qub  # numbers of qubit
    q_reg_size_c = num_qub  # numbers of classic bit

    # quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="DM_hess")

    # quantum circuit: "DM (manual) for hessian"
    dm_hessian_circuit(bc, sh1, sh2, num_qub, num_l, val_g, shift, kk, ent_gate)

    # parameters vector
    initial_values = initial_parameter
    param_dict = dict(zip(bc.parameters, initial_values))
    circ = bc.bind_parameters(param_dict) 

    # measuring lists
    qubit_index = list(range(num_qub - 1, -1, -1))
    qubit_index2 = list(range(num_qub))

    # measure the system qubits
    circ.measure(qubit_index, qubit_index2)

    # simulator
    simulator = AerSimulator()

    # run quantum circuit
    sim_result = simulator.run(circ, shots=shots).result()
    data = sim_result.get_counts(circ)

    # extract counts    
    minus = 0.
    plus = 0.
    
    for l in data.keys():
        one_counter = sum(1 for jk, char in enumerate(l) if char == '1' and kk[num_qub - jk - 1] != 'I')
        
        if one_counter % 2 != 0:
            minus += data[l] / shots
        else:
            plus += data[l] / shots

    # Hessian estimation cost function for a single pauli string
    mean_val = -plus + minus
    H_real_dm[sh1][sh2] += np.real(cps[k]) * mean_val / (4 * sin(shift))  # probability of |0> in the system state
    
    return None 

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation DM 
def dm_hessian(pars, H_real_dm, spop, num_qub, num_l, ent_gate, shift, shots, val_g):
    """Calculate second order derivatives with Direct Measurement method. \n

    # Parameters: \n
    pars -- Parameters of the rotational gates \n
    H_real_dm -- Empty array where the function puts the Hessians information \n
    spop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    shift -- Value of shift of 'Parameter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    """            
    
    for sh1 in range(len(val_g)):
        for k, kk in enumerate(spop.paulis):
            for sh2 in range(len(val_g)):
                for shift_sign in range(2):
                    if shift_sign == 1:
                        shift *= -1 
                        
                    dm_derivative_hessian(pars, sh1, sh2, num_qub, num_l, ent_gate, shift, kk, val_g, shots, H_real_dm, spop.coeffs, k)
                shift = abs(shift)

    return None 
