#!/usr/bin/python3


from math import sin
from qiskit import QuantumCircuit, execute, BasicAer
from qiskit import Aer,QuantumRegister, ClassicalRegister, execute,transpile,QuantumCircuit
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeLondonV2, FakeManilaV2, FakeJakarta
import random


#---------------------------------------------------------------------------------------------
#import QNDM package

from qndm.core import *
from qndm.tool.q_gate_count import q_counter

#---------------------------------------------------------------------------------------------
#####################################
#                                   #
#        Gradient QNDM              #      
#                                   #
#####################################

#//////////#
#   main   #
#//////////#

#QNDM
def qndm_derivative(p_deco, initial_pameter,shift_position,newspop,num_qub,num_l, ent_gate,shift,G_real_qndm,gates_tot_qndm,shots,noise,val_g, simp): 

    # Setup qubit register
    q_reg_size = num_qub + 1 #numbers of qubit (sys + det)
    q_reg_size_c = 1 #numbers of classic bit
    detect_index = num_qub #index number of the detector qubit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="QNDM")

    #quantum circuit: "QNDM for gradient"
    qndm_gradient_circuit(bc, shift_position,newspop,num_qub,num_l,val_g,detect_index,shift, simp,ent_gate)
    
    #parameters initialization
    initial_values = [p_deco/2,p_deco/2]
    for i in range(len(initial_pameter)):
        initial_values.append(initial_pameter[i])

    param_dict = dict(zip(bc.parameters, initial_values))
    circ=bc.bind_parameters(param_dict)

    # measure the detector qubit
    circ.measure(detect_index, 0)

    # Quantum gates Counter
    q_counter(gates_tot_qndm,circ.decompose().decompose().decompose())
 
    #simulated NOISE
    #backend = Aer.get_backend('aer_simulator_stabilizer')
    backend = BasicAer.get_backend('qasm_simulator')
    coupling_map = None
    basis_gates = None

    if noise == True:
        backend = FakeJakarta()
        coupling_map = backend.configuration().coupling_map
        noise_model = NoiseModel.from_backend(backend)
        basis_gates = noise_model.basis_gates

    #run quantum circuit with noise simulator
    job = execute(circ, backend=backend, shots=shots)
    result = job.result()
    data = result.get_counts(circ)

    p0,p1 = 0,0
    #extract counts
    for l in data.keys():
        if l == '0': # counts for detector in zero state
           p0 += data[l]/shots # probability of |0> in the detector state

        elif l == '1': # counts for detector in one state
           p1 += data[l]/shots # probability of |1> in the detector state

    #calculate value of derivative in one direction  
    from math import sqrt,asin,sin
    G_real_qndm[shift_position] = asin(2*p1-1)/(2*p_deco)

    return None 

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Gradient calculation QNDM
def qndm_gradient(lambda1, in_par ,G_real_qndm, newspop, num_qub, num_l, ent_gate,shift,gates_tot_qndm,noise,shots,val_g, simp = False, deri = False):

    """Calculate first order derivatives with QNDM method. \n

    #Paramenters: \n
    lambda1 -- Value of coupling paramenter
    in_par -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the derivatives information \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    gates_tot_qndm -- Empty array where the function put the cost of derivatives information \n
    noise -- If noise active noise == True \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    simp -- if True the sim method is active \n
    deri = False -- If True the function calcualte only a random element of gradient \n
    """
    for shift_position in range(len(val_g)):
        qndm_derivative(lambda1, in_par ,shift_position, newspop, num_qub, num_l, ent_gate, shift, G_real_qndm,gates_tot_qndm,shots,noise,val_g, simp)


#---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################
#                                   #
#        Gradient DM Manual         #      
#                                   #
#####################################

#//////////#
#   main   #
#//////////#

#main DM Manual
def dm_derivative(initial_pameter,shift_position,num_qub,num_l, ent_gate, shift,kk,val_g,shots,gates_tot_dm2,G_real_dm2,noise,cps,k): 


    # Setup qubit register
    q_reg_size = num_qub #numbers of qubit (sys + det)
    q_reg_size_c= num_qub #numbers of classic bit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="DM")

    #quantum circuit: "DM (manual) for gradient"
    DM_manual(bc, shift_position,num_qub,num_l,val_g,shift,kk,ent_gate)
    

    #parameters initialization
    initial_values = []
    for i in range(len(initial_pameter)):
      initial_values.append(initial_pameter[i])

    param_dict = dict(zip(bc.parameters, initial_values))
    circ=bc.bind_parameters(param_dict)
        

    # measure the system qubits
    qubit_index = []
    qubit_index2 = []

    for i in range(num_qub):
      qubit_index.append(num_qub-i-1)
      qubit_index2.append(i)


    circ.measure(qubit_index, qubit_index2)

   
    # Quantum gates Counter
    q_counter(gates_tot_dm2,circ.decompose().decompose().decompose())

        #simulated NOISE
    backend = BasicAer.get_backend('qasm_simulator')
    coupling_map = None
    basis_gates = None

    if noise == True:
        backend = FakeJakarta()
        coupling_map = backend.configuration().coupling_map
        noise_model = NoiseModel.from_backend(backend)
        basis_gates = noise_model.basis_gates
 
    #run quantum circuit with noise simulator
    job = execute(circ, backend=backend, shots=shots)
    result = job.result()
    data = result.get_counts(circ)

    #extract counts


    
    minus = 0.
    plus = 0.
    for l in data.keys():
        
        one_counter = 0 
        for jk,char in enumerate(l):
            if char == '1'and kk[num_qub-jk-1]!='I':
          
                one_counter += 1
        if one_counter % 2 != 0:
       
                minus += data[l]/shots

        else:
                plus += data[l]/shots


    mean_val = -plus+minus
    G_real_dm2[shift_position] += np.real(cps[k])*mean_val/(2*sin(shift)) # probability of |0> in the system state
    
    return None 


#---------------------------------------------------------------------------------------------------------------------------------------------------
# Gradient calculation DM manual
def dm_gradient(in_par ,G_real_dm2, spop, num_qub, num_l, ent_gate,shift,gates_tot_dm2,noise,shots,val_g, deri = False):
    """Calculate first order derivatives with Direct Measurament method. \n

    #Paramenters: \n
    in_par -- Parameters of the rotational gates \n
    G_real_dm2 -- Empty array where the function put the derivatives information \n
    spop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    gates_tot_dm2 -- Empty array where the function put the cost of derivatives information \n
    noise -- If noise active noise == True \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    deri = False -- If True the function calcualte only a random element of gradient \n
    """            
    
    for shift_position in range(len(val_g)):
        for k,kk in enumerate(spop.paulis):

            for shift_sign in range(2):
                if shift_sign == 1:
                    shift *=-1 
                        
                dm_derivative(in_par,shift_position, num_qub, num_l, ent_gate,shift, kk,val_g,shots,gates_tot_dm2,G_real_dm2,noise,spop.coeffs,k)
            shift = abs(shift)

    return None 
