#!/usr/bin/python3


from math import sin
from qiskit import QuantumCircuit, execute, BasicAer
from qiskit import Aer,QuantumRegister, ClassicalRegister, execute,transpile,QuantumCircuit
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeLondonV2, FakeManilaV2, FakeJakarta
from qndm.derivatives.gradient.qndm import qndm_gradient_circuit
from qndm.derivatives.gradient.dm import dm_gradient_circuit
from qndm.derivatives.hessian.qndm import qndm_hessian_circuit
from qndm.derivatives.hessian.dm import dm_hessian_circuit
from qndm.hamiltonians.normalization.coeff_norm import get_coeffs_norm

import random
import numpy as np
from math import asin


#---------------------------------------------------------------------------------------------
#import QNDM package
#from qndm.tool.q_gate_count import q_counter

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

    #hamiltonian normalization
    norm_coeff = get_coeffs_norm(newspop) 

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
    #q_counter(gates_tot_qndm,circ.decompose().decompose().decompose())
 
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
    G_real_qndm[shift_position] = norm_coeff*asin(2*p1-1)/(2*p_deco)

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
#            Gradient DM            #      
#                                   #
#####################################

def dm_derivative(initial_pameter,shift_position,num_qub,num_l, ent_gate, shift,kk,val_g,shots,gates_tot_dm2,G_real_dm2,noise,cps,k): 


    # Setup qubit register
    q_reg_size = num_qub #numbers of qubit (sys + det)
    q_reg_size_c= num_qub #numbers of classic bit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="DM")

    #quantum circuit: "DM (manual) for gradient"
    dm_gradient_circuit(bc, shift_position,num_qub,num_l,val_g,shift,kk,ent_gate)
    

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
    #    q_counter(gates_tot_dm2,circ.decompose().decompose().decompose())

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
# Gradient calculation DM 
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



#---------------------------------------------------------------------------------------------
#####################################
#                                   #
#        Hessian QNDM               #      
#                                   #
#####################################

#//////////#
#   main   #
#//////////#

def qndm_derivative_hessian(p_deco, initial_pameter,shift_position,sh2,pm,num_qub,num_l,ent_gate,shift,G_real_qndm,H_real_qndm,gates_tot_qndm_hess,shots,noise,val_g,simp, gradient_calc): 

    # Setup qubit register
    q_reg_size = num_qub + 1 #numbers of qubit
    q_reg_size_c = 1 #numbers of classic bit
    detect_index=num_qub #index number of the detector qubit


    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")

    if gradient_calc == True:
        #gradient
        bc = QuantumCircuit(q_reg, c_reg, name="QNDM")
        #quantum circuit: "QNDM for gradient"
        qndm_hessian_circuit(bc, shift_position,pm,num_qub,num_l,val_g,detect_index,shift, simp,ent_gate)
        #parameters initialization
        initial_values = [p_deco/2,p_deco/2]

        for i in range(len(initial_pameter)):
            initial_values.append(initial_pameter[i])

            param_dict = dict(zip(bc.parameters, initial_values))
            circ=bc.bind_parameters(param_dict)

            # measure the detector qubit
            circ.measure(detect_index, 0)

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

            p0,p1 = 0,0

            #extract counts
            for l in data.keys():
                if l == '0': # counts for detector in zero state
                    p0 += data[l]/shots # probability of |0> in the detector state
                    
                elif l == '1': # counts for detector in one state
                    p1 += data[l]/shots # probability of |1> in the detector state
        
            #calculate value of hessian   
            G_real_qndm[shift_position] = asin(2*p1-1)/(2*sin(shift)*p_deco)



    #hessian
    bc_hess = QuantumCircuit(q_reg, c_reg, name="QNDM_hess")
    #quantum circuit: "QNDM for hessian"
    qndm_hessian_circuit(bc_hess, shift_position,sh2,pm,num_qub,num_l,detect_index,shift,val_g, simp,ent_gate)
    
    #parameters initialization
    initial_values_hess = [p_deco/2,p_deco/2,p_deco/2,p_deco/2]

    for i in range(len(val_g)):
      initial_values_hess.append(initial_pameter[i])

    param_dict = dict(zip(bc_hess.parameters, initial_values_hess))
    circ_hess=bc_hess.bind_parameters(param_dict)

        
    # measure the detector qubit
    circ_hess.measure(detect_index, 0)


    # Quantum gates Counter
    #q_counter(gates_tot_qndm_hess,circ_hess.decompose().decompose().decompose())

    #print    
    #print(circ_hess) 
    #print(circ_hess.decompose().decompose())


    #simulated NOISE
   
    backend = BasicAer.get_backend('qasm_simulator')
    coupling_map = None
    basis_gates = None

    if noise == True:
        backend = FakeJakarta()
        coupling_map = backend.configuration().coupling_map
        noise_model = NoiseModel.from_backend(backend)
        basis_gates = noise_model.basis_gates



  
    job = execute(circ_hess, backend=backend, shots=shots)
    result = job.result()
    data = result.get_counts(circ_hess)   

    p0_hess,p1_hess = 0,0

    #extract counts
    for l in data.keys():
        if l == '0': # counts for detector in zero state
           p0_hess += data[l]/shots # probability of |0> in the detector state
             
        elif l == '1': # counts for detector in one state
           p1_hess += data[l]/shots # probability of |1> in the detector state
              
    #calculate value of hessian   
    H_real_qndm[shift_position,sh2] = asin(2*p1_hess-1)/(4*sin(shift)**2*p_deco)
    

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation QNDM
def qndm_hessian(lambda1, in_par ,G_real_qndm, H_real_qndm, newspop, num_qub, num_l, ent_gate, shift,gates_tot_qndm_hess,noise,shots,val_g,simp = False, gradient_calc = False, deri = False):
    """Calculate Second order derivatives with QNDM method. \n

    #Paramenters: \n
    lambda1 -- Value of coupling paramenter
    in_par -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the first derivatives information \n
    H_real_qndm -- Empty array where the function put the second derivatives information \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    gates_tot_qndm_hess -- Empty array where the function put the cost of derivatives information \n
    noise -- If noise active noise == True \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    simp -- If True the sim method is active \n
    gradient_calc -- If True the function given in output also the gradient
    deri = False -- If True the function calcualte only the first element of gradient \n
    """
 
    for sh1 in range(len(val_g)):
        for sh2 in range(len(val_g)):
            qndm_derivative_hessian(lambda1, in_par,sh1,sh2,newspop,num_qub,num_l,ent_gate,shift,G_real_qndm,H_real_qndm,gates_tot_qndm_hess,shots,noise,val_g,simp, gradient_calc)
            

#---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################
#                                   #
#            Hessian DM             #      
#                                   #
#####################################

#//////////#
#   main   #
#//////////#

def dm_derivative_hessian(initial_pameter,sh1,sh2,num_qub,num_l,ent_gate,shift1,shift2,kk,val_g,shots,gates_tot_dm2,G_real_dm, H_real_dm2,noise,cps,k, gradient_calc): 


    # Setup qubit register
    q_reg_size = num_qub #numbers of qubit (sys + det)
    q_reg_size_c= num_qub #numbers of classic bit
    #gradient

    
    if gradient_calc == True and sh2 == 0 :

        #quantum circuit
        q_reg = QuantumRegister(q_reg_size, "q")
        c_reg = ClassicalRegister(q_reg_size_c, "c")
        bc = QuantumCircuit(q_reg, c_reg, name="DM")

        #quantum circuit: "DM for gradient"
        dm_gradient_circuit(bc, sh1,num_qub,num_l,val_g,shift1,kk,ent_gate)
        

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
        #print(circ)
        #print(circ.decompose().decompose().decompose())
      

            #simulated NOISE
        backend = Aer.get_backend('aer_simulator_stabilizer')
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
        #here there is a division by four, this is correct because I am froce to iterate the gradient 2 times
        G_real_dm[sh1] += cps[k]*mean_val/(4*sin(shift1)) # probability of |0> in the system state
  

    #hessian

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc_hess = QuantumCircuit(q_reg, c_reg, name="DM")

    #quantum circuit: "DM (manual) for gradient"

    dm_hessian_circuit(bc_hess,sh1,sh2,num_qub,num_l,val_g,shift1,shift2,kk,ent_gate)

    #parameters initialization
    initial_values = []
    for i in range(len(initial_pameter)):
      initial_values.append(initial_pameter[i])

    param_dict = dict(zip(bc_hess.parameters, initial_values))
    circ_hess=bc_hess.bind_parameters(param_dict)
        

    # measure the system qubits
    qubit_index = []
    qubit_index2 = []

    for i in range(num_qub):
      qubit_index.append(num_qub-i-1)
      qubit_index2.append(i)


    circ_hess.measure(qubit_index, qubit_index2)
    #print(circ_hess)
    #print(circ_hess.decompose().decompose().decompose())
    # Quantum gates Counter
   
    #q_counter(gates_tot_dm2,circ_hess.decompose().decompose().decompose())

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
    job = execute(circ_hess, backend=backend, shots=shots)
    result = job.result()
    data = result.get_counts(circ_hess)

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
    if shift1 <0 and shift2 <0:
        H_real_dm2[sh1,sh2] -= cps[k]*mean_val/(4*(sin(shift1))**2) # probability of |0> in the system state

    if shift1 <0 and shift2 >0:
        H_real_dm2[sh1,sh2] += cps[k]*mean_val/(4*(sin(shift1))**2) # probability of |0> in the system state

    if shift1 >0 and shift2 <0:
        H_real_dm2[sh1,sh2] += cps[k]*mean_val/(4*(sin(shift1))**2) # probability of |0> in the system state

    if shift1 >0 and shift2 >0:    
        H_real_dm2[sh1,sh2] -= cps[k]*mean_val/(4*(sin(shift1))**2) # probability of |0> in the system state
   
    return
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation QNDM
def dm_hessian(in_par ,G_real_dm, H_real_dm2, PS, num_qub, num_l,ent_gate, shift,gates_tot_dm2,noise,shots,val_g,cps, gradient_calc = False, deri = False):

    """Calculate Second order derivatives with DM method. \n

    #Paramenters: \n

    in_par -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the first derivatives information \n
    H_real_qndm -- Empty array where the function put the second derivatives information \n
    PS -- Pauli string \n 
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    gates_tot_dm_hess -- Empty array where the function put the cost of derivatives information \n
    noise -- If noise active noise == True \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    cps -- Coefficient Pauli strings \n
    gradient_calc -- If True the function given in output also the gradient
    deri = False -- If True the function calcualte only the first element of gradient \n"""


    for sh1 in range(len(val_g)):
        for sh2 in range(len(val_g)):
            for k,kk in enumerate(PS):
                #sign of shift
                for shift_sign in range(4):
                    shift1 = shift
                    shift2 = shift   

                    if shift_sign == 1:
                        shift1 = -shift

                    elif shift_sign == 2:
                        shift2 = -shift

                    elif shift_sign == 3:
                        shift1 = -shift
                        shift2 = -shift


                    dm_derivative_hessian(in_par,sh1,sh2,num_qub,num_l,ent_gate,shift1,shift2,kk,val_g,shots,gates_tot_dm2,G_real_dm, H_real_dm2,noise,cps,k, gradient_calc)