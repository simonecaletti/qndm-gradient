#!/usr/bin/python3


from math import sin
from qiskit import QuantumCircuit, execute, BasicAer
from qiskit import QuantumRegister, ClassicalRegister, execute,QuantumCircuit
from qndm.derivatives.gradient.qndm import qndm_gradient_circuit
from qndm.derivatives.gradient.dm import dm_gradient_circuit
from qndm.derivatives.hessian.qndm import qndm_hessian_circuit
from qndm.derivatives.hessian.dm import dm_hessian_circuit
from qndm.hamiltonians.normalization.lam_balancing import get_lambda_balancing

import numpy as np
from math import asin



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
def qndm_derivative(lambda1, initial_pameter,shift_position,newspop,num_qub,num_l, ent_gate,shift,G_real_qndm,shots,val_g): 

    # Setup qubit register
    q_reg_size = num_qub + 1 #numbers of qubit (sys + det)
    q_reg_size_c = 1 #numbers of classic bit
    detect_index = num_qub #index number of the detector qubit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="QNDM")

    #balacing of lambda in function of hamiltonian
    lambda1 = get_lambda_balancing(newspop,lambda1)  

    #quantum circuit: "QNDM for gradient"
    qndm_gradient_circuit(bc, shift_position,newspop,num_qub,num_l,val_g,detect_index,shift,ent_gate)
    
    #parameters vector
    initial_values = [lambda1/2,lambda1/2]
    for i in range(len(initial_pameter)):
        initial_values.append(initial_pameter[i])

    #parameters initialization
    param_dict = dict(zip(bc.parameters, initial_values))
    circ=bc.bind_parameters(param_dict)

    # measure the detector qubit
    circ.measure(detect_index, 0)
 
    #backend
    backend = BasicAer.get_backend('qasm_simulator')

    #run quantum circuit 
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

    #derivative in the direction e_(shift_position)
    G_real_qndm[shift_position] = asin(2*p1-1)/(2*(lambda1))

    return None 

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Gradient calculation QNDM
def qndm_gradient(lambda1, pars ,G_real_qndm, newspop, num_qub, num_l, ent_gate,shift,shots,val_g):

    """Calculate first order derivatives with QNDM method. \n

    #Paramenters: \n
    lambda1 -- Value of coupling paramenter
    pars -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the derivatives information \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    """
    
    for shift_position in range(len(val_g)):
        qndm_derivative(lambda1, pars ,shift_position, newspop, num_qub, num_l, ent_gate, shift, G_real_qndm,shots,val_g)


#---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################
#                                   #
#            Gradient DM            #      
#                                   #
#####################################

def dm_derivative(initial_pameter,shift_position,num_qub,num_l, ent_gate, shift,kk,val_g,shots,G_real_dm,cps,k): 

    # Setup qubit register
    q_reg_size = num_qub #numbers of qubit
    q_reg_size_c= num_qub #numbers of classic bit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc = QuantumCircuit(q_reg, c_reg, name="DM")

    #quantum circuit: "DM (manual) for gradient"
    dm_gradient_circuit(bc, shift_position,num_qub,num_l,val_g,shift,kk,ent_gate)
    
    #parameters vector
    initial_values = []
    for i in range(len(initial_pameter)):
      initial_values.append(initial_pameter[i])

    #parameters initialization
    param_dict = dict(zip(bc.parameters, initial_values))
    circ=bc.bind_parameters(param_dict) 

    #measuring lists
    qubit_index = []
    qubit_index2 = []
    for i in range(num_qub):
      qubit_index.append(num_qub-i-1)
      qubit_index2.append(i)

    # measure the system qubits
    circ.measure(qubit_index, qubit_index2)

    #backend
    backend = BasicAer.get_backend('qasm_simulator')
 
    #run quantum circuit
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

    #estimation cost function for a single pauli string
    mean_val = -plus+minus

    #derivative in the direction e_(shift_position) 
    G_real_dm[shift_position] += np.real(cps[k])*mean_val/(2*sin(shift)) # probability of |0> in the system state
    
    return None 


#---------------------------------------------------------------------------------------------------------------------------------------------------
# Gradient calculation DM 
def dm_gradient(pars ,G_real_dm, spop, num_qub, num_l, ent_gate,shift,shots,val_g):
    """Calculate first order derivatives with Direct Measurament method. \n

    #Paramenters: \n
    pars -- Parameters of the rotational gates \n
    G_real_dm -- Empty array where the function put the derivatives information \n
    spop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    """            
    
    for shift_position in range(len(val_g)):
        for k,kk in enumerate(spop.paulis):

            for shift_sign in range(2):
                if shift_sign == 1:
                    shift *=-1 
                        
                dm_derivative(pars,shift_position, num_qub, num_l, ent_gate,shift, kk,val_g,shots,G_real_dm,spop.coeffs,k)
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

def qndm_derivative_hessian(lambda1, initial_pameter,sh1,sh2,newspop,num_qub,num_l,ent_gate,shift,G_real_qndm,H_real_qndm,shots,val_g, gradient_calc): 

    # Setup qubit register
    q_reg_size = num_qub + 1 #numbers of qubit (sys + det)
    q_reg_size_c = 1 #numbers of classic bit
    detect_index=num_qub #index number of the detector qubit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")

    #balacing of lambda in function of hamiltonian
    lambda1 = get_lambda_balancing(newspop,lambda1) 

    #hessian
    bc_hess = QuantumCircuit(q_reg, c_reg, name="QNDM_hess")

    #quantum circuit: "QNDM for hessian"
    qndm_hessian_circuit(bc_hess, sh1,sh2,newspop,num_qub,num_l,detect_index,shift,val_g,ent_gate)
    
    #parameters vector
    initial_values_hess = [lambda1/2,lambda1/2,lambda1/2,lambda1/2]
    for i in range(len(val_g)):
      initial_values_hess.append(initial_pameter[i])

    #parameters initialization
    param_dict = dict(zip(bc_hess.parameters, initial_values_hess))
    circ_hess=bc_hess.bind_parameters(param_dict)

    # measure the detector qubit
    circ_hess.measure(detect_index, 0)

    #backend
    backend = BasicAer.get_backend('qasm_simulator')

    #run quantum circuit
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
              
    #hessian in the direction e_(sh1) e e_(sh2)   
    H_real_qndm[sh1,sh2] = asin(2*p1_hess-1)/(4*sin(shift)**2*lambda1)

    #if true the code calculates also the gradient
    if gradient_calc == True:
        qndm_derivative(lambda1, initial_pameter ,sh1, newspop, num_qub, num_l, ent_gate, shift, G_real_qndm,shots,val_g)
        

    return

#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation QNDM
def qndm_hessian(lambda1, pars ,G_real_qndm, H_real_qndm, newspop, num_qub, num_l, ent_gate, shift,shots,val_g, gradient_calc = False):
    """Calculate Second order derivatives with QNDM method. \n

    #Paramenters: \n
    lambda1 -- Value of coupling paramenter
    pars -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the first derivatives information \n
    H_real_qndm -- Empty array where the function put the second derivatives information \n
    newspop -- Hamiltonian \n
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    ent_gate -- Type of entanglement gate (0=CNOT, 1=SWAP) \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    gradient_calc -- If True the function given in output also the gradient
    """
 
    for sh1 in range(len(val_g)):
        for sh2 in range(len(val_g)):
            qndm_derivative_hessian(lambda1, pars,sh1,sh2,newspop,num_qub,num_l,ent_gate,shift,G_real_qndm,H_real_qndm,shots,val_g, gradient_calc)
            

#---------------------------------------------------------------------------------------------------------------------------------------------------
#####################################
#                                   #
#            Hessian DM             #      
#                                   #
#####################################

#//////////#
#   main   #
#//////////#

def dm_derivative_hessian(initial_pameter,sh1,sh2,num_qub,num_l,ent_gate,shift1,shift2,kk,val_g,shots,G_real_dm, H_real_dm,cps,k, gradient_calc): 


    # Setup qubit register
    q_reg_size = num_qub #numbers of qubit 
    q_reg_size_c= num_qub #numbers of classic bit

    #quantum circuit
    q_reg = QuantumRegister(q_reg_size, "q")
    c_reg = ClassicalRegister(q_reg_size_c, "c")
    bc_hess = QuantumCircuit(q_reg, c_reg, name="DM")

    #quantum circuit: "DM (manual) for gradient"
    dm_hessian_circuit(bc_hess,sh1,sh2,num_qub,num_l,val_g,shift1,shift2,kk,ent_gate)

    #paramenters vector
    initial_values = []
    for i in range(len(initial_pameter)):
      initial_values.append(initial_pameter[i])

    #parameters initialization
    param_dict = dict(zip(bc_hess.parameters, initial_values))
    circ_hess=bc_hess.bind_parameters(param_dict)

    
    #measuring lists
    qubit_index = []
    qubit_index2 = []
    for i in range(num_qub):
      qubit_index.append(num_qub-i-1)
      qubit_index2.append(i)

    # measure the system qubits
    circ_hess.measure(qubit_index, qubit_index2)

    #backend
    backend = BasicAer.get_backend('qasm_simulator')
 
    #run quantum circuit 
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

    #estimation cost function for a single pauli string
    mean_val = -plus+minus

    #hessian in the direction e_(sh1) e e_(sh2)   
    if shift1 <0 and shift2 <0:
        H_real_dm[sh1,sh2] -= np.real(cps[k]*mean_val/(4*(sin(shift1))**2)) # probability of |0> in the system state

    if shift1 <0 and shift2 >0:
        H_real_dm[sh1,sh2] += np.real(cps[k]*mean_val/(4*(sin(shift1))**2)) # probability of |0> in the system state

    if shift1 >0 and shift2 <0:
        H_real_dm[sh1,sh2] += np.real(cps[k]*mean_val/(4*(sin(shift1))**2)) # probability of |0> in the system state

    if shift1 >0 and shift2 >0:    
        H_real_dm[sh1,sh2] -= np.real(cps[k]*mean_val/(4*(sin(shift1))**2)) # probability of |0> in the system state
   

    #if true the code calculates also the gradient
    if gradient_calc == True and sh2 == 0 :

        dm_derivative(initial_pameter,sh1, num_qub, num_l, ent_gate,shift1, kk,val_g,shots,G_real_dm,cps/2,k)

    return
#---------------------------------------------------------------------------------------------------------------------------------------------------
# Hessian calculation QNDM
def dm_hessian(pars ,G_real_dm, H_real_dm, spop, num_qub, num_l,ent_gate, shift,shots,val_g, gradient_calc = False):
    """Calculate Second order derivatives with DM method. \n

    #Paramenters: \n

    pars -- Parameters of the rotational gates \n
    G_real_qndm -- Empty array where the function put the first derivatives information \n
    H_real_qndm -- Empty array where the function put the second derivatives information \n
    PS -- Pauli string \n 
    num_qub -- Number of Qubits \n
    num_l -- Number of Layer \n
    shift -- Value of shit of 'Paramenter shift rule' \n
    shots -- Number of Shots \n
    val_g -- serial description of the rotational gates  \n
    cps -- Coefficient Pauli strings \n
    gradient_calc -- If True the function given in output also the gradient"""

    for sh1 in range(len(val_g)):
        for sh2 in range(len(val_g)):
            for k,kk in enumerate(spop.paulis):
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


                    dm_derivative_hessian(pars,sh1,sh2,num_qub,num_l,ent_gate,shift1,shift2,kk,val_g,shots,G_real_dm, H_real_dm,spop.coeffs,k, gradient_calc)
