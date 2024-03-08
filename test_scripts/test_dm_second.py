#!/usr/bin/python3

#Code for calculation second derivatives:  DM method
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.


#--------------------------------------------------------------------------------------------

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi,sin,sqrt
from qiskit import Aer
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeLondonV2, FakeManilaV2, FakeJakarta

#---------------------------------------------------------------------------------------------
#import QNDM package

from qndm.hamiltonians.examples import add_detector, get_SparsePauliOp
from qndm.core import *

from qndm.hamiltonians.examples import get_hamiltonian
from qndm.hamiltonians.hydrogen import get_model



#---------------------------------------------------------------------------------------------

#number of qubit of the quantum register
num_qub = 2

############################
#                          #
#      Hamiltonian M       #
#                          #
############################

hamlib_ = False

if hamlib_ == True:
    # Input to use with the qndm.hamiltonians.hydrogen package
    shape = 'linear' # linear, pyramid, ring, sheet
    r = 0.6 # between 0.5 and 2.0 (step 0.1)
    key = "/ham_BK/"
    
    #get the model (we need for first because sets nqubit)
    model = get_model(num_qub, shape, r, key)

    PS = model["PS"]
    cps = model["cps"]
    pauli_string = len(cps)
    num_qub = model["nqubit"]
  

else:

    #select the pauli string number
    pauli_string = 2

    #hamiltonians M
    PS, cps = get_hamiltonian(num_qub, pauli_string)


spop = get_SparsePauliOp(PS, cps) #spop = sparse pauli operator

#hamiltonians for QNDM: here we add the detector operator equal to Z
PS_QNDM, cps_QNDM = add_detector(PS, cps)
newspop = get_SparsePauliOp(PS_QNDM, cps_QNDM) #After adding the detector

############################
#                          #
#     Quantum Circuit      #
#                          #
############################


#layer = rotational layers + entanglement layer
#number of layers = rotational layers + entanglement layer
num_l = 2

#inside a layer: number of rotational layers
lay_u = 1

#entanglement layer
ent_gate = 0
# if ent_gate = 0 ---> CNOT
# if ent_gate = 1 ---> SWAP

#total number of parameters per qubit
pars=lay_u*num_l*num_qub 

#Rotational array: here there are the gates information to implent unitary trasformation U
#code: rx = 1, ry = 2, rz = 3
val_g = np.random.randint(1, 4, size=pars) #val_g = [1,1,2,2,3,3]
val_g = np.array(val_g)

#Parameters array: here there are the parameters information for each gates in U
lon=10
cas = np.random.rand(lon,pars)




#############################
#                           #
#  Derivative Paramenters   #
#                           #
#############################


# Input to use with the qndm.hamiltonians.example package

shots = 500 #number shots for a single evaluation

#shift (paramenter shift rule)
shift = pi/2 

#noise
noise = False 
#if noise = False --> Simulator = FakeSimulator
#if noise = True --> Simulator = Aer

#coupling parameter QNDM
lambda1 = 0.1 


#--------------------------------------------------------------------------------------------


#####R U N - C A R D#####

print(f"Writing the RunCard...", end="")


# Save the data to a .txt file inside the directory
output_path = f'./output_path'
run = os.path.join(output_path,'RunCard_Der_2.txt' )

with open(run, "w") as f:

    f.write(" -------------------- \n")
    f.write("|                    |\n")
    f.write("|  RUNCARD (Der run) |\n")
    f.write("|                    |\n")
    f.write(" -------------------- \n")
    f.write("\n")

    #number of qubits
    f.write("Number of qubits = {} \n".format(num_qub))
    f.write("\n")

    #number of layers 
    f.write("Layers = {} \n".format(num_l))

    f.write("Rotation U array: {} \n".format(val_g))
    f.write('Info: rx = 1, ry = 2, rz = 3 \n')

    #entanglement gates
    # if ent_gate = 0 ---> CNOT
    # if ent_gate = 1 ---> SWAP
    if ent_gate == 0:
        f.write("Entaglment Gates are CNOTS\n")
    if ent_gate == 1:
        f.write("Entaglment Gates are SWAPS\n")

    #shots
    f.write("Shots = {} \n".format(shots))


    f.write("Hamiltonian {} \n".format(spop))
    f.write("\n")
    f.write("\n")

    #shift (paramenter shift rule)
    f.write("Shift of parameter shift rule: {}\n".format(shift))



    #simulated NOISE
    backend = Aer.get_backend('aer_simulator_stabilizer')
    coupling_map = None
    basis_gates = None
    noise_model = None

    if noise == True:
        backend = FakeJakarta()
        coupling_map = backend.configuration().coupling_map
        noise_model = NoiseModel.from_backend(backend)
        basis_gates = noise_model.basis_gates

    # Select the QasmSimulator from the Aer provider

    #simulator = FakeJakartaV2()
    f.write("noise model: {} \n".format(noise_model))
    f.write("basis gates: {}\n".format(basis_gates))
    f.write("couplingmap: {}\n".format(coupling_map))


    f.close()
print("done!")

#------------------------------------------------------------------

print("Into the derivatives process...", end="")


#array to save hessian results 
gd_med_DM = np.zeros((lon,pars))
hess_med_DM2 = np.zeros((lon,pars,pars))

#gate counter
gates_tot_dm2 = np.zeros(12)


for j,in_par in enumerate(cas):

    #Gradient DM
    G_real_dm = np.zeros(pars)
    #Hessian DM
    H_real_dm2 = np.zeros((pars,pars))
   

    dm_hessian(in_par ,G_real_dm, H_real_dm2, PS, num_qub, num_l, ent_gate,shift,gates_tot_dm2,noise,shots,val_g,cps)
    
   
    gd_med_DM[j,:] = G_real_dm
    hess_med_DM2[j,:,:] = H_real_dm2
    

#error calculation:
error_DM = np.zeros((lon,pars,pars))
error_DM[:,:,:] = (1/(sqrt(shots)*sin(shift)))*sum(np.real(np.asarray(spop.coeffs)))


print(" done!")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


####  W R I T I N G - O U T P U T #####


#reshaping
hess_DM = hess_med_DM2.reshape((10, -1))

#dataframes
df_DM = pd.DataFrame(hess_DM)
data_input = pd.DataFrame(cas)


for i in range(pars*pars):
    df_DM = df_DM.rename(columns={i: f'der_{i}'})


df_DM = pd.concat([data_input, df_DM], axis = 1)


# Adding a column with quantum info
info_qs = np.array([num_qub,num_l,pauli_string,shots,lambda1,lon])
data_col = np.zeros(lon)
data_col[:len(info_qs)] = info_qs

df_DM['info'] = data_col


# Save the dataframes to a CSV files
file_path_DM = os.path.join(output_path, f"DM_der_second.csv")

df_DM.to_csv(file_path_DM, index=False)
 
