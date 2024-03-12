#!/usr/bin/python3

#Code for calculation second derivatives:  DM method
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.


#--------------------------------------------------------------------------------------------


import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, sin
from qiskit import Aer
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeLondonV2, FakeManilaV2, FakeJakarta


#---------------------------------------------------------------------------------------------
#import QNDM package


from qndm.hamiltonians.examples import add_detector, get_SparsePauliOp
from qndm.core import *

from qndm.hamiltonians.examples import get_hamiltonian
from qndm.hamiltonians.hydrogen import get_model
from qndm.tools.error import get_qndm_error_second, get_dm_error_second
from qndm.tools.runcard import print_runcard


#---------------------------------------------------------------------------------------------

############################
#                          #
#      Hamiltonian M       #
#                          #
############################

hamlib_ = False

if hamlib_ == True:
    # Input to use with the qndm.hamiltonians.hydrogen package
    n = 2 # atoms count
    shape = 'linear' # linear, pyramid, ring, sheet
    r = 0.6 # between 0.5 and 2.0 (step 0.1)
    key = "/ham_BK/"
    
    #get the model (we need for first because sets nqubit)
    model = get_model(n, shape, r, key)

    PS = model["PS"]
    cps = model["cps"]
    pauli_string = len(cps)
    num_qub = model["nqubit"]
  

else:

    #number of qubit of the quantum register
    num_qub = 2

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
num_l = 4

#inside a layer: number of rotational layers
lay_u = 1

#entanglement layer
ent_gate = 0
# if ent_gate = 0 ---> CNOT
# if ent_gate = 1 ---> SWAP

#total number of parameters per qubit
n_pars=lay_u*num_l*num_qub 

#Rotational array: here there are the gates information to implent unitary trasformation U
#code: rx = 1, ry = 2, rz = 3
val_g = np.random.randint(1, 4, size=n_pars) #val_g = [1,1,2,2,3,3]

#Parameters array: here there are the parameters information for each gates in U
pars = np.random.rand(n_pars)



#############################
#                           #
#  Derivative Paramenters   #
#                           #
#############################


# Input to use with the qndm.hamiltonians.example package

shots = 10000 #number shots for a single evaluation

#shift (paramenter shift rule)
shift = pi/2 

#noise
noise = False 
#if noise = False --> Simulator = FakeSimulator
#if noise = True --> Simulator = Aer


#--------------------------------------------------------------------------------------------

#R U N - C A R D#

print_runcard(num_qub, num_l, val_g, newspop, shots, ent_gate=0, noise=False, output_path="./output_test")

#------------------------------------------------------------------

print("Into the derivatives process...", end="")


#Gradient DM
G_real_dm = np.zeros(n_pars)
#Hessian DM
H_real_dm = np.zeros((n_pars,n_pars))
   

dm_hessian(pars ,G_real_dm, H_real_dm, spop, num_qub, num_l, ent_gate,shift,noise,shots,val_g)

    

#error calculation:
error_DM = get_dm_error_second(n_pars, spop, shots, shift)

print(" done!")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


####  W R I T I N G - O U T P U T #####


#reshaping


#dataframes
H_real_dm = H_real_dm.ravel()
error_DM = error_DM.ravel()

#DM dataframe
DM_data = {
    'Derivatives': H_real_dm,
    'Errors': error_DM
}
df_DM = pd.DataFrame(DM_data)

output_path = "./output_test"
# Save the dataframes to a CSV files
DM_path = os.path.join(output_path,'DM_der_second.csv' )

df_DM.to_csv(DM_path, index=False)