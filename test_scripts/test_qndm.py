#!/usr/bin/python3

#Code for gradient evaluation
#using the QNDM method.
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
sys.path.insert(1, '/Users/giovanniminuto/Documents/My_codes/QNDM_release/qndm-gradient')

from qndm.hamiltonians.examples import add_detector, get_SparsePauliOp
from qndm.core import *

from qndm.hamiltonians.examples import get_hamiltonian
from qndm.hamiltonians.hydrogen import get_model
from qndm.tools.error import get_qndm_error
from qndm.tools.runcard import print_runcard



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
lay_u = 2

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
cas = np.random.rand(pars)




#############################
#                           #
#  Derivative Paramenters   #
#                           #
#############################


# Input to use with the qndm.hamiltonians.example package

shots = 50000 #number shots for a single evaluation

#shift (paramenter shift rule)
shift = pi/2 

#noise
noise = False 
#if noise = False --> Simulator = FakeSimulator
#if noise = True --> Simulator = Aer

#coupling parameter QNDM
lambda1 = 0.1 


#--------------------------------------------------------------------------------------------
#R U N - C A R D#

print_runcard(num_qub, num_l, val_g, newspop, shots, lambda1, ent_gate=0, noise=False, output_path="./output_test")
#------------------------------------------------------------------

print("Into the derivatives process...", end="")


#array to save gradient results 
G_real_qndm = np.zeros(pars)

#gate counter
gates_tot_qndm=np.zeros(12)

#gradient with qndm
qndm_gradient(lambda1, cas ,G_real_qndm, newspop, num_qub, num_l, ent_gate, shift,gates_tot_qndm,noise,shots,val_g)

#error calculation:
error_QNDM = get_qndm_error(pars, G_real_qndm, lambda1, shots, shift)


print(" done!")



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####  W R I T I N G - O U T P U T #####


#QNDM dataframe
QNDM_data = {
    'Parameters': cas,
    'Derivatives': G_real_qndm,
    'Errors': error_QNDM
}
df_QNDM = pd.DataFrame(QNDM_data)


# Save the dataframes to a CSV files
output_path = "./output_test"
QNDM_path = os.path.join(output_path,'QNDM_der_0.csv' )


df_QNDM.to_csv(QNDM_path, index=False)
 
