#!/usr/bin/python3

#Code for gradient evaluation
#using the QNDM method.
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.

#--------------------------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from math import pi



#---------------------------------------------------------------------------------------------
#import QNDM package

from qndm.hamiltonians.examples import add_detector, get_SparsePauliOp
from qndm.core import *
from qndm.hamiltonians.examples import get_hamiltonian
from qndm.hamiltonians.hydrogen import get_model
from qndm.utils.error import get_qndm_error



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
    PS, cps = get_hamiltonian(num_qub, pauli_string, sel=10)


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

shots = 50 #number shots for a single evaluation

#shift (paramenter shift rule)
shift = pi/2 

#coupling parameter QNDM
lambda1 = 0.1 

#--------------------------------------------------------------------------------------------
#R U N - C A R D#

#print_run_card(output_dir="./output_test", n_qubits = num_qub, n_layers = num_l, val_g = val_g, spop =  newspop, n_shots = shots, lambda1 = lambda1, ent_gate=0,)
#------------------------------------------------------------------


print("Into the derivatives process...", end="")


#array to save gradient results 
G_real_qndm = np.zeros(n_pars)

#gate counter
gates_tot_qndm=np.zeros(12)

#gradient with qndm
if __name__ == '__main__':
    G_real_qndm = qndm_gradient(lambda1=lambda1, 
              pars=pars,
              newspop = newspop,
              num_qub = num_qub,
              num_l = num_l,
              ent_gate = ent_gate,
              shift = shift,            
              shots = shots,
              val_g = val_g)

    #error calculation:
error_QNDM = get_qndm_error(pars, G_real_qndm, lambda1, shots, shift)


print(" done!")



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

####  W R I T I N G - O U T P U T #####


#QNDM dataframe
QNDM_data = {
    'Parameters': pars,
    'Derivatives': G_real_qndm,
    'Errors': error_QNDM
}
df_QNDM = pd.DataFrame(QNDM_data)


# Save the dataframes to a CSV files
output_path = "./output_test"
QNDM_path = os.path.join(output_path,'QNDM_der.csv' )


df_QNDM.to_csv(QNDM_path, index=False)
 
