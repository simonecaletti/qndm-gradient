#!/usr/bin/python3

#Code for gradient optimization
#using the QNDM method.
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 and 2404.02245 [quant-ph] in case you use 
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
from qiskit.primitives import Estimator
from qiskit.circuit import Parameter,ParameterVector

#---------------------------------------------------------------------------------------------
#import QNDM package

from qndm.hamiltonians.examples import add_detector, get_SparsePauliOp
from qndm.core import *

from qndm.hamiltonians.examples import get_hamiltonian
from qndm.hamiltonians.hydrogen import get_model
from qndm.tools.error import get_qndm_error
from qndm.tools.runcard import print_runcard

from qndm.layers.unitaries_gradient import *



#---------------------------------------------------------------------------------------------

#number of qubit of the quantum register
num_qub = 2

############################
#                          #
#      Hamiltonian M       #
#                          #
############################

hamlib_ = True

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
num_l = 3

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

#learning rate
lr = 0.2

#iteration number
iter = 50


#--------------------------------------------------------------------------------------------
#R U N - C A R D#

print_runcard(num_qub, num_l, val_g, newspop, shots, lambda1, ent_gate=0, noise=False, output_path="./output_test")

#to add: opt paramenters
#------------------------------------------------------------------

print("Into the optimization process...", end="")


#cost function evaluation
def evaluate_expectation(cas,num_qub, num_l,lay_u,val_g,shift,ent_gate,spop,shots):

    #circuit initialization
    circuit = QuantumCircuit(num_qub)
    params = ParameterVector("theta", length=num_qub*num_l*lay_u)

    #unitary trasformation
    l_d=U1(val_g,params,num_qub,num_l,shift,10000,ent_gate)
    qubits = []
    for i in range(num_qub):
        qubits.append(i)
    circuit.compose(l_d, qubits=qubits, inplace=True)

    #evaluation
    estimator = Estimator()
    expectation_value = estimator.run(circuit, spop, cas, shots = shots).result().values.real

    return expectation_value


############GRADIENT DESCENT QNDM###############

save_par=np.zeros((iter,pars))
save_cost=np.zeros((iter,1))

for i in range(iter):
    #parameters saving 
    save_par[i,:]=cas
    #Cost function saving 
    save_cost[i] = evaluate_expectation(cas,num_qub, num_l,lay_u,val_g,shift,ent_gate,spop,shots)

    G_real_qndm = np.zeros(pars)

    #gradient 
    qndm_gradient(lambda1, cas ,G_real_qndm, newspop, num_qub, num_l, ent_gate, shift,noise,shots,val_g)

    #gradient optimization
    cas=cas-lr*G_real_qndm

print(" done!")


#------------------------------------------------------------------

####  W R I T I N G - O U T P U T #####

print("Printing output data...", end="")

out = open("output.dat", "w")
out.write("# OUTPUT\n")
out.write("# QNDM Cost Function Minimization \n")
out.write("\n")

#cost function plot
cir_evals_qndm = np.array([i for i in range(iter)])


headline = "Circ Evals, Cost Func"
for j in range(pars):
    headline += ", Param{}".format(j)
headline += "\n"
out.write(headline)

for i in range(iter):
    line = "{}, {}".format(cir_evals_qndm[i], save_cost[i][0])
    for j in range(pars):
        line += ", {}".format(save_par[i, j])
    line += "\n"
    out.write(line)

print(" done!")
