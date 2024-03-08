#!/usr/bin/python3

#Code for gradient evaluation
#using the QNDM method.
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.

import os
import sys
import pandas as pd
import numpy as np

#####R U N - C A R D#####

def print_runcard(num_qub, num_l, val_g, ent_gate=0, shots, shift="np.pi/2", lambda1, noise=False, output_path="./output_test"):
    print("Writing the RunCard...", end="")

    # Save the data to a .txt file inside the directory
    run = os.path.join(output_path,'RunCard_Der.txt' )

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
            f.write("Entaglment Gates are CNOTS \n")
        if ent_gate == 1:
            f.write("Entaglment Gates are SWAPS \n")
    
        #shots
        f.write("Shots = {} \n".format(shots))
    

        f.write("Hamiltonian {} \n".format(spop))
        f.write("\n")
        f.write("\n")

        #shift (paramenter shift rule)
        f.write("Shift of parameter shift rule: {}\n".format(shift))

        # fixed lambda QNDM
        f.write("Lambda (QNDM coupling) = {}\n".format(lambda1))


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
