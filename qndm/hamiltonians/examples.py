#Library of hamiltonians taken from [arxiv:number]
#for the QNDM optimization code.
#Written by: S. Caletti and G. Minuto
#Contacts: simone.caletti@ge.infn.it

#-------------------------------------------------------------------------------

from qiskit.quantum_info import SparsePauliOp

#-------------------------------------------------------------------------------

def get_SparsePauliOp(PS, cps): #create the SparsePauliOp from a the list
    spop = []
    for pauli, coeff in zip(PS, cps):
        couple = (pauli, coeff)
        spop.append(couple)
    return SparsePauliOp.from_list(spop)

def add_detector(PS, cps): #add the detector to create hamiltonian for qndm
    newPS = []
    for pauli in PS:
        pauli += "Z"
        newPS.append(pauli)
    return newPS, cps

def get_hamiltonian(n_qub, sel): #gerate a random hamitonian op 


            import random

            # Specify the values to include in the string
            values = ['X', 'Y', 'Z']
            PS = []
            cps = []

            #gauss info
            mu = 1
            sigma = 0.5

            # Set the length of the random string
            string_length = n_qub  # Change this to your desired length

            # Generate a sum of random string 
            for _ in range(int(sel)):

                random_string = ''.join(random.choice(values) for _ in range(string_length))
                prob_string = random.gauss(mu, sigma)

                PS.append(random_string)
                cps.append(prob_string)
              

            return PS, cps
