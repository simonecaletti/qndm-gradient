#Library of hamiltonians taken from [arxiv:number]
#for the QNDM optimization code.
#Written by: S. Caletti and G. Minuto
#Contacts: simone.caletti@ge.infn.it

#-------------------------------------------------------------------------------

from qiskit.quantum_info import SparsePauliOp
import random
import numpy as np
from qiskit.quantum_info import Pauli
from qiskit.quantum_info import SparsePauliOp

#-------------------------------------------------------------------------------


#==========================================================#
#
# function to generate a observable where the eigenvalues are gaussian values
#
#==========================================================#


def binary_representation(n, num_bits):
    """Return the binary representation of n with num_bits bits."""
    return format(n, '0{}b'.format(num_bits))

def calculate_spop(diagonal_elements):
    """Write the Pauli strings for a given diagonal observable."""
    num_qubits = int(np.log2(len(diagonal_elements)))


    PS,cps = [],[]

    for i in range(2**num_qubits):
        binary_i = binary_representation(i, num_qubits)
        coefficient = 0
        for k in range(2**num_qubits):
            binary_k = binary_representation(k, num_qubits)
            exponent = sum(int(binary_i[j]) * int(binary_k[j]) for j in range(num_qubits))
            coefficient += (-1)**exponent * diagonal_elements[k]

        coefficient /= 2**num_qubits
        if abs(coefficient) < 10**(-2):
            continue 
        pauli_string = ''.join('Z' if bit == '1' else 'I' for bit in binary_i)

        PS.append(pauli_string)
        cps.append(np.real(coefficient))

    

    return PS,cps

#-------------------------------------------------------------------------------

def get_SparsePauliOp(PS, cps):  # Create the SparsePauliOp from the list
    spop = []
    for pauli, coeff in zip(PS, cps):
        if isinstance(pauli, str) and isinstance(coeff, (int, float)):
            spop.append((pauli, coeff))
        else:
            raise ValueError(f"Elemento non valido in PS o cps: {pauli}, {coeff}")
    return SparsePauliOp.from_list(spop)




def add_detector(PS, cps): #add the detector to create hamiltonian for qndm
    newPS = []
    for pauli in PS:
        pauli += "Z"
        newPS.append(pauli)
    return newPS, cps






def get_hamiltonian(n_qub, n_ps, sel=None, mu=None, sigma=None):
    """
    Generates a Hamiltonian in terms of Pauli strings and coefficients.

    Parameters:
        n_qub: int
            Number of qubits.
        n_ps: int
            Number of Pauli strings.
        sel: str or int, optional
            Selection mode for generating the Hamiltonian.
            'cerezo'     -> Cerezo's method.
            'ising'      -> Modified Ising Hamiltonian.
            'ising_mod'  -> Alternative Ising Hamiltonian.
            'ising_zero' -> Ising Hamiltonian with near-zero coefficients.
            int          -> Random Observable with specified number of Pauli strings.
        mu: float, optional
            Mean value for the normal distribution (default depends on `sel`).
        sigma: float, optional
            Standard deviation for the normal distribution (default depends on `sel`).

    Returns:
        PS:  list of str   -> Pauli strings
        cps: list of float -> coefficients of the Pauli strings
    """

    PS  = []
    cps = []

    # Generate an observable as per Cerezo's method
    if sel == 'cerezo':
        identities = 'I' * n_qub
        PS.append(identities)
        cps.append(-0.5)

        for j in range(n_qub):
            z_id = ''.join('Z' if i == j else 'I' for i in range(n_qub))
            PS.append(z_id)
            cps.append(-0.5 / n_qub)

    # Generate Ising Hamiltonian variants
    elif sel in {'ising', 'ising_mod', 'ising_zero'}:
        strings = ['X', 'Y', 'Z', 'I']

        # Default values if mu and sigma are not provided
        default_params = {
            'ising': (1, 0.1),
            'ising_mod': (1, 0.1),
            'ising_zero': (0, 0.01)
        }
        mu = mu if mu is not None else default_params[sel][0]
        sigma = sigma if sigma is not None else default_params[sel][1]

        for _ in range(n_ps):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            sign = random.choice([-1, 1]) if sel == 'ising_mod' else 1
            prob_string = sign * random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)

    # Generate a random observable
    elif isinstance(sel, int):
        strings = ['X', 'Y', 'Z', 'I']
        mu = mu if mu is not None else 10
        sigma = sigma if sigma is not None else 0.1

        for _ in range(sel):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            prob_string = random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)

    else:
        raise ValueError("Invalid value for 'sel'. Must be 'cerezo', 'ising', 'ising_mod', 'ising_zero' or an integer.")

    return PS, cps




















'''

def get_hamiltonian(n_qub, n_ps, sel=None): 
    """
    Generates a Hamiltonian in terms of Pauli strings and coefficients.
    
    Parameters:
        n_qub: int
            Number of qubits.
        n_ps: int
            Number of Pauli strings.
        sel: str or int, optional
            Selection mode for generating the Hamiltonian.
            'cerezo' -> Cerezo's method.
            'ising' -> Modified Ising Hamiltonian.
            int -> Random Observable with specified number of Pauli strings.
    
    Returns:
        PS:  list of str   -> Pauli strings
        cps: list of float -> coefficients of the Pauli strings
    """

    
    PS = []
    cps = []

    # Generate an observable as per Cerezo's method
    if sel == 'cerezo':
        identities = ''.join('I' for _ in range(n_qub))
        PS.append(identities)
        cps.append(-(1/2))

        for j in range(n_qub):
            z_id = ''.join('Z' if i == j else 'I' for i in range(n_qub))
            PS.append(z_id)
            cps.append((-1/2) / n_qub)

    # Generate Ising Hamiltonian
    elif sel == 'ising':
        strings = ['X', 'Y', 'Z', 'I']
        mu = 1
        sigma = 0.1

        for _ in range(n_ps):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            sign = 1 # random.choice([-1, 1])
            prob_string = sign * random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)


    # Generate a modified Ising Hamiltonian
    elif sel == 'ising_mod':
        strings = ['X', 'Y', 'Z', 'I']
        mu = 1
        sigma = 0.1

        for _ in range(n_ps):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            sign = random.choice([-1, 1])
            prob_string = sign * random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)




    # Generate a modified Ising Hamiltonian
    elif sel == 'ising_zero':
        strings = ['X', 'Y', 'Z', 'I']
        mu = 0
        sigma = 0.01

        for _ in range(n_ps):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            sign = 1 #random.choice([-1, 1])
            prob_string = sign * random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)



    # Generate a random observable
    elif isinstance(sel, int):
        strings = ['X', 'Y', 'Z', 'I']
        mu = 10
        sigma = 0.1

        for _ in range(sel):
            random_string = ''.join(random.choice(strings) for _ in range(n_qub))
            prob_string = random.gauss(mu, sigma)
            PS.append(random_string)
            cps.append(prob_string)

    else:
        raise ValueError("Invalid value for 'sel'. Must be 'cerezo', 'ising', or an integer.")

    return PS, cps
'''






