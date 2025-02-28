#==========================================================#
#
# QISKIT code to calculate gradient with zero-order method
#==========================================================#

from qiskit import QuantumCircuit
import numpy as np
import concurrent.futures
import time
import os
import sys
from qiskit.quantum_info import Statevector
from qiskit.circuit import ParameterVector
from pprint import pprint
import pickle

#==========================================================#
#
# Functions for parameters and gates_var
#==========================================================#

def read_parameters(output_dir, read_pars : str, n_qubits : int, n_layers : int):
 
    # Parametri
    parameters = None

    if read_pars == "True":
        try:
            # Leggi i parametri dal file
            param_filename = f'{output_dir}/parameters_n_{n_qubits}_qubits.txt'
            parameters = np.loadtxt(param_filename, dtype = 'float')  # Specifica il tipo come float
            
            # Controlla che i dati possano essere rimodellati correttamente
            if parameters.size == n_layers * n_qubits:
                parameters = parameters.reshape((n_layers, n_qubits))
                print(f'Parametri caricati da {param_filename} e rimodellati a matrice di dimensioni', parameters.shape)
            else:
                raise ValueError(f"I dati nel file non sono compatibili con una matrice {n_layers}x{n_qubits}")
        except FileNotFoundError:
            print(f'File {param_filename} non trovato. Generazione di parametri casuali...')
        except ValueError as e:
            print(e)
            print("Errore con il file dei parametri. Generazione di parametri casuali...")

    if parameters is None:
        # Genera i parametri casualmente se non sono stati caricati correttamente
        parameters = 2 * np.pi * np.random.rand(n_layers, n_qubits)
        print('Parametri generati casualmente.')


    return parameters


def read_gates_var(output_dir, read_gates: str, n_qubits: int, n_layers: int):
    
    if read_gates == 'True':
        gates_var = np.loadtxt(f'{output_dir}/gates_var_n_{n_qubits}_qubits.txt', dtype=int)
        print(f'gates_var letto da {output_dir}')
        
    elif read_gates == 'False':
        gates_var = np.random.randint(0, 3, size=(n_layers, n_qubits))
        print('Rotazioni generate casualmente.')
    
    return gates_var




#==========================================================#
#
# Run Card
#
#==========================================================#

def print_run_card(output_dir : str, n_qubits : int, n_layers : int, parameters, val_g, n_shots, spop, ent_gate, lr, method: str, lambda1 = 0, shift=np.pi/2):

    # Controllo se output_dir è vuota
    if output_dir == '':
        output_dir = os.getcwd()  # Scegli la directory corrente
        print("Warning: output_dir is empty. Using current directory:", output_dir)


    # Controlla che la directory esista, altrimenti la crea
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    run_card_filename = f'run_card_{method}.txt'
    
    run_card_path = os.path.join(output_dir, run_card_filename)

    with open(run_card_path, "w") as f:
        
        f.write(" -------------------- \n")
        f.write("|                    |\n")
        f.write("|  RUNCARD           |\n")
        f.write("|                    |\n")
        f.write(" -------------------- \n")
        f.write("\n")
        f.write(f"Method           = {method}\n")
        f.write(f"Number of qubits = {n_qubits}\n")
        f.write(f"Number of Layers = {n_layers}\n")
        f.write(f"Number of Shots  = {n_shots}\n")
        f.write(f"Shift            = {shift}\n")
        f.write(f"Learning rate    = {lr}\n")

        f.write(f"Entangling gate  = {ent_gate}\n")
        if ent_gate == 0:
            f.write("Entaglment Gates are CNOTS \n")
        if ent_gate == 1:
            f.write("Entaglment Gates are SWAPS \n")

        # fixed lambda QNDM, if lambda1=0 -> DM 
        if lambda1 == 0:
            f.write("No Lambda coupling (DM run)\n")
        else:
            f.write("Lambda (QNDM coupling) = {}\n".format(lambda1))

        f.write("\n")
        f.write(f"Hamiltonian: {spop}\n\n")
        f.write(f"Parameters: {parameters}\n\n")
        f.write("Rotation U array: {} \n".format(val_g))
        
    print(f"runcard salvata in {run_card_path}")




#==========================================================#
#
# Printing section
#
#==========================================================#

def save_circuit_parameters(output_dir, parameters, gates_var):
    # Controlla che la directory esista, altrimenti la crea
    if output_dir == '':
        output_dir = os.getcwd()
    
    # Controlla che la directory esista, altrimenti la crea
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    parameters_filename = 'parameters.txt'
    gates_var_filename  = 'gates_var.txt'

    # Percorso completo del file
    parameters_path = os.path.join(output_dir, parameters_filename)
    gates_var_path  = os.path.join(output_dir, gates_var_filename)


    # Scrivi il gradiente nel file
    with open(parameters_path, 'w') as file:
        for par in parameters:
            file.write(f"{par}\n")

    # Scrivi il gradiente nel file
    with open(gates_var_path, 'w') as file:
        for gates in gates_var:
            file.write(f"{gates}\n")

    print(f"Parametri salvati in {parameters_path}")
    print(f"gates var salvati in {gates_var_path}")


def save_gradient(output_dir : str, gradient, method: str):
    # Controlla che la directory esista, altrimenti la crea
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Imposta il nome del file in base al metodo
    if method == 'psr':
        filename = 'gradient_psr.txt'
    elif method == 'zero_order':
        filename = 'gradient_zero_order.txt'
    else:
        raise ValueError("Metodo non riconosciuto: scegli 'psr' o 'zero_order'")
    
    # Percorso completo del file
    file_path = os.path.join(output_dir, filename)
    
    # Scrivi il gradiente nel file
    with open(file_path, 'w') as file:
        for grad in gradient:
            # Formatta ogni riga del gradiente con sei cifre decimali, separando i valori con spazi
            formatted_grad = " ".join(f"{g:.6f}" for g in grad)
            file.write(f"{formatted_grad}\n")

    print(f"Gradienti salvati in {file_path}")





#==========================================================#
#
# Optimization file
#
#==========================================================#

def print_opt_file(output_dir, iteration, cost, parameters, method: str, loop = True):

    if loop == True:
        filename = f"{output_dir}/output_{method}.dat"

        with open(filename, 'a') as f:  # 'a' per appendere i dati
            f.write(f"Iteration:    {iteration}\n")
            f.write(f"Cost function: {cost}    \n")
            f.write(f"Parameters:              \n")
            f.write(f"{parameters}")
            f.write(f"\n\n")

    elif loop == False:

        filename = f"{output_dir}/output_final_{method}.dat"

        with open(filename, 'a') as f:  # 'a' per appendere i dati
            f.write(f"Final cost value: {cost}    \n")
            f.write(f"Final parameters:           \n")
            f.write(f"{parameters}")
            f.write(f"\n\n")




#==========================================================#
#
# Optimization file
#
#==========================================================#

def print_file(output_dir, n_iterations, cost_list, time_list, method: str):


    # Controllo se output_dir è vuota
    if output_dir == '':
        output_dir = os.getcwd()  # Scegli la directory corrente
        print("Warning: output_dir is empty. Using current directory:", output_dir)


    # Crea la directory se non esiste
    os.makedirs(output_dir, exist_ok=True)

    filename = f"{output_dir}/output_{method}.txt"


    with open(filename, 'a') as f:  # 'a' to append data
        f.write("Iteration\tCost function value\tIteration time (s)\n")
        for i in range(n_iterations):
            f.write(f"{i + 1}\t{cost_list[i]}\t{time_list[i]:.4f}\n")





#==========================================================#
#
# Save and load a SparsePauliOp hamiltonian into a pickle
# file
#
#==========================================================#

def save_hamiltonian(sparse_op, filename):
    with open(filename, 'wb') as f:
        pickle.dump(sparse_op, f)

def load_hamiltonian(folder, filename):
    filepath = os.path.join(folder, filename)
    with open(filepath, 'rb') as f:
        return pickle.load(f)



#==========================================================#
#
# Recover the optimization data from the log file. The 
# output has the same shape of the output_dm.txt and 
# output_qndm.txt
#
#==========================================================#

import re

def recover_form_log(file_path, output_path):
    """
    Extract the optimization data from the log file and then it stores them into a .txt file.
    This file has the same structure of the final optimization file output_dm or output_qndm 
    
    Args:
        file_path (str): Percorso del file di log da analizzare.
        output_path (str): Percorso del file di output per salvare i dati estratti.

    Returns:
        None
    """
    # opening of the input log file
    with open(file_path, 'r') as infile:
        lines = infile.readlines()

    # Lista per memorizzare i risultati
    extracted_data = []

    # Regex to extract data from the row. All the rows we want to extraxct information are of the form
    # Iteration :  13%|█▎        | 126/1000 [4:49:25<33:38:11, 138.55s/it]127/1000, -73.97042260614565

    pattern = r"Iteration\s*:\s*\d+\%\|.*?(\d+)/\d+\s*\[(.*?s/it)\].*?,\s*(-?\d+\.\d+)"

    # add the header to the output file 
    extracted_data.append("Iteration\tCost function value\tIteration time (s)\n")

    # Cicla attraverso le righe e applica la regex
    for line in lines:
        match = re.search(pattern, line)
        if match:
            # Estrai il numero di iterazione, il tempo e il valore della funzione di costo
            iteration = match.group(1)
            time_spent = match.group(2)

            # Estrai solo la parte numerica del tempo (prima della virgola)
            time_spent = time_spent.split(',')[1].strip().replace('s/it', '')  # Rimuovi 's/it' e prendi solo la parte numerica
            cost_value = match.group(3)

            # Aggiungi i dati estratti alla lista in formato tabellare
            extracted_data.append(f"{iteration}\t{cost_value}\t{time_spent}\n")

    # Scrivi i dati estratti nel file di output
    with open(output_path, 'w') as outfile:
        outfile.writelines(extracted_data)
