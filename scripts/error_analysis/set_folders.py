#!/usr/bin/env python3

#==========================================================#
#
# QISKIT code to generate optimization folders with
# circuits and parameters.
#==========================================================#

import os
import sys
import numpy as np
import shutil

# Carica le variabili d'ambiente dal file .env
load_dotenv()

# Legge la variabile e la aggiunge a sys.path
project_path = os.getenv("QNDM_PATH")
if project_path and project_path not in sys.path:
    sys.path.append(project_path)


# 
# QNDM packages importation
#----------------------------------------------------------#

from qndm.core import *

# 
# utils packages importation
#----------------------------------------------------------#

from utils.tools import *

#==========================================================#
#
# Utility functions
#==========================================================#

def get_base_directory():
    """Returns the base directory where the script is located."""
    return os.path.dirname(os.path.abspath(__file__))

def create_directory(path):
    """Creates a directory if it does not already exist."""
    os.makedirs(path, exist_ok=True)
    print(f"Directory '{path}' created or already exists.")

def save_to_file(data, filename, fmt):
    """Saves data to a file with the specified format."""
    np.savetxt(filename, data, fmt=fmt)
    print(f"File saved: {filename}")

#==========================================================#
# Parameter and gate generation
#==========================================================#

def generate_parameters(output_dir, n_qubits, n_layers):
    """Generates random circuit parameters and saves them."""
    parameters = 2 * np.pi * np.random.rand(n_layers, n_qubits)
    param_filename = os.path.join(output_dir, f'parameters_n_{n_qubits}_qubits.txt')
    save_to_file(parameters.flatten(), param_filename, fmt='%.18e')
    return parameters

def generate_gates_var(output_dir, n_qubits, n_layers):
    """Generates random gate variables and saves them."""
    gates_var = np.random.randint(1, 4, size=n_layers * n_qubits)
    gates_var_filename = os.path.join(output_dir, f'gates_var_n_{n_qubits}_qubits.txt')
    save_to_file(gates_var, gates_var_filename, fmt='%d')
    return gates_var


#==========================================================#
# Script creation
#==========================================================#

def create_script(output_dir, script_name, script_code, variables):
    """
    Creates a bash script for execution.
    """
    
    script_path = os.path.join(output_dir, script_name)
    with open(script_path, "w") as script_file:
        script_file.write("#!/bin/bash\n\n")
        script_file.write("# Deactivate any active virtual environment\n")
        script_file.write("deactivate 2>/dev/null || true\n\n")
        script_file.write("# Activate Conda environment\n")
        script_file.write("source ~/anaconda3/bin/activate\n")
        script_file.write("conda activate /home/dmelegari/anaconda3\n\n")
        script_file.write("# Debug information\n")
        script_file.write('echo "Python path: $(which python3)"\n')
        script_file.write('echo "Python version: $(python3 --version)"\n')
        script_file.write('echo "Pip path: $(which pip3)"\n')
        script_file.write("pip3 show qiskit\n\n")
        script_file.write("# Add project directory to PYTHONPATH\n")
        script_file.write("export PYTHONPATH=$PYTHONPATH:$(pwd)\n")
        script_file.write('echo "PYTHONPATH: $PYTHONPATH"\n\n')
        script_file.write(f"cd {output_dir} || {{ echo 'Directory not found'; exit 1; }}\n\n")
        script_file.write(f"echo 'Running {script_name} script...'\n")
        script_file.write(script_code.format(**variables))
    os.chmod(script_path, 0o755)
    
    print(f"Script created: {script_path}")


#==========================================================#
# Hamiltonian pre-generation
#==========================================================#

def pregenerate_hamiltonians(base_dir, n_qubits, pauli_strings_list):
    """
    Pre-generates Hamiltonians for each Pauli string count and saves them.
    """
    hamiltonians = {}
    for n_ps in pauli_strings_list:
        hamiltonian_dir = os.path.join(base_dir, f"hamiltonians/{n_ps}_pauli_strings")
        create_directory(hamiltonian_dir)

        # Generate DM Hamiltonian
        PS_dm, cps_dm = get_hamiltonian(n_qubits, n_ps, sel = str('ising_zero'), mu = mu, sigma = sigma)
        spop_dm = get_SparsePauliOp(PS_dm, cps_dm)
        dm_path = os.path.join(hamiltonian_dir, f"hamiltonian_{n_ps}_pauli_strings_dm.pkl")
        save_hamiltonian(spop_dm, dm_path)

        # Generate QNDM Hamiltonian
        PS_QNDM, cps_QNDM = add_detector(PS_dm, cps_dm)
        spop_qndm = get_SparsePauliOp(PS_QNDM, cps_QNDM)
        qndm_path = os.path.join(hamiltonian_dir, f"hamiltonian_{n_ps}_pauli_strings_qndm.pkl")
        save_hamiltonian(spop_qndm, qndm_path)

        # Calculate lambda value
        one_norm = np.sum(np.abs(cps_dm))
        lambda_value = 1 / np.sqrt(one_norm) # lambda1 #lambda1 = 1 / np.sqrt(one_norm) #1 / one_norm # lambda1 #1 / np.sqrt(one_norm)
        lambda_path = os.path.join(hamiltonian_dir, f"lambda_{n_ps}_pauli_strings.txt")
        with open(lambda_path, "w") as lambda_file:
            lambda_file.write(f"{lambda_value:.18e}\n")

        # Store paths and lambda
        hamiltonians[n_ps] = {
            "dm": dm_path,
            "qndm": qndm_path,
            "lambda": lambda_value
        }
        print(f"Hamiltonians for {n_ps} Pauli strings pre-generated and saved.")
    return hamiltonians

#==========================================================#
# Creates the optimization sets and shares pre-generated 
# Hamiltonians
#==========================================================#

def create_sets(base_dir, n_sets, n_qubits, n_layers, scripts, pauli_strings_list, hamiltonians):

    all_set_paths = []

    for i in range(n_sets):
        set_dir = os.path.join(base_dir, f'set_{i + 1}')
        create_directory(set_dir)

        for n_ps in pauli_strings_list:
            pauli_dir = os.path.join(set_dir, f"{n_ps}_pauli_strings")
            create_directory(pauli_dir)
            all_set_paths.append(pauli_dir)

            # Copy pre-generated Hamiltonians
            pregen_hamiltonian = hamiltonians[n_ps]
            os.symlink(pregen_hamiltonian["dm"], os.path.join(pauli_dir, os.path.basename(pregen_hamiltonian["dm"])))
            os.symlink(pregen_hamiltonian["qndm"], os.path.join(pauli_dir, os.path.basename(pregen_hamiltonian["qndm"])))

            # Save lambda value
            lambda_path = os.path.join(pauli_dir, f"lambda_{n_ps}_pauli_strings.txt")
            with open(lambda_path, "w") as lambda_file:
                lambda_file.write(f"{pregen_hamiltonian['lambda']:.18e}\n")

            # Generate and save parameters and gates
            generate_parameters(pauli_dir, n_qubits, n_layers)
            generate_gates_var(pauli_dir, n_qubits, n_layers)

            # Update scripts configuration for this set
            for script_name, script_details in scripts.items():
                variables = script_details["variables"].copy()
                variables["output_dir"] = pauli_dir
                variables["n_pauli_strings"] = n_ps
                variables["lambda1"] = pregen_hamiltonian["lambda"]

                create_script(pauli_dir, script_name, script_details["code"], variables)

            print(f"Subdirectory with {n_ps} Pauli strings created in '{pauli_dir}'.")

        print(f"Set {i + 1} created in '{set_dir}'.")

    # Write all paths to .folder_0 file without overwriting and add a trailing comma
    folder_0_path = os.path.join(get_base_directory(), ".folder_0")
    with open(folder_0_path, "a") as folder_0_file: 
        folder_0_file.write(", ".join(all_set_paths) + ", ")
    print(f"Paths appended to file: {folder_0_path}")


if __name__ == "__main__":
    
    n_sets             = 10
    n_qubits           = 10
    n_layers           = 15
    lay_u              = 1
    ent_gate           = 0
    shift              = np.pi / 2
    noise              = False
    lr                 = 0.05
    n_gradients        = 25
    pauli_strings_list = [10, 20, 50, 100, 250, 500, 750, 1000]
    read_pars          = "True"
    read_gates         = "True"
    lambda1            = 1
    mu                 = 0
    sigma              = 5

    shots_list = [500]
    
    for shots in shots_list:
        # Define paths
        base_dir = os.path.join(get_base_directory(), f"500_shots_lambda_sqrt_norm/results_dm_and_qndm_{shots}_shots")
        create_directory(base_dir)
        
        # Pre-generate Hamiltonians
        hamiltonians = pregenerate_hamiltonians(base_dir, n_qubits, pauli_strings_list)
        
        # Scripts configuration
        scripts = {
            "run.sh": {
                "code": "python3 {dm_code} {n_qubits} {n_layers} {shots} {lay_u} {ent_gate} {shift} {noise} {lambda1} {n_gradients} {n_pauli_strings} {read_pars} {read_gates} {output_dir}\n",
                "variables": {
                    "dm_code": os.path.join(get_base_directory(), "main.py"),
                    "n_qubits": n_qubits,
                    "n_layers": n_layers,
                    "shots": shots,
                    "lay_u": lay_u,
                    "ent_gate": ent_gate,
                    "shift": shift,
                    "noise": noise,
                    "lambda1": lambda1,
                    "n_gradients": n_gradients,
                    "n_pauli_strings": "{n_pauli_strings}",
                    "read_pars": read_pars,
                    "read_gates": read_gates,
                    "output_dir": "{output_dir}"
                }
            }
        }
        
        # Generate sets 
        create_sets(base_dir, n_sets, n_qubits, n_layers, scripts, pauli_strings_list, hamiltonians)
    
    print("All folders created.")






