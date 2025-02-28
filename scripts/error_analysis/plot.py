import os
import numpy as np
import matplotlib.pyplot as plt

def find_valid_directories(base_dir):
    valid_dirs = []
    for root, dirs, _ in os.walk(base_dir):
        if "hamiltonians" in root:
            continue  # Salta la cartella hamiltonians
        
        set_dirs = [d for d in dirs if d.startswith("set_")]
        for set_dir in set_dirs:
            valid_dirs.append(os.path.join(root, set_dir))
    
    return valid_dirs

def load_data(file_path):
    data = np.loadtxt(file_path, skiprows=1)  # Salta l'intestazione
    return data[:, 0], data[:, 1]  # N_ps, valori corrispondenti

def plot_data(N_ps, values1, values2, label1, label2, title, filename):
    plt.figure()
    plt.plot(N_ps, values1, marker='o', linestyle='-', label=label1)
    plt.plot(N_ps, values2, marker='s', linestyle='-', label=label2)
    plt.xlabel('J')
    plt.ylabel(title)
    plt.legend()
    plt.grid()
    plt.title(title)
    plt.savefig(filename)
    plt.close()

def plot_with_errorbars(N_ps, mean1, mean2, err1, err2, label1, label2, title, filename):
    # Filtra i dati per rimuovere i punti con N_ps = 10 e 20
    mask = ~np.isin(N_ps, [10, 20])
    
    N_ps_filtered = np.array(N_ps)[mask]
    mean1_filtered = np.array(mean1)[mask]
    mean2_filtered = np.array(mean2)[mask]
    err1_filtered = np.array(err1)[mask]
    err2_filtered = np.array(err2)[mask]

    plt.figure()
    plt.errorbar(N_ps_filtered, mean1_filtered, yerr=err1_filtered, fmt='o', color='red', label=label1, capsize=5)
    plt.errorbar(N_ps_filtered, mean2_filtered, yerr=err2_filtered, fmt='s', color='blue', label=label2, capsize=5)
    plt.xlabel(r'$J$')
    plt.ylabel(r'$\mu(g_{l}^{i}(J))  \pm    \mu(MSE^{i}(J))$')
    plt.xticks([50, 100, 250, 500, 750, 1000])  # Imposta i ticks personalizzati
    plt.legend()
    plt.grid(False)
    # plt.title(title)
    plt.savefig(filename)
    plt.close()




def process_and_plot(base_directory):
    directories = find_valid_directories(base_directory)
    
    for directory in directories:
        files = {
            "var_dm": os.path.join(directory, "var_dm.txt"),
            "var_qndm": os.path.join(directory, "var_qndm.txt"),
            "bias_dm": os.path.join(directory, "bias_dm.txt"),
            "bias_qndm": os.path.join(directory, "bias_qndm.txt"),
            "MSE_dm": os.path.join(directory, "MSE_dm.txt"),
            "MSE_qndm": os.path.join(directory, "MSE_qndm.txt"),
            "oom_dm": os.path.join(directory, "oom_dm.txt"),
            "oom_qndm": os.path.join(directory, "oom_qndm.txt"),
            "oom_sv": os.path.join(directory, "oom_sv.txt"),
            "mean_dm": os.path.join(directory, "mean_dm.txt"),
            "mean_qndm": os.path.join(directory, "mean_qndm.txt"),
        }
        
        # Caricamento dati
        N_ps, var_dm = load_data(files["var_dm"])
        _, var_qndm  = load_data(files["var_qndm"])
        _, bias_dm   = load_data(files["bias_dm"])
        _, bias_qndm = load_data(files["bias_qndm"])
        _, MSE_dm    = load_data(files["MSE_dm"])
        _, MSE_qndm  = load_data(files["MSE_qndm"])
        _, oom_dm    = load_data(files["oom_dm"])
        _, oom_qndm  = load_data(files["oom_qndm"])
        _, oom_sv    = load_data(files["oom_sv"])
        _, mean_dm   = load_data(files["mean_dm"])
        _, mean_qndm = load_data(files["mean_qndm"])
        
        # Creazione dei plot
        plot_data(N_ps, var_dm,  var_qndm,  "var_dm",  "var_qndm",  "Variance", os.path.join(directory, "var.png"))
        plot_data(N_ps, bias_dm, bias_qndm, "bias_dm", "bias_qndm", "Bias",     os.path.join(directory, "bias.png"))
        plot_data(N_ps, MSE_dm,  MSE_qndm,  "MSE_dm",  "MSE_qndm",  "MSE",      os.path.join(directory, "MSE.png"))
        
        # MSE normalizzato
        MSE_norm_dm   = np.sqrt(MSE_dm)   / np.abs(oom_sv)
        MSE_norm_qndm = np.sqrt(MSE_qndm) / np.abs(oom_sv)
        plot_data(N_ps, MSE_norm_dm, MSE_norm_qndm, "MSE_dm / OoM_dm", "MSE_qndm / OoM_qndm", "MSE Normalized", os.path.join(directory, "MSE_normalized.png"))
        
        # Plot con errorbar (media ± MSE)
        plot_with_errorbars(N_ps, np.abs(mean_dm), np.abs(mean_qndm), MSE_dm, MSE_qndm, "i=DM", "i=QNDM", "", os.path.join(directory, "mean_with_error.png"))

if __name__ == "__main__":
    base_directory = "/home/mele/cacchio/qndm/test_scripts/check_error/ising_mu0_sigma5"
    process_and_plot(base_directory)
