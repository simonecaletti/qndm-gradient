#!/usr/bin/python3

#-------------------------------------------------------------

import zipfile
import requests
from io import BytesIO
import openfermion as of
import h5py
import sys

#-------------------------------------------------------------

# url = "url_to_hamlib"
# filename, opfname = make_request(url)
# key = "/ham_BK/"
# num_qub = get_nqubit(opfname, key)
# PS, cps = get_hamiltonian(opfname, key, num_qub)
# spop = get_SparsePauliOp(PS, cps)

def make_request(url):

    r = requests.get(url, stream = True) 

    if r.status_code == 200:
        z = zipfile.ZipFile(BytesIO(r.content))
        filename = z.namelist()[0]
        opfname = z.open(filename, "r")

    elif r.status_code == 404:
        filename = ""
        opfname = ""

    return filename, opfname, r.status_code

def get_nqubit(op_fname, key):
    
    with h5py.File(op_fname, "r") as f:
        ham_str = f[key][()].decode("utf-8") #This object is a string
 
    protolist = ham_str.split("+\n")
    for pl in protolist:
        x = pl.split(" ", 1)
        new_str = x[1].strip()
        new_str = new_str[1:-1]
        PS_list = new_str.split(" ")

        max = 0
        for gate in PS_list:
            if gate == "": 
                continue
            else:
                q_index = int(gate[-1])
                if q_index > max: max = q_index

    nqubit = max + 1

    return nqubit

def get_hamiltonian(op_fname, key, num_qub):

    with h5py.File(op_fname, "r") as f:
        ham_str = f[key][()].decode("utf-8") #This object is a string
    
    PS, cps = str_to_ham(ham_str, num_qub)

    return PS, cps

def combine_PS(PS_list, num_qub):
    
    PS = ""
    insert_identity = True

    if PS_list[0] == "" and len(PS_list) == 1:
        for q in range(num_qub):
            PS += "I"
    
    else:
        for q in range(num_qub):

            for el in PS_list:
                if int(el[-1]) == q:
                    PS += el[0]
                    insert_identity = False

            if insert_identity:
                PS += "I" 

            insert_identity = True
    
    return PS

def str_to_PS(PS_str, num_qub):

    new_str = PS_str.strip()
    new_str = new_str[1:-1]
    PS_list = new_str.split(" ")
    PS = combine_PS(PS_list, num_qub)

    return PS

def str_to_coeff(cps_str):
    
    if cps_str[0] == "(":
        new_str = cps_str[1:-1]
        new_cps = new_str.split("+")
        real = new_cps[0]
        img = new_cps[1]
    
        if img != "0j":
            print("Error: cps has an imaginary part.")
            sys.exit()
        else:
            cps = float(real)

    else:
        cps = float(cps_str)

    return cps

def str_to_ham(ham_str, num_qub):
    
    protolist = ham_str.split("+\n")

    PS = []
    cps = []
    for pl in protolist:
        x = pl.split(" ", 1)
        xcps = str_to_coeff(x[0])
        xPS = str_to_PS(x[1], num_qub)
        
        cps.append(xcps)
        PS.append(xPS)

    return PS, cps
