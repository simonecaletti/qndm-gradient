#!/usr/bin/python3

#--------------------------------------------------

from qndm.hamiltonians.interface import make_request, get_hamiltonian, get_nqubit
from qndm.hamlib.structure import get_hdf5_keys
import sys
from string import Template 

#-------------------------------------------------

def get_model_lithium(selected_key):

    new_url = "https://portal.nersc.gov/cfs/m888/dcamps/hamlib/chemistry/electronic/standard/Li2.hdf5.zip"

    
    # Call the newurl
    # (we use the status_code property to check if it is valid)
    model = {}    
    print("Request sent to: {}".format(new_url))
    fname, opfname, status = make_request(new_url)
    
    if status == 200:
        model["url"] = new_url
        model["fname"] = fname
        model["opfname"] = opfname
        #print(model)
        
    
    elif status == 404:
        print("Error: the provided url is not valid.")
        sys.exit()
    
    # Add nqubit, key and hamiltonian info
    print("Adding nqubit and key info...", end="")
    keys = get_hdf5_keys(model["opfname"])

    keys = get_hdf5_keys(model["opfname"])
    if selected_key in keys:
        model["key"] = selected_key
        model["nqubit"] = get_nqubit(model["opfname"], selected_key)
        model["PS"], model["cps"] = get_hamiltonian(model["opfname"], model["key"], model["nqubit"])
    else:
        print("Error: selected key not available for the model.")
        sys.exit()
    print("done!")
   
    print("done!")

    return model
    
