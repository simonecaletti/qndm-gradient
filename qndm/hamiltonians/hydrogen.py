#!/usr/bin/python3

#--------------------------------------------------

from qndm.hamiltonians.interface import make_request, get_hamiltonian, get_nqubit
from qndm.hamlib.structure import get_hdf5_keys
import sys
from string import Template

#-------------------------------------------------

def get_model(n, shape, r, selected_key):

    init_url = "https://portal.nersc.gov/cfs/m888/dcamps/hamlib/chemistry/electronic/hydrogen_data/"

    #Cases explored by hand (the number is n)
    end_url2  = "H2_linear/ES_H2_linear_R${r}_sto-6g.hdf5.zip"
    end_url4  = "H4_linear/ES_H4_linear_R${r}_sto-6g_ham.hdf5.zip"
    end_url6  = "H6_linear/ES_H6_linear_R${r}_sto-6g_ham.hdf5.zip"
    end_url8  = "H8_linear/ES_H8_linear_R${r}_sto-6g_ham.hdf5.zip"
    end_url10 = "H10_${s}/ES_H10_${s}_R${r}_sto-6g_ham.hdf5.zip"
    end_url12 = "H12_${s}/ES_H12_${s}_R${r}_sto-6g_ham.hdf5.zip"
    end_url14 = "H14_${s}/ES_H14_${s}_R${r}_sto-6g_ham.hdf5.zip"
    end_url16 = "H16_${s}/ES_H16_${s}_R${r}_sto-6g_ham.hdf5.zip"
    end_url18 = "H18_linear/ES_H18_linear_R${r}_sto-6g.hdf5.zip"
    end_url20 = "H20_linear/ES_H20_linear_R${r}_sto-6g.hdf5.zip"
    end_url24 = "H24_linear/ES_H24_linear_R${r}_sto-6g.hdf5.zip"
    end_url28 = "H28_linear/ES_H28_linear_R${r}_sto-6g.hdf5.zip"
    end_url32 = "H32_linear/ES_H32_linear_R${r}_sto-6g.hdf5.zip"
    end_url36 = "H36_linear/ES_H36_linear_R${r}_sto-6g.hdf5.zip"
    end_url40 = "H40_linear/ES_H40_linear_R1.0_sto-6g.hdf5.zip"
    end_url50 = "H50_linear/ES_H50_linear_R1.0_sto-6g.hdf5.zip"
    end_url60 = "H60_linear/ES_H60_linear_R1.0_sto-6g.hdf5.zip"
    end_url70 = "H70_linear/ES_H70_linear_R1.0_sto-6g.hdf5.zip"
    end_url80 = "H80_linear/ES_H80_linear_R1.0_sto-6g.hdf5.zip"
    
    # Create all the possible url combination
    if n in [10, 12, 14, 16]: # Cases explored by hand
        if shape in ["linear", "pyramid", "ring", "sheet"]:
            if r in [x / 10 for x in range(5, 21)] :
                
                if n == 10: t = Template(end_url10)
                elif n == 12: t = Template(end_url12)
                elif n == 14: t = Template(end_url14)
                elif n == 16: t = Template(end_url16)
    
                end_url = t.safe_substitute(s=shape, r=r)
                #print(end_url)
                newurl = init_url + end_url
            
            else:
                print("Error: r value not valid.")
                sys.exit()

        else:
            print("Error: shape info not valid.")
            sys.exit()

    elif n in [2, 4, 6, 8, 18, 20, 24, 28, 32, 36, 40, 50, 60, 70, 80] and shape == "linear":
        if r in [x / 10 for x in range(5, 21)]:
                
            if n == 2: t = Template(end_url2)
            elif n == 4: t = Template(end_url4)
            elif n == 6: t = Template(end_url6)
            elif n == 8: t = Template(end_url8)
            elif n == 18: t = Template(end_url18)
            elif n == 20: t = Template(end_url20)
            elif n == 24: t = Template(end_url24)
            elif n == 28: t = Template(end_url28)
            elif n == 32: t = Template(end_url32)
            elif n == 36: t = Template(end_url36)
            elif n == 40: t = Template(end_url40)
            elif n == 50: t = Template(end_url50)
            elif n == 60: t = Template(end_url60)
            elif n == 70: t = Template(end_url70)
            elif n == 80: t = Template(end_url80)

            end_url = t.safe_substitute(r=r)
            #print(end_url)
            newurl = init_url + end_url

        else:
            print("Error: r value not valid.")
            sys.exit()

    else:
        print("Error: n value or shape not valid.")
        sys.exit()
        
    # Call the newurl
    # (we use the status_code property to check if it is valid)
    model = {}    
    print("Request sent to: {}".format(newurl))
    fname, opfname, status = make_request(newurl)
    
    if status == 200:
        model["url"] = newurl
        model["fname"] = fname
        model["opfname"] = opfname
        #print(model)
        
    
    elif status == 404:
        print("Error: the provided url is not valid.")
        sys.exit()
    
    # Add nqubit, key and hamiltonian info
    print("Adding nqubit and key info...", end="")
    keys = get_hdf5_keys(model["opfname"])
    if selected_key in keys:
        model["key"] = selected_key
        model["nqubit"] = get_nqubit(model["opfname"], selected_key)
        model["PS"], model["cps"] = get_hamiltonian(model["opfname"], model["key"], model["nqubit"])
    else:
        print("Error: selected key not available for the model.")
        sys.exit()
    print("done!")

    return model
    
