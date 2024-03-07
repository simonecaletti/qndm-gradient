#!/usr/bin/python3

# Contains snippets provided in Appendix A.1
# of paper arXiv:2306.13126
#
##-------------------------------------------------------------------------------
import networkx as nx
import mat2qubit as m2q
import openfermion as of
import h5py
import numpy as np

#---------------------------------------------------------------------------------

def parse_through_hdf5(func):
    """ Decorator function that iterates through an HDF5 file and performs
    the action specified by ‘ func ‘ on the internal and leaf nodes in the HDF5 file . """
    def wrapper(obj ,path = "/" ,key = None):
        if type(obj) in [h5py._hl.group.Group, h5py._hl.files.File]:
            for ky in obj.keys():
                func(obj ,path ,key=ky ,leaf = False)
                wrapper(obj = obj[ky], path = path+ky+"/", key = ky)
        elif type(obj) == h5py._hl.dataset.Dataset:
            func(obj ,path ,key=None, leaf = True)
    return wrapper

def print_hdf5_structure(fname_hdf5 : str):
    """ Print the path structure of the HDF5 file .
    Args
    ----
    fname_hdf5 ( str ) : full path where HDF5 file is stored
    """
    @parse_through_hdf5
    def action(obj, path = "/", key = None, leaf = False):
        if key is not None :
            print((path.count ( "/") -1) * "\t" , "-" ,key , ":" ,path + key + "/")
        if leaf :
            print ((path.count( "/") -1) * "\t" , "[^^DATASET^^]")
    with h5py.File(fname_hdf5, "r") as f:
        action(f["/"])


def get_hdf5_keys (fname_hdf5 : str):
    """ Get a list of keys to all datasets stored in the HDF5 file .
    Args
    ----
    fname_hdf5 ( str ) : full path where HDF5 file is stored
    """
    all_keys = []
    @parse_through_hdf5
    def action (obj ,path = "/" , key = None , leaf = False):
        if leaf is True:
            all_keys.append(path)
            
    with h5py.File(fname_hdf5 ,"r") as f:
        action(f["/"])
    
    return all_keys

