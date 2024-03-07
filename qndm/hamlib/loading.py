#!/usr/bin/python3

# Contains snippets provided in Appendix A.2
# of paper arXiv:2306.13126

#-------------------------------------------------------------------------------
import networkx as nx
import mat2qubit as m2q
import openfermion as of
import h5py
import numpy as np
import zipfile
import requests
from io import BytesIO

#---------------------------------------------------------------------------------

def read_graph_hdf5( fname_hdf5 : str , key : str ):
    """ Read networkx graphs from HDF5 file at specified key . Returns a single networkx
    graph .
    """
    with h5py.File(fname_hdf5, "r") as f:
        G = nx.Graph(list(np.array(f[key])))
    
    return G

def read_gridpositions_hdf5(fname_hdf5 : str , key : str ):
    """ Read grid positions , stored as attribute of each networkx graph from HDF5 file
    at specified key . Returns grid positions of nodes associated with a single graph .
    """
    with h5py.File(fname_hdf5, "r") as f:
        dataset = f[key]
        gridpositions_dict = dict(dataset.attrs.items())
    
    return gridpositions_dict

def read_openfermion_hdf5 (fname_hdf5 : str , key : str , optype = of.QubitOperator):
    """ Read any openfermion operator object from HDF5 file at specified key .
    ’ optype ’ is the op class , can be of . QubitOperator or of . Fe rm io n Op er at o r .
    """
    with h5py.File( fname_hdf5 , "r" , libver = "latest") as f:
        op = optype(f[key][()].decode("utf-8"))
    
    return op

def read_mat2qubit_hdf5(fname_hdf5 : str ,key : str ):
    """ Returns mat2qubit ’s qSymbOp operator from HDF5 file at specified key . """
    with h5py.File(fname_hdf5, "r" ) as f:
        op = m2q.qSymbOp(f[key][()].decode("utf-8"))
    
    return op

def read_clause_list_hdf5(fname_hdf5 : str ,key : str):
    """ Read clause list from HDF5 file at specified key . Returns clause list in DIMACS
    format . """
    clause_list = []
    with h5py.File(fname_hdf5, "r") as f:
        for clause in list(np.array(f[key])):
            clause_list.append([v for v in clause])
    
    return clause_list


