#!/usr/bin/python3

#Code for gradient evaluation
#using the QNDM method.
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.

import os
import sys
import pandas as pd
import numpy as np

#QNDM error 
def get_qndm_error(pars, G_real_qndm, lambda1, shots, shift=np.pi/2):
    error = np.zeros(pars)
    
    for i in range(pars):
        error[i] = (1/(2*np.sqrt(shots)))*(1/(np.sqrt(1-(np.sin(lambda1*np.sin(shift)*G_real_qndm[i]))**2)))*(1/(lambda1*np.sin(shift)))
    
    return error


#DM error
def get_dm_error(pars, spop, shots, shift=np.pi/2):
    error = np.zeros(pars)
    error[:] = (1/(np.sqrt(2*shots)*np.sin(shift)))*sum(np.real(np.asarray(spop.coeffs)))
    
    return error

