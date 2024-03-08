#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_hessian import *
import numpy as np
from math import cos, sin

def dm_hessian_circuit(circ,sh1,sh2,num_qub,num_l,val_g,shift1,shift2,kk,ent_gate):
  #initialization of detector
  params = ParameterVector("theta", length=len(val_g))

  qubits = []
  for i in range(num_qub):
    qubits.append(i)

  l_d = U1_hess(val_g,params,num_qub,num_l,shift1,shift2,sh1,sh2,ent_gate)
  
  circ.compose(l_d, qubits=qubits, inplace=True)
  
  #observable
  gate_dm = [gate_dm for gate_dm in kk]
  for i,val_h in enumerate(gate_dm):
    if val_h == 'X':
      circ.h(num_qub-1-i)
      
    elif val_h == 'Y':
      circ.sdg(num_qub-1-i)
      circ.h(num_qub-1-i)

    elif val_h == 'Z':
      pass
    elif val_h == 'I':
      pass
    else:
      raise ValueError("You are not reading well spop!")
      

  return None
