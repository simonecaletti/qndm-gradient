#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_gradient import *
import numpy as np
from math import cos, sin

#Quantum Circuit
def dm_gradient_circuit(circ,shift_position,num_qub,num_l,val_g,shift,kk,ent_gate):
  #initialization of detector
  params = ParameterVector("theta", length=len(val_g))

  qubits = []
  for i in range(num_qub):
    qubits.append(i)


  l_d=U1(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(l_d, qubits=qubits, inplace=True)
  
  #observable
  gate_dm = [gate_dm for gate_dm in kk]
  print(gate_dm)
  for i,val_h in enumerate(gate_dm):

    if str(val_h) == 'X':
      circ.h(i)
      
    elif str(val_h) == 'Y':
      circ.sdg(i)
      circ.h(i)

    elif str(val_h) == 'Z':
      pass
    elif str(val_h) == 'I':
      pass
    else:
      raise ValueError("you are not reading well spop")
      

  return None
