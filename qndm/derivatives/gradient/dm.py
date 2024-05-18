#!/usr/bin/python3

from qiskit.circuit import ParameterVector
from qndm.layers.unitaries_gradient import *


#Quantum Circuit
def dm_gradient_circuit(circ,shift_position,num_qub,num_l,val_g,shift,kk,ent_gate):

  #initialization of paramenters \theta vector 
  params = ParameterVector("theta", length=len(val_g))

  #List with the qubits position information to compose circ with the U 
  qubits_U = []
  for i in range(num_qub):
    qubits_U.append(i)

  #Unitary trasformation: U1:|00...0>->|\psi(\theta +- shift*e_(shift_position))
  unitary1=U1(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(unitary1, qubits=qubits_U, inplace=True)
  
  #Rotation of the state in function of Pauli String
  gate_dm = [gate_dm for gate_dm in kk]
  
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
