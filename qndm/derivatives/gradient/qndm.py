#!/usr/bin/python3
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_gradient import *


#Quantum Circuit
def qndm_gradient_circuit(circ,shift_position,pm,num_qub,num_l,val_g,q_d,shift, simp,ent_gate):

  #initialization of detector
  circ.h(q_d)

  #initialization of paramenters \theta vector 
  params = ParameterVector("theta", length=len(val_g))

  #Lists with the qubits position information to compose circ with the U and the PauliEvolutionGates
  qubits_U = []
  qubits_exp = []
  for i in range(num_qub):
    qubits_U.append(i)
    qubits_exp.append(i)
  qubits_exp.insert(0,num_qub)

  #first unitary trasformation: U1:|00...0>->|\psi(\theta - shift*e_(shift_position))
  unitary1=U1(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(unitary1, qubits=qubits_U, inplace=True)

  #first coupling interation
  evo_time = Parameter('p_deco')
  trotterized_op = PauliEvolutionGate(pm,evo_time)
  circ.append(trotterized_op, qubits_exp)

 
  #second unitary trasformation: U1_dag|\psi(\theta - shift*e_(shift_position))->|00...0>
  unitary1_dag=U1_dag(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(unitary1_dag, qubits=qubits_U, inplace=True)

  #third unitary trasformation: U2:|00...0>->|\psi(\theta + shift*e_(shift_position))
  unitary2=U2(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(unitary2, qubits=qubits_U, inplace=True)


  #second coupling interation
  evo_time2 = Parameter('p_deco2')
  trotterized_op2 = PauliEvolutionGate(pm,-evo_time2)
  circ.append(trotterized_op2, qubits_exp)
  

  #qubit for measure the detector
  circ.h(q_d)
  circ.s(q_d)
  circ.h(q_d)

  return None

