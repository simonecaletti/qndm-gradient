#!/usr/bin/python3

from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_hessian import *

def qndm_hessian_circuit(circ,sh1,sh2,pm,num_qub,num_l,q_d,shift,val_g,ent_gate):


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


  #first unitary trasformation: U1:|00...0>->|\psi(\theta - shift*e_(sh1)+shift*e_(sh2))
  unitary1=U1_hess(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary1, qubits=qubits_U, inplace=True)


  #first coupling interation
  evo_time = Parameter('p_deco')
  trotterized_op = PauliEvolutionGate(pm,evo_time)
  circ.append(trotterized_op, qubits_exp)                            

  #second unitary trasformation: U1_dag|\psi(\theta - shift*e_(sh1)+shift*e_(sh2))->|00...0>
  unitary1_dag=U1_hess_dag(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary1_dag, qubits=qubits_U, inplace=True)

  #third unitary trasformation: U2:|00...0>->|\psi(\theta - shift*e_(sh1)-shift*e_(sh2))
  unitary2=U2_hess(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary2, qubits=qubits_U, inplace=True)

  #second coupling interation
  evo_time2 = Parameter('p_deco2')
  trotterized_op2 = PauliEvolutionGate(pm,-evo_time2)
  circ.append(trotterized_op2, qubits_exp)

  #fourth unitary trasformation: U2_dag|\psi(\theta - shift*e_(sh1)-shift*e_(sh2))->|00...0>
  unitary2_dag=U2_hess_dag(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary2_dag, qubits=qubits_U, inplace=True)

  #fifth unitary trasformation: U3:|00...0>->|\psi(\theta + shift*e_(sh1)-shift*e_(sh2))
  unitary3=U3_hess(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary3, qubits=qubits_U, inplace=True)

  #third coupling interation
  evo_time3 = Parameter('p_deco3')
  trotterized_op3 = PauliEvolutionGate(pm,evo_time3)
  circ.append(trotterized_op3, qubits_exp)

  #sixth unitary trasformation: U3_dag|\psi(\theta + shift*e_(sh1)-shift*e_(sh2))->|00...0>
  unitary3_dag=U3_hess_dag(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary3_dag, qubits=qubits_U, inplace=True)

  #seventh unitary trasformation: U3:|00...0>->|\psi(\theta + shift*e_(sh1)+shift*e_(sh2))
  unitary4=U4_hess(val_g,params,num_qub,num_l,sh1,sh2,shift,ent_gate)
  circ.compose(unitary4, qubits=qubits_U, inplace=True)

  #fourth coupling interation
  evo_time4 = Parameter('p_deco4')
  trotterized_op4 = PauliEvolutionGate(pm,-evo_time4)
  circ.append(trotterized_op4, qubits_exp)
  
  #qubit for measure the detector
  circ.h(q_d)
  circ.s(q_d)
  circ.h(q_d)