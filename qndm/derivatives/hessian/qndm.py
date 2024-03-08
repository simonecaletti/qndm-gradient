#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_hessian import *
import numpy as np
from math import cos, sin

def qndm_hessian_circuit(circ,shift_position,sh2,pm,num_qub,num_l,q_d,shift,val_g,simp,ent_gate):


  #initialization of detector
  circ.h(q_d)
  params = ParameterVector("theta", length=len(val_g))

  qubits = []
  qubits_ = []
  for i in range(num_qub):
    qubits.append(i)
    qubits_.append(i)
  qubits_.insert(0,num_qub)

  lay_qn=U1_hess(val_g,params,num_qub,num_l,shift_position,sh2,shift,ent_gate)
 

  circ.compose(lay_qn, qubits=qubits, inplace=True)


  #first exponential interation
  
  evo_time = Parameter('p_deco')
  trotterized_op = PauliEvolutionGate(pm,evo_time)
  circ.append(trotterized_op, qubits_)                            

  
  #simplification part
  div = 0
  if simp == True:
    if shift_position <= sh2:
      div = int(shift_position/num_qub)
    if shift_position > sh2:
      div = int(sh2/num_qub) 
    #if we are above the first layer we have a part of U1+ that can be simplified with U2
    if div>=1:
      for i in range(div):
        val_g =  val_g[num_qub:]

  l_2=U2_hess(val_g,params,num_qub,num_l,shift_position,sh2,shift,ent_gate)
  circ.compose(l_2, qubits=qubits, inplace=True)

  #second exponential interation
  evo_time2 = Parameter('p_deco2')
  trotterized_op2 = PauliEvolutionGate(pm,-evo_time2)
  circ.append(trotterized_op2, qubits_)

  #Update U3
  l_3=U3_hess(val_g,params,num_qub,num_l,shift_position,sh2,shift,ent_gate)
  circ.compose(l_3, qubits=qubits, inplace=True)

  #third exponential interation
  evo_time3 = Parameter('p_deco3')
  trotterized_op3 = PauliEvolutionGate(pm,evo_time3)
  circ.append(trotterized_op3, qubits_)

  #Update U4
  l_4=U4_hess(val_g,params,num_qub,num_l,shift_position,sh2,shift,ent_gate)
  circ.compose(l_4, qubits=qubits, inplace=True)

  #fourth exponential interation
  evo_time4 = Parameter('p_deco4')
  trotterized_op4 = PauliEvolutionGate(pm,-evo_time4)
  circ.append(trotterized_op4, qubits_)
  
  #qubit for measure the detector
  circ.h(q_d)
  circ.s(q_d)
  circ.h(q_d)