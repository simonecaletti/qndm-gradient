#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
from qndm.layers.unitaries_gradient import *
import numpy as np
from math import cos, sin

#Quantum Circuit
def qndm_gradient_circuit(circ,shift_position,pm,num_qub,num_l,val_g,q_d,shift, simp,ent_gate):

  #initialization of detector
  circ.h(q_d)
  params = ParameterVector("theta", length=len(val_g))

  qubits = []
  qubits_ = []
  for i in range(num_qub):
    qubits.append(i)
    qubits_.append(i)
  qubits_.insert(0,num_qub)

  lay_qn=U1(val_g,params,num_qub,num_l,shift,shift_position,ent_gate)
  circ.compose(lay_qn, qubits=qubits, inplace=True)

  #first exponential interation
  
  evo_time = Parameter('p_deco')
  trotterized_op = PauliEvolutionGate(pm,evo_time)
  circ.append(trotterized_op, qubits_)

  #U2 

  #Update U2
  val_g2 = val_g
  div = 0
  if simp == True:
    div = int(shift_position/num_qub)
    #if we are above the first layer we have a part of U1+ that can be simplified with U2
        if div>=1:
      for i in range(div):
        val_g2 =  val_g2[num_qub:]


  l_d=U1_dag(val_g2,params,num_qub,num_l,shift,shift_position,div,ent_gate)
  circ.compose(l_d, qubits=qubits, inplace=True)
  l_2=U2(val_g2,params,num_qub,num_l,shift,shift_position,div,ent_gate)
  circ.compose(l_2, qubits=qubits, inplace=True)


  #second exponential interation
  evo_time2 = Parameter('p_deco2')
  trotterized_op2 = PauliEvolutionGate(pm,-evo_time2)
  circ.append(trotterized_op2, qubits_)
  

  #qubit for measure the detector
  circ.h(q_d)
  circ.s(q_d)
  circ.h(q_d)

  return None

