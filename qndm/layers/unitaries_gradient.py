#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
import numpy as np
from math import cos, sin

#-------------------------------------------------------------------------------

#####################################
#                                   #
#            Gradient               #      
#                                   #
#####################################

#U1
def U1(val_g,par_var,num_qu,num_l,shift,shift_pos,ent_gate):
  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")

#rotation part
  for k in range(val_g.shape[0]): 
    if val_g[k]==1:
      if k == shift_pos:
        circ.rx(par_var[k]-shift,k%num_qu)
      else:  
        circ.rx(par_var[k],k%num_qu)
    elif val_g[k]==2:
      if k == shift_pos:
        circ.ry(par_var[k]-shift,k%num_qu)
      else: 
        circ.ry(par_var[k],k%num_qu)
    elif val_g[k]==3:
      if k == shift_pos:
        circ.rz(par_var[k]-shift,k%num_qu)
      else: 
        circ.rz(par_var[k],k%num_qu)
   
#entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
          #circ.cx(ii*2, 2*ii+1)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            #circ.cx(2*j+1, 2*j+2)
            ent[ent_gate](2*j+1)
  return circ

#U1 dagger
def U1_dag(val_g,par_var,num_qu,num_l,shift,shift_p,div,ent_gate):
  
  val_to_p = 0
  if div >= 1:
    shift_pos = shift_p - num_qu*div
    val_to_p = num_qu * div

  else:
    shift_pos  = shift_p
  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")
  
  for k in reversed(range(val_g.shape[0])): 
#entanglement part
    if  k !=0 and (k+1)%((len(val_g)/(num_l))) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
            #circ.cx(2*j+1, 2*j+2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
          #circ.cx(ii*2, 2*ii+1)
        
#rotation part     
    if val_g[k]==1:
      if k == shift_pos:
        circ.rx(-par_var[k+val_to_p]+shift,k%num_qu)
      else:  
        circ.rx(-par_var[k+val_to_p],k%num_qu)
    elif val_g[k]==2:
      if k == shift_pos:
        circ.ry(-par_var[k+val_to_p]+shift,k%num_qu)
      else: 
        circ.ry(-par_var[k+val_to_p],k%num_qu)
    elif val_g[k]==3:
      if k == shift_pos:
        circ.rz(-par_var[k+val_to_p]+shift,k%num_qu)
      else: 
        circ.rz(-par_var[k+val_to_p],k%num_qu)

  return circ

#U2
def U2(val_g,par_var,num_qu,num_l,shift,shift_p,div,ent_gate):
  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")
  val_to_p = 0
  if div >= 1:
    shift_pos = shift_p - num_qu*div
    val_to_p = num_qu * div
  else:
    shift_pos = shift_p
    
  for k in range(val_g.shape[0]): 
#rotation part
    if val_g[k]==1:
      if k == shift_pos:
        circ.rx(par_var[k+val_to_p]+shift,k%num_qu)
      else:  
        circ.rx(par_var[k+val_to_p],k%num_qu)
    elif val_g[k]==2:
      if k == shift_pos:
        circ.ry(par_var[k+val_to_p]+shift,k%num_qu)
      else: 
        circ.ry(par_var[k+val_to_p],k%num_qu)
    elif val_g[k]==3:
      if k == shift_pos:
        circ.rz(par_var[k+val_to_p]+shift,k%num_qu)
      else: 
        circ.rz(par_var[k+val_to_p],k%num_qu)

#entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
          #circ.cx(ii*2, 2*ii+1)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            #circ.cx(2*j+1, 2*j+2) 
            ent[ent_gate](2*j+1)

  return circ
  


