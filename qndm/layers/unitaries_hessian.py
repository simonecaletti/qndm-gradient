#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Parameter,ParameterVector
from qiskit.circuit.library import PauliEvolutionGate
import numpy as np
from math import cos, sin

####################################
#                                   #
#           Hessian QNDM            #      
#                                   #
#####################################

#U1 hessian 
def U1_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):
    q_reg = QuantumRegister(num_qu, "q")
    circ = QuantumCircuit(q_reg, name="layer-")

    #rotation part
    for k in range(val_g.shape[0]): 
        sig = 1
        if val_g[k]==1:
            if sh2 == shift_pos and k == sh2:
                circ.rx(par_var[k],k%num_qu)
            elif k == shift_pos:
                circ.rx(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rx(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.rx(par_var[k],k%num_qu)

        elif val_g[k]==2:
            if sh2 == shift_pos and k == sh2:
                circ.ry(par_var[k],k%num_qu)
            elif k == shift_pos:
                circ.ry(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                circ.ry(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.ry(par_var[k],k%num_qu)

        elif val_g[k]==3:
            if sh2 == shift_pos and k == sh2:
                circ.rz(par_var[k],k%num_qu)
            elif k == shift_pos:
                circ.rz(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rz(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.rz(par_var[k],k%num_qu)

        #entanglement part
        if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
            ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
            if num_qu!=1:
                r_cx = int(num_qu/2)
                for ii in range(r_cx):
                    ent[ent_gate](2*ii)
                if num_qu>2:
                    for j in range(r_cx-(num_qu+1)%2):
                        ent[ent_gate](2*j+1)

    return circ



#U1 hessian dagger
def U1_hess_dag(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):
  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")
  val_to_p = 0
  if div >= 1:
    
    sh2 = sh2 - num_qu*div
    shift_pos = shift_pos - num_qu*div
    val_to_p = num_qu * div
  
  
  for k in reversed(range(val_g.shape[0])): 
    sig = 1

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
        
    #rotation part     
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
        circ.rx(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
        circ.ry(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
        circ.rz(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k+val_to_p],k%num_qu)

  return circ 

#U2 hessian
def U2_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):
    q_reg = QuantumRegister(num_qu, "q")
    circ = QuantumCircuit(q_reg, name="layer-")
    val_to_p = 0
    if div >= 1:
    
        sh2 = sh2 - num_qu*div
        shift_pos = shift_pos - num_qu*div
        val_to_p = num_qu * div

    for k in range(val_g.shape[0]): 
    #rotation part
        sig = -1
        if val_g[k]==1:
            if sh2 == shift_pos and k == sh2:
                circ.rx(par_var[k+val_to_p]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
            else:  
                circ.rx(par_var[k+val_to_p],k%num_qu)


        elif val_g[k]==2:
            if sh2 == shift_pos and k == sh2:
                circ.ry(par_var[k+val_to_p]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
            else:  
                circ.ry(par_var[k+val_to_p],k%num_qu)

        elif val_g[k]==3:
            if sh2 == shift_pos and k == sh2:
                circ.rz(par_var[k+val_to_p]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
            else:  
                circ.rz(par_var[k+val_to_p],k%num_qu)

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0:
        ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
        if num_qu!=1:
            r_cx = int(num_qu/2)
            for ii in range(r_cx):
                ent[ent_gate](2*ii)
            if num_qu>2:
                for j in range(r_cx-(num_qu+1)%2):
                    ent[ent_gate](2*j+1)

         
    return circ


#U2 hessian dag 
def U2_hess_dag(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):

  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")

  val_to_p = 0
  if div >= 1:
    
    sh2 = sh2 - num_qu*div
    shift_pos = shift_pos - num_qu*div
    val_to_p = num_qu * div
  
  for k in reversed(range(val_g.shape[0])): 
    sig = -1

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
          
        
    #rotation part     
    if val_g[k]==1:

         
      if sh2 == shift_pos and k == sh2:
          circ.rx(-par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
          circ.ry(-par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
          circ.rz(-par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k+val_to_p],k%num_qu)

  return circ 


#U3 hessian 
def U3_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):

  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")

  val_to_p = 0
  if div >= 1:

    sh2 = sh2 - num_qu*div
    shift_pos = shift_pos - num_qu*div
    val_to_p = num_qu * div
 
  for k in range(val_g.shape[0]): 
    #rotation part
    sig = -1
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
          circ.rx(par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.rx(par_var[k+val_to_p],k%num_qu)


    elif val_g[k]==2:
      if sh2 == shift_pos and k == sh2:
          circ.ry(par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = +1
        circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.ry(par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==3:
      if sh2 == shift_pos and k == sh2:
          circ.rz(par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.rz(par_var[k+val_to_p],k%num_qu)

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)

  return circ


#U3 hessian dag 
def U4_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):

  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")

  val_to_p = 0
  if div >= 1:
   
      sh2 = sh2 - num_qu*div
      shift_pos = shift_pos - num_qu*div
      val_to_p = num_qu * div
  
  for k in reversed(range(val_g.shape[0])): 
    sig = -1

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)

        
    #rotation part     
    if val_g[k]==1:

         
      if sh2 == shift_pos and k == sh2:
          circ.rx(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
          circ.ry(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.ry(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
          circ.rz(-par_var[k+val_to_p],k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(-par_var[k+val_to_p]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k+val_to_p],k%num_qu)

  return circ 


#U4 hessian 
def U4_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,div,ent_gate):

  q_reg = QuantumRegister(num_qu, "q")
  circ = QuantumCircuit(q_reg, name="layer-")

  val_to_p = 0
  if div >= 1:
   
      sh2 = sh2 - num_qu*div
      shift_pos = shift_pos - num_qu*div
      val_to_p = num_qu * div


  for k in range(val_g.shape[0]): 
    #rotation part
    sig = 1
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
          circ.rx(par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.rx(par_var[k+val_to_p],k%num_qu)


    elif val_g[k]==2:
      if sh2 == shift_pos and k == sh2:
          circ.ry(par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.ry(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.ry(par_var[k+val_to_p],k%num_qu)

    elif val_g[k]==3:
      if sh2 == shift_pos and k == sh2:
          circ.rz(par_var[k+val_to_p]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(par_var[k+val_to_p]-sig*shift,k%num_qu)
      else:  
        circ.rz(par_var[k+val_to_p],k%num_qu)

    #entanglement part
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]
      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
    return circ
