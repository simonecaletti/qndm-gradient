#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit

#-------------------------------------------------------------------------------

#entaglement type gate

def ent(circ, x, ent_gate):
  if ent_gate == 0:
    return circ.cx(x, x+1)
  if ent_gate == 1:
    return circ.swap(x, x+1)
  



#####################################
#                                   #
#            Gradient               #      
#                                   #
#####################################

#Unitary trasformation: U1:|00...0>->|\psi(\theta - shift*e_(shift_position))
def U1(val_g, par_var, num_qu, num_l,  shift, shift_pos, ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)

  for k in range(val_g.shape[0]): 
    
  #parametrized single gate layer
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
   
    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      #the condition is satisfied when the code complete the parametrized single gate layer


      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent(circ,2*ii,ent_gate)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent(circ,2*j+1,ent_gate)
  return circ

#second unitary trasformation: U1_dag|\psi(\theta - shift*e_(shift_position))->|00...0>
def U1_dag(val_g,par_var,num_qu,num_l,shift,shift_pos,ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)
  
  for k in reversed(range(val_g.shape[0])): 
    #no parametrized entanglement layer
    if  k !=0 and (k+1)%((len(val_g)/(num_l))) == 0: 

      if num_qu!=1:

        r_cx = int(num_qu/2)
        if num_qu>2:

          for j in range(r_cx-(num_qu+1)%2):
            ent(circ,2*j+1,ent_gate)

        for ii in range(r_cx):
          ent(circ,2*ii,ent_gate)
        
    #parametrized single gate layer
    if val_g[k]==1:
      if k == shift_pos:
        circ.rx(-par_var[k]+shift,k%num_qu)
      else:  
        circ.rx(-par_var[k],k%num_qu)

    elif val_g[k]==2:
      if k == shift_pos:
        circ.ry(-par_var[k]+shift,k%num_qu)
      else: 
        circ.ry(-par_var[k],k%num_qu)
    
    elif val_g[k]==3:
      if k == shift_pos:
        circ.rz(-par_var[k]+shift,k%num_qu)
      else: 
        circ.rz(-par_var[k],k%num_qu)

  return circ

#third unitary trasformation: U2:|00...0>->|\psi(\theta + shift*e_(shift_position))
def U2(val_g,par_var,num_qu,num_l,shift,shift_pos,ent_gate):


  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)


  for k in range(val_g.shape[0]): 

  #parametrized single gate layer
    if val_g[k]==1:
      if k == shift_pos:
        circ.rx(par_var[k]+shift,k%num_qu)
      else:  
        circ.rx(par_var[k],k%num_qu)

    elif val_g[k]==2:
      if k == shift_pos:
        circ.ry(par_var[k]+shift,k%num_qu)
      else: 
        circ.ry(par_var[k],k%num_qu)

    elif val_g[k]==3:
      if k == shift_pos:
        circ.rz(par_var[k]+shift,k%num_qu)
      else: 
        circ.rz(par_var[k],k%num_qu)

    #no parametrized entanglement layer   
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
    #the condition is satisfied when the code complete the parametrized single gate layer

      if num_qu!=1:
        r_cx = int(num_qu/2)

        for ii in range(r_cx):
          ent(circ,2*ii,ent_gate)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent(circ,2*j+1,ent_gate)

  return circ
  


