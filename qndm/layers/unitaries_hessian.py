#!/usr/bin/python3

from qiskit import QuantumRegister, QuantumCircuit

####################################
#                                   #
#           Hessian QNDM            #      
#                                   #
#####################################

#first unitary trasformation: U1:|00...0>->|\psi(\theta - shift*e_(sh1)+shift*e_(sh2))
def U1_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):
    
    #quantum register
    q_reg = QuantumRegister(num_qu, "q")
    #quantum circuit
    circ = QuantumCircuit(q_reg)

    #parametrized single gate layer
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

        #no parametrized entanglement layer
        if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
            #the condition is satisfied when the code complete the parametrized single gate layer

            ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

            if num_qu!=1:
                r_cx = int(num_qu/2)
                for ii in range(r_cx):
                    ent[ent_gate](2*ii)

                if num_qu>2:
                    for j in range(r_cx-(num_qu+1)%2):
                        ent[ent_gate](2*j+1)

    return circ


#second unitary trasformation: U1_dag|\psi(\theta - shift*e_(sh1)+shift*e_(sh2))->|00...0>
def U1_hess_dag(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):
  
  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)

  
  for k in reversed(range(val_g.shape[0])): 
    sig = 1

    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 

      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
        
    #parametrized single gate layer
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
        circ.rx(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
        circ.ry(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
        circ.rz(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k],k%num_qu)

  return circ 

#third unitary trasformation: U2:|00...0>->|\psi(\theta - shift*e_(sh1)-shift*e_(sh2))
def U2_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):
  
  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)


  for k in range(val_g.shape[0]): 
  #parametrized single gate layer
      sig = -1
      if val_g[k]==1:
            if sh2 == shift_pos and k == sh2:
                circ.rx(par_var[k]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.rx(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rx(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.rx(par_var[k],k%num_qu)


      elif val_g[k]==2:
            if sh2 == shift_pos and k == sh2:
                circ.ry(par_var[k]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.ry(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.ry(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.ry(par_var[k],k%num_qu)

      elif val_g[k]==3:
            if sh2 == shift_pos and k == sh2:
                circ.rz(par_var[k]-2*shift,k%num_qu)
            elif k == shift_pos:
                circ.rz(par_var[k]-sig*shift,k%num_qu)
            elif k == sh2:
                sig = -1
                circ.rz(par_var[k]-sig*shift,k%num_qu)
            else:  
                circ.rz(par_var[k],k%num_qu)

      #no parametrized entanglement layer
      if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0:
        #the condition is satisfied when the code complete the parametrized single gate layer

        ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

        if num_qu!=1:
            r_cx = int(num_qu/2)
            for ii in range(r_cx):
                ent[ent_gate](2*ii)
            if num_qu>2:
                for j in range(r_cx-(num_qu+1)%2):
                    ent[ent_gate](2*j+1)

         
  return circ


#fourth unitary trasformation: U2_dag|\psi(\theta - shift*e_(sh1)-shift*e_(sh2))->|00...0>
def U2_hess_dag(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)

  for k in reversed(range(val_g.shape[0])): 
    sig = -1

  
    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 

      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)
          
        
    #parametrized single gate layer 
    if val_g[k]==1:

         
      if sh2 == shift_pos and k == sh2:
          circ.rx(-par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
          circ.ry(-par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
          circ.rz(-par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = -1
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k],k%num_qu)

  return circ 


#fifth unitary trasformation: U3:|00...0>->|\psi(\theta + shift*e_(sh1)-shift*e_(sh2))
def U3_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)

 
  for k in range(val_g.shape[0]): 
  #parametrized single gate layer
    sig = -1
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
          circ.rx(par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rx(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.rx(par_var[k],k%num_qu)


    elif val_g[k]==2:
      if sh2 == shift_pos and k == sh2:
          circ.ry(par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.ry(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = +1
        circ.ry(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.ry(par_var[k],k%num_qu)

    elif val_g[k]==3:
      if sh2 == shift_pos and k == sh2:
          circ.rz(par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rz(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.rz(par_var[k],k%num_qu)

    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      #the condition is satisfied when the code complete the parametrized single gate layer

      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)

  return circ


#sixth unitary trasformation: U3_dag|\psi(\theta + shift*e_(sh1)-shift*e_(sh2))->|00...0>
def U3_hess_dag(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)

  
  for k in reversed(range(val_g.shape[0])): 
    sig = -1

    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 

      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)

        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
        for ii in range(r_cx):
          ent[ent_gate](2*ii)

        
    #parametrized single gate layer   
    if val_g[k]==1:

         
      if sh2 == shift_pos and k == sh2:
          circ.rx(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rx(-par_var[k],k%num_qu)

    elif val_g[k]==2:

      if sh2 == shift_pos and k == sh2:
          circ.ry(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.ry(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.ry(-par_var[k],k%num_qu)

    elif val_g[k]==3:

      if sh2 == shift_pos and k == sh2:
          circ.rz(-par_var[k],k%num_qu)
      elif k == shift_pos:
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(-par_var[k]+sig*shift,k%num_qu)
      else:  
        circ.rz(-par_var[k],k%num_qu)

  return circ 


#seventh unitary trasformation: U3:|00...0>->|\psi(\theta + shift*e_(sh1)+shift*e_(sh2))
def U4_hess(val_g,par_var,num_qu,num_l,shift_pos,sh2,shift,ent_gate):

  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)


  for k in range(val_g.shape[0]): 
    
    #parametrized single gate layer
    sig = 1
    if val_g[k]==1:
      if sh2 == shift_pos and k == sh2:
          circ.rx(par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rx(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rx(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.rx(par_var[k],k%num_qu)


    elif val_g[k]==2:
      if sh2 == shift_pos and k == sh2:
          circ.ry(par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.ry(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.ry(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.ry(par_var[k],k%num_qu)

    elif val_g[k]==3:
      if sh2 == shift_pos and k == sh2:
          circ.rz(par_var[k]+2*shift,k%num_qu)
      elif k == shift_pos:
        circ.rz(par_var[k]-sig*shift,k%num_qu)
      elif k == sh2:
        sig = 1
        circ.rz(par_var[k]-sig*shift,k%num_qu)
      else:  
        circ.rz(par_var[k],k%num_qu)

    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      #the condition is satisfied when the code complete the parametrized single gate layer

      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)

        for ii in range(r_cx):
          ent[ent_gate](2*ii)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
  return circ

#DM


#Unitary trasformation: U1:|00...0>->|\psi(\theta +- shift*e_(sh1)+-shift*e_(sh2))
def U1_hess_DM(val_g,par_var,num_qu,num_l,shift1,shift2,sh1,sh2,ent_gate):
  
  #quantum register
  q_reg = QuantumRegister(num_qu, "q")
  #quantum circuit
  circ = QuantumCircuit(q_reg)


    #parametrized single gate layer
  for k in range(val_g.shape[0]): 
    if val_g[k]==1:
      if k == sh1 and k == sh2:
        circ.rx(par_var[k]-(shift1+shift2),k%num_qu)
      elif k == sh1:
        circ.rx(par_var[k]-shift1,k%num_qu)

      elif k == sh2:
        circ.rx(par_var[k]-shift2,k%num_qu) 

      else:  
        circ.rx(par_var[k],k%num_qu)

    elif val_g[k]==2:

      if k == sh1 and k == sh2:
        circ.ry(par_var[k]-(shift1+shift2),k%num_qu)

      elif k == sh1:
        circ.ry(par_var[k]-shift1,k%num_qu)

      elif k == sh2:
        circ.ry(par_var[k]-shift2,k%num_qu)

      else: 
        circ.ry(par_var[k],k%num_qu)
    elif val_g[k]==3:

      if k == sh1 and k == sh2:
        circ.rz(par_var[k]-(shift1+shift2),k%num_qu)

      elif k == sh1:
        circ.rz(par_var[k]-shift1,k%num_qu)

      elif k == sh2:
        circ.rz(par_var[k]-shift2,k%num_qu)

      else: 
        circ.rz(par_var[k],k%num_qu)

    
    #no parametrized entanglement layer
    if  k !=0 and (k+1)%(len(val_g)/(num_l)) == 0: 
      ent = [lambda x:circ.cx(x, x+1),lambda x:circ.swap(x, x+1)]

      if num_qu!=1:
        r_cx = int(num_qu/2)
        for ii in range(r_cx):
            ent[ent_gate](ii*2)
        if num_qu>2:
          for j in range(r_cx-(num_qu+1)%2):
            ent[ent_gate](2*j+1)
  return circ