#!/usr/bin/python

#Code for gradient evaluation
#using the QNDM method.
#Written by: G. Minuto and S. Caletti
#Contacts: giovanni.minuto@uniroma1.it
#Cite 2301.07128 [quant-ph] in case you use 
#these code or part of it.

#gates counter for the circuit
def q_counter(gates_t,circ):

  gat=dict(circ.count_ops()) 


  for l in gat.keys():
        if l == 'rz': 
           gates_t[0] = gat[l] + gates_t[0]
             
        elif l == 'cx': 
           gates_t[1] = gat[l] + gates_t[1] 

        elif l == 'rx': 
           gates_t[2] = gat[l] + gates_t[2]

        elif l == 'ry': 
           gates_t[3] = gat[l] + gates_t[3]

        elif l == 'h': 
           gates_t[4] = gat[l] + gates_t[4]

        elif l == 'u1': 
           gates_t[5] = gat[l] + gates_t[5]
        
        elif l == 's': 
           gates_t[6] = gat[l] + gates_t[6]

        elif l == 'sdg': 
           gates_t[7] = gat[l] + gates_t[7]

        elif l == 'u2': 
           gates_t[8] = gat[l] + gates_t[8]

        elif l == 'u3': 
           gates_t[9] = gat[l] + gates_t[9]

        elif l == 'measure': 
           continue
         
        elif l == 'u': 
           gates_t[10] = gat[l] + gates_t[10]

        else:
            print(l)
            raise(ValueError('there are other gates'))
