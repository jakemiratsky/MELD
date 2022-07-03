import mdtraj as md
from matplotlib import pyplot as plt
import sys
import numpy as np
from itertools import chain
traj=md.load_pdb('native_free.pdb')
protein = []
A = range(0,10)
B = range(20,30)
m1 = chain(A,B)
for i in m1:
    protein.append(i)
pair=[]
for i in range(len(protein)):
    for j in range(i+1, len(protein)):
        pair.append([protein[i], protein[j]])

dist = md.compute_contacts(traj, pair, scheme='CA')
print(len(dist[1]))
dis=np.reshape(dist[0],(len(dist[1]),1))

dist_list=dis.tolist()

file = open("protein_contacts.dat","w")
accpt=[]
for i in zip(pair, dist_list):
    j  = [val for sublist in i for val in sublist]
    
    if j[2] <= 0.8:
        accpt.append([j[0],j[1],j[2]])
        #file.write('%d %d\n'%(j[0], j[1]))
        file.write( "{} CA {} CA {}\n".format(j[0]+1, j[1]+1,float(j[2])))
        file.write("\n")
#print(accpt)
