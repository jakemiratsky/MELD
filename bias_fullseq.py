import mdtraj as md
import sys
import numpy as np
import pandas as pd
from itertools import chain

# Load the distance matrix from a file using pandas
df1 = pd.read_csv('CRBN_fullseq_new1_cmap.txt', sep='\t')
df2 = pd.read_csv('CRBN_fullseq_new2_cmap.txt', sep='\t')
df3 = pd.read_csv('CRBN_fullseq_new3_cmap.txt', sep='\t')
df4 = pd.read_csv('CRBN_fullseq_new4_cmap.txt', sep='\t')
df5 = pd.read_csv('CRBN_fullseq_new5_cmap.txt', sep='\t')

# Convert the dataframe to a numpy array
matrix1 = df1.to_numpy()
matrix2 = df2.to_numpy()
matrix3 = df3.to_numpy()
matrix4 = df4.to_numpy()
matrix5 = df5.to_numpy()

# Creating the average matrix
ave_matrix = (matrix1+matrix2+matrix3+matrix4+matrix5)/5.0

## Create a new array with 1's where the original array has 1's, and 0's elsewhere
new_array = np.where(ave_matrix == 1, 1, 0)

# Getting the indices where value == 1
# Use numpy.where() to find the indices where the values are 1
indices = np.transpose(np.where(new_array == 1))

# Eliminate reflections of pairs across the diagonal
pairs = [pair for pair in indices if pair[0] <= pair[1]]

'''# if(536 == 0):
#     complex = []
#     a = range(29-1,227)
#     for i in a:
#         complex.append(i)
# else:   
#     complex = []
#     a = range(29-1,227)
#     b = range(536-1,734)
#     m = chain(a,b)
#     for i in m:
#         complex.append(i)
# pair=[]
# for i in range(len(complex)):
#     for j in range(i+1, len(complex)):
#         pair.append([complex[i],complex[j]])'''


traj=md.load_pdb('TEMPLATES/clean_model.pdb')

dist = md.compute_contacts(traj, pairs, scheme="CA")
print(len(dist[1]))
dis=np.reshape(dist[0],(len(dist[1]),1))

dist_list=dis.tolist()

file = open('bias_fullseq.dat','w')
accpt=[]
for i in zip(pairs, dist_list):
    j  = [val for sublist in i for val in sublist]
    
    if 0 < j[2] <= 1: # 10 A
        accpt.append([j[0],j[1],j[2]])
        file.write( "{} CA {} CA {}\n".format(j[0]+1, j[1]+1,float(j[2])))
        file.write("\n")

        