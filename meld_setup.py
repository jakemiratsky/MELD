#!/usr/bin/env python 

#======================================================

#       File Name : meld_setup.py

#       Author : Jacob Miratsky

#       Affiliation : ASU, School of Molecular Sciences

#       Contact : jacoba1@asu.edu

#       Creation Date : 30-09-2022

#       Usage : In main.py... 
#from meld_setup import *
#bias_dict = {
#    'Dist1': {'a1':'1', 'a2':'1', 'b1':'1', 'b2':'1', 
#              'distance':'4', 'confidence':'1.00'}, 
#    'Cart1': {'a1':'1', 'a2':'1', 'delta':'1'} 
#            }
#
#For 1 distance group, set 'b1':0
#Dist1 distance in A
#Cart1 delta in A
#
#meld_dict = {
#        'reps':'4', 'steps':'5000', 'pdb':'.pdb', 
#        'psf':'.psf', 'output':'meld.out'}
#
#model = setup_scripts()
#
#model.setup(meld_dict, bias_dict, dist_rest, cart_rest)
#model.Exec('clean_model.pdb', dist_rest)
#
#model.setup(meld_dict, bias_dict, dist_rest, cart_rest)
#model.Exec('clean_model.pdb', dist_rest)
#
#
#                    ADVANDED USE 
#
#from meld_setup import *
#k = ['.90', '.80', '.70', '.60', '.50']
#for i in k:
#   bias_dict = {
#    'Dist1': {'a1':'1', 'a2':'1', 'b1':'1', 'b2':'1', 
#              'distance':'4', 'confidence':i}, 
#    'Cart1': {'a1':'1', 'a2':'1', 'delta':'1'} 
#            }
#   meld_dict = {
#        'reps':'4', 'steps':'5000', 'pdb':'.pdb', 
#        'psf':'.psf', 'output':'meld.out'}
#    os.mkdir('sim_%sconfidence'%i)
#    shutil.copy(f'meld_dict["pdb"]', f'Simulation/meld_dict["pdb"]')
#    os.chdir('sim_%sconfidence'%i)
#    model = setup_scripts()
#    model.setup(meld_dict, bias_dict, dist_rest, cart_rest)
#    model.Exec('clean_model.pdb', 3)
#======================================================
import os 
import shutil
import time
def bash(bashCommand):
    import subprocess
    # print("\nBashCommand:\n%s\n" % (bashCommand))
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    output = output.decode("utf-8")
    if not error == None:
        error = error.decode("utf-8")
    for line in output:
        if "ERROR" in line:
            print(line)
    print("StandardOutput:\n%s\nError:%s\n" % (output, error))
    return output, error
def new_dir(dire):
    if os.path.exists(dire):
        print('%s exists... recreating now...'%dire)
        shutil.rmtree(dire)
        os.makedirs(dire)
    else:
        os.makedirs(dire)
def mv_files(file, dire, copy=False):
    files = [f for f in file]
    if os.path.exists(dire):
        print('%s exists... recreating now...' % dire)
        shutil.rmtree(dire)
        os.makedirs(dire)
    else:
        os.mkdir(dire)
    if copy == False:
        for f in files:
            if os.path.exists(f):
               shutil.move(f, dire)
            else:
                print('No files were moved to %s...' % dire)
    elif copy == True:
        for f in files:
            if os.path.exists(f):
                ini = '%s' % f
                fin = '%s/%s' % (dire, f)
                shutil.copyfile(ini, fin)
            else:
                print('No files were moved to %s...' % dire)
class setup_scripts():
    def __init__(self):
        None
    def setup(self, meld_dict, bias_dict, dist_rest=None, cart_rest=None, multi_job=False):
        new_dir('Simulation')
        shutil.copy(f'{meld_dict["pdb"]}', f'Simulation/{meld_dict["pdb"]}')
        os.chdir('Simulation')
        temp_fin = float(float(meld_dict['reps']) * 10 + 300)
        if os.path.exists('TEMPLATES'):
            os.rmdir('TEMPLATES')
            os.mkdir('TEMPLATES')
        else:
            os.mkdir('TEMPLATES')
        file=open('setup_MELD.py', 'w')
        file.write(
f'''#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
import meld.system.montecarlo as mc
from meld.system.restraints import LinearRamp,ConstantRamp
from collections import namedtuple
import glob as glob
import os

N_REPLICAS = {meld_dict['reps']}
N_STEPS = {meld_dict['steps']}
BLOCK_SIZE = 500

def gen_state_templates(index, templates):                                      
    n_templates = len(templates)
    print((index,n_templates,index%n_templates))
    a = system.ProteinMoleculeFromPdbFile(templates[index%n_templates])
    b = system.SystemBuilder(forcefield="ff14sbside")
    c = b.build_system_from_molecules([a])
    pos = c._coordinates
    c._box_vectors=np.array([0.,0.,0.])
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy,c._box_vectors)

def get_cartesian_restraints(s, scaler, residues, delta=None, k = 250.):
    cart = []
    backbone = ['CA']
    for i in residues:
        for b in backbone:
            atom_index = s.index_of_atom(i,b) - 1
            x,y,z = s.coordinates[atom_index]/10
            rest = s.restraints.create_restraint('cartesian',scaler,
            LinearRamp(0,15,0,1), res_index = i, atom_name = b, x=x, 
            y=y, z=z, delta=delta, force_const=k)
            cart.append(rest)
    return cart

def get_dist_restraints(filename, s, scaler):            
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1)) 
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4])                          

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),     
                                              r1=0, r2=0, r3=dist, r4=dist+0.5, k=350,
                                              atom_1_res_index=i, atom_2_res_index=j,
                                              atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def setup_system():
    templates = glob.glob('TEMPLATES/*.pdb')
    
    # build the system
    p = system.ProteinMoleculeFromPdbFile(templates[0])
    b = system.SystemBuilder(forcefield="ff14sbside")
    s = b.build_system_from_molecules([p])
    s.temperature_scaler = system.GeometricTemperatureScaler(0, 1.0, 300, {temp_fin})
    n_res = s.residue_numbers[-1]

    
    const_scaler = s.restraints.create_scaler('constant')   
    # nonconst_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)  

''') 
        file.close()
        for i in range(1, dist_rest+1):
            file=open('bias%s.py' % i, 'w')
            file.write(
f'''import mdtraj as md
import sys
import numpy as np
from itertools import chain

traj=md.load_pdb('TEMPLATES/clean_model.pdb')
if({bias_dict['Dist%s'%i]['b1']} == 0):
    complex = []
    a = range({bias_dict['Dist%s'%i]['a1']}-1,{bias_dict['Dist%s'%i]['a2']})
    for i in a:
        complex.append(i)
else:   
    complex = []
    a = range({bias_dict['Dist%s'%i]['a1']}-1,{bias_dict['Dist%s'%i]['a2']})
    b = range({bias_dict['Dist%s'%i]['b1']}-1,{bias_dict['Dist%s'%i]['b2']})
    m = chain(a,b)
    for i in m:
        complex.append(i)
pair=[]
for i in range(len(complex)):
    for j in range(i+1, len(complex)):
        pair.append([complex[i],complex[j]])

dist = md.compute_contacts(traj, pair, scheme='CA')
print(len(dist[1]))
dis=np.reshape(dist[0],(len(dist[1]),1))

dist_list=dis.tolist()

file = open('bias%s.dat','w')
accpt=[]
for i in zip(pair, dist_list):
    j  = [val for sublist in i for val in sublist]
    
    if 0 < j[2] <= {bias_dict['Dist%s'%i]['distance']}/10:
        accpt.append([j[0],j[1],j[2]])
        file.write( "{{}} CA {{}} CA {{}}\\n".format(j[0]+1, j[1]+1,float(j[2])))
        file.write("\\n")
''' % i)
            file.close()
        if dist_rest == 0:
            None 
        elif dist_rest !=0:
            for i in range(1, dist_rest+1):
                file=open('setup_MELD.py', 'a')
                file.write(
f'''
    distance_restraints_%s = get_dist_restraints('bias%s.dat',s,scaler=const_scaler) 
    s.restraints.add_selectively_active_collection(distance_restraints_%s, int(len(distance_restraints_%s)*{bias_dict['Dist%s'%i]['confidence']}))
''' % (i, i, i, i))
                file.close()
        if cart_rest == 0:
            None
        elif cart_rest != 0:
            for i in range(1, cart_rest+1):
                file=open('setup_MELD.py', 'a')
                file.write(
f'''
    cartesian_restraints_%s = list(range({bias_dict['Cart%s'%i]['a1']},{bias_dict['Cart%s'%i]['a2']}))
    cart_%s = get_cartesian_restraints(s, const_scaler, cartesian_restraints_%s, delta={bias_dict['Cart%s'%i]['delta']}/10)
    s.restraints.add_as_always_active_list(cart_%s)
''' % (i, i, i, i))
                file.close()
        file=open('setup_MELD.py', 'a')
        file.write(
f'''
   #create the options
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.use_big_timestep = False
    options.use_bigger_timestep = True
    options.cutoff = 1.8

    options.use_amap = False
    options.amap_alpha_bias = 1.0
    options.amap_beta_bias = 1.0
    options.timesteps = 11111
    options.minimize_steps = 500000
    options.min_mc = None
    options.run_mc = None
    for i in range(30):
        #print "!!!DO NOT USE MC MIN!!!"
        print ("Heads up! using MC minimizer!")
    #options.min_mc = sched

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=100)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)

    remd_runner = leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    #states = [gen_state(s, i) for i in range(N_REPLICAS)]
    states = [gen_state_templates(i,templates) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()''')
        file.close()
        file=open('tleap.in', 'w')
        file.write(
f'''set default PBradii mbondi3
source leaprc.protein.ff14SB
source leaprc.DNA.bsc1
source leaprc.RNA.OL3
mol = loadpdb clean_model.pdb 
saveamberparm mol system.top system.mdcrd
quit
''')
        file.close()
        file=open('AMBER.sh', 'w')
        file.write(
f'''#!/bin/bash
#SBATCH -p asinghargpu1
#SBATCH -q asinghargpu1
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH -n 2
#SBATCH -o AMBER.out

module purge
module load anaconda/py3
source activate ambertools
pdb4amber -i ../{meld_dict['pdb']} -y -o clean_model.pdb
echo Cleaning model for AMBER
tleap -f tleap.in
''')
        file.close()
        mv_files(['AMBER.sh', 'tleap.in'], 'AMBER')
        if multi_job == False:
            file=open('meld.sh', 'w')
            file.write(
f'''#!/bin/bash
#SBATCH -p asinghargpu1
#SBATCH -q asinghargpu1
#SBATCH -N {meld_dict['reps']}
#SBATCH -t 7-00:00:00
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -o {meld_dict['output']}

module purge
module load anaconda/py3
source activate meld
python setup_MELD.py
mpirun launch_remd
''' )
            file.close()
        elif multi_job == True:
            file=open('meld.sh', 'w')
            file.write(
f'''#!/bin/bash
#SBATCH -p gpu
#SBATCH -q wildfire 
#SBATCH -N {meld_dict['reps']}
#SBATCH -t 7-00:00:00
#SBATCH -c 8
#SBATCH --gres=gpu:4
#SBATCH -o {meld_dict['output']}

module purge
module load anaconda/py3
source activate meld
python setup_MELD.py
mpirun launch_remd
''' )
            file.close()
        file=open('Bias.sh', 'w')
        num_rest = dist_rest + cart_rest
        file.write(
f'''#!/bin/bash
#SBATCH -p asinghargpu1
#SBATCH -q asinghargpu1
#SBATCH -N 1
#SBATCH -t 0-01:00:00
#SBATCH -n 8
#SBATCH --gres=gpu:1
#SBATCH -o Bias_runs.out

module load anaconda/py3
source activate mdtraj

for i in {{1..%d}}; do
    for f in bias"$i".py; do python3 "$f"; done
    done
exit 
''' % num_rest)
        file.close()
    def Exec(self, file_path, num_bias):
        os.chdir('AMBER')
        bash('sbatch AMBER.sh')
        while not os.path.exists(file_path):
            time.sleep(5)
        if os.path.isfile(file_path):
            os.chdir('../')
            shutil.move('AMBER/clean_model.pdb', 'TEMPLATES/clean_model.pdb')
            bash('sbatch Bias.sh')
        else:
            raise ValueError('%s is not a file!'%file_path)

        while not os.path.exists('bias%d.dat'%num_bias):
            time.sleep(1)
        if os.path.isfile('bias%d.dat'%num_bias):
            bash('sbatch meld.sh')
        else:
            raise ValueError('bias%d.dat is not a file!'%num_bias)

