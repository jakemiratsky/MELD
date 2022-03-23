# MELD
Model Employing Limited Data
**1. Downloading MELD

MELD is a plugin from OpenMM, a high performance toolkit for molecular simulations, and can be downloaded using the following steps: 

*module load anaconda/py3 

*conda create -n MELD 

*source activate MELD

*conda install -c conda-forge openmm openmpi meld

This should ensure that the MELD environment has been created on your system. To test the functionality use: 

*(meld) python -m simtk.testInstallation

Assuming you are directly interacting with a compute node via interactive, you should see 4 platforms.

**2. Initialize system

The next step once the MELD environment has been successfully implemented, is to initialize your system. For this, please refer to the setup_MELD.py script. 
Note that the restraints provided in the aforementioned script will not coincide with your system and the script should be changed accordingly. One requirement is that youur system has the corresponding topology and parameter files. For this we can use the following command from AMBER:

*tleap -f tleap.in

This will produce the necessary parameter files from use in the setup_MELD.py script

**3. Run REMD

Once the correct restraints have been set up for your system, all that is left is to begin the replica exchange (See the job.sh script).

**4. Analysis

Once the job has completed, there are some further command line operations within MELD that we can run. First off, say that you wanted to determine if your replica exchange was sufficient, i.e. there was good exchange: 

*(meld) analyze_remd visualize_trace 
*(meld) analyze_remd visualize_fup

Moreover, we can isolate one of the replicas from REMD by: 

*(meld) extract_trajectory extract_traj_dcd --replica 0 traj0.dcd 

