# cm-cRBM
Code package for the paper "Correlation-enhanced Neural Networks as Interpretable Variational Quantum States"


We here provide a minimal package containing the sampling of a cRBM of several quantities. 
The stored cRBM examples correspond to eigenstates of the toric code with different applied magnetic fields. 

To run a simulation, the user might choose between computation of cRBM energies 
(contained in the folder 'compute_energies), the calculation of Wilson loop expectation values 
(contained in the folder 'Wilsonloop_expvalue') or the computation of the second Renyi entropies 
(contained in the folder 'Renyi_entropy').

In order for the simulations to be able to run on a standard laptop in less than 5 minutes, 
we provided the optimized weights for a small lattice (N=18 spins) in the folder 'weights' for different 
field directions, including the weights of excited states.

For minimal usage, change directory into either 'compute_energies', 'Wilsonloop_expvalue' or 'Renyi_entropy'.


####################### COMPILING THE CODE ###########################################################

The provided codes ending with ".cc" use the standard c++ library and no additional packages. They can
be compiled using the command "make" in the folder 'example_simulation_MonteCarlo'.
The executable "run_MC" is produced.

A c++11 compliant compiler is assumed. Parallel programming using OpenMP is used.
These codes are used to calculate expectation values via Monte Carlo sampling.

######################################################################################################


########################## RUNNING A SIMULATION ######################################################

Minimal usage:
Once compiled, the simulation can be performed by executing './run'.

The simulation should run in less than 5 minutes on a standard laptop.

The simulation then performs a computation of the chosen quantity of the ground state
along the field direction h(1,0,1) (self-dual line). The outcome containing two columns
(first column: field strength, second column: computed quantity) is saved into the file
- 'energies.txt' if the folder energy_compare is chosen
- 'expvalues_AsBp.txt' if the folder Wilsonloops_expvalue is chosen
- 'entropies.txt' if the folder Renyi_entropy is chosen

######################################################################################################


---------------------------------- additional options ------------------------------------------------

########################### Changing the field direction #############################################

The field direction can be chosen by opening the file 'main.cc' and setting the parameter
'field_direction' to either 1, 2 or 3. In particular, the values correspond to the field directions
- 'field_direction=1': field in (1,0,1)-direction (self-dual line)
- 'field_direction=2': field in (1,0.2,0.5)-direction 
- 'field_direction=3': field (h,0.2,h) (shifted self-dual line by h_y=0.2)

To run the simulation with changed field-direction, the code should be first re-compiled using 'make'.
Then, executing './run_MC' runs the simulation.

######################################################################################################


######################### Excited states #############################################################

The cRBM weights of the first excited state obtained via minimization of the cost function with added
orthogonality constraints are given for the field in (1,0.2,0.5)-direction (field_direction=2).

To compute the given quantities for the first excited state, open the file 'main.cc' and set 'exc=true'
together with 'field_direction=2'.

To run the simulation with changed field-direction, the code should be first re-compiled using 'make'.
Then, executing './run_MC' runs the simulation.

######################################################################################################

######################### Comparing with exact solution ##############################################

The lowest four energies computed via exact diagonalization are given as comparison contained in the
folder 'energy_compare'.

In particular, the exact energies for the field directions 1,2,3 are given in the files 
'exact_energies_fielddir$i.txt', where i=1,2,3. The first column corresponds to the field strength,
the next 4 columns to the energies of the four lowest-lying eigenstates.

#####################################################################################################

################################## additional minor changes #########################################

- The number of conducted Monte Carlo sweeps can be changed by opening the file 'sampler.cc' 
(or 'sampler_renyi.cc') and changing the parameter 'nsweeps'.

- The size of the bipartition considered in the calculation of the Renyi entropy is specified by
 the member variable 'A_' in the class 'Sampler_Renyi' and is set to half of the lattice size.
It can be changed by opening the file 'sampler_renyi.cc' and setting the variable 'A_' to a different 
value fulfilling '0<A_<2*L_*L_-1'.

- The simulations can be parallelized by opening the Makefile and enabling OpenMp by uncommenting the flap '-fopenmp' and uncommenting the line starting with '#pragma omp parallel...' in the file 'main.cc' (using 50 or 60 threads depending on the choice of simulation).

#####################################################################################################
