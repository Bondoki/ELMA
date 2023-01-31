# ELMA
ELMA - Extensions of the LeMonADE library
Repository with updates, analyzers, and projects for sharing BFM stuff related to various topics.


## Installation

* Clone and Install `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install cmake (minimum version 2.8)
* Any instance of OpenMPI
* Just do for standard compilation:
 
````sh
    # generates the projects
    mkdir build
    cd build
    cmake -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/ ..
    make
````

## Test WL-BLENDER single machine
````sh
mpiexec --mca orte_base_help_aggregate 0 -n 6 ./WangLandauSimulatorExtendedShellInteractionBLENDER  -i LinearChain_N32_PerXYZ256_EShell_rC2.451_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 10000000 -r 500000 -b 500000 --min-win -45.0 --max-win -20.2  --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.005 --read-in-BFM 0 --walker 2 --CZero 1000000.0 --OneOverN 0.1 > /dev/null 2>&1  
mpiexec --mca orte_base_help_aggregate 0 -n 6 ./WangLandauSimulatorNextNeighborBLENDER  -i LinearChain_N32_PerXYZ256_NNShell_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 4000000 -r 500000 -b 500000  --min-win -97.0 --max-win -49.0 --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.005 --read-in-BFM 0 --walker 2 --CZero 1000000.0 --OneOverN 0.1 > /dev/null 2>&1
````

## Test WL-BLENDER multiple machine
````sh
# in total 5*39 = 195 number of processes
#SBATCH --nodes=5  # number of nodes
#SBATCH --ntasks-per-node=39 # number jobs on node
mpirun ./WangLandauSimulatorNextNeighborBLENDER  -i LinearChain_N32_PerXYZ256_NNShell_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 5000000000 -r 500000 -b 500000 --min-win -97.0 --max-win -49.0 --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.005 --read-in-BFM 1 --walker 5 > /dev/null 2>&1  
````
or  

````sh
# in total 3*39 = 117 number of processes
#SBATCH --nodes=3  # number of nodes
#SBATCH --ntasks-per-node=39 # number jobs on node
mpirun ./WangLandauSimulatorExtendedShellInteractionBLENDER  -i LinearChain_N32_PerXYZ256_EShell_rC2.451_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 10000000000 -r 500000 -b 500000 --min-win -45.0 --max-win -20.2  --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.005 --read-in-BFM 1 --walker 3 > /dev/null 2>&1  
````

NOTE1: specify the number of walker within one energy window (--walker)  
NOTE2: the number of energy windows (number of processes/number of walker) has to be odd.  
NOTE3: by using Wang-Landau-BLENDER algorithm (10.1103/physreve.102.063304) you have to set the simulation time tsim instead of modification factor threshold  
NOTE4: by default --CZero 1000000.0 --OneOverN 0.1  
NOTE5: by using Wang-Landau-BLENDER algorithm avoid all input parameter "--mod-factor", "--threshold-mod-factor", or "--threshold-mod-factor-1t" as they are meaningless. Otherwise the programm will terminate.  

## ToDo
* add implementation for Rg^2 calculation for NextNeighbor
* add implementation for Rg^2 calculation for ExtendedShellInteraction

## License

See the LICENSE in the root directory.
