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

## Test WL-MPI single machine
````sh
mpiexec --mca orte_base_help_aggregate 0 -n 6 ./WangLandauSimulatorNextNeighborMPI  -i LinearChain_N32_PerXYZ256_NNShell_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 100000000000 -r 10000 -b 50000 -f 1.00001163461111 --threshold-mod-factor 1.000000010 --min-win -30.2 --max-win 0.4 --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.005 --read-in-BFM 0 --walker 2 > /dev/null 2>&1
````

## Test WL-MPI multiple machine
````sh
# in total 5*39 = 195 number of processes
#SBATCH --nodes=5  # number of nodes
#SBATCH --ntasks-per-node=39
mpirun ./WangLandauSimulatorNextNeighborMPI -i LinearChain_N384_PerXYZ256_NNShell_SelfAttraction_E-0.40.bfm --min -10000.2 --max 0.6 --bins 25002 -m 100000000000 -r 100000 -b 500000 -f 1.01 --threshold-mod-factor 1.000000001 --min-win -1514.0 --max-win -525.0 --HGLnDOS guess.dat --dump 0 --overlap 0.66 --length-increase 0.001 --read-in-BFM 0 --walker 5 > /dev/null 2>&1
````

NOTE1: specify the number of walker within one energy window (--walker)  
NOTE2: the number of energy windows (number of processes/number of walker) has to be odd.

## ToDo
* add implementation for Rg^2 calculation for NextNeighbor
* add WL-MPI implementation for ExtendedShellInteraction
* add implementation for Rg^2 calculation for ExtendedShellInteraction

## License

See the LICENSE in the root directory.
