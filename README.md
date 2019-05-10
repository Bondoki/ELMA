# ELMA
ELMA - Extensions of the LeMonADE library
Repository with updates, analyzers, and projects for sharing BFM stuff related to various topics.


## Installation

* Clone and Install `git clone https://github.com/LeMonADE-project/LeMonADE.git`
* Install cmake (minimum version 2.8)
* Just do for standard compilation:
 
````sh
    # generates the projects
    mkdir build
    cd build
    cmake -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE-library/include/ -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE-library/lib/ ..
    make
````

## Test WL-OpenMP
````sh
./WangLandauSimulatorNextNeighborOpenMP -i LinearChain_N32_PerXYZ128_NNS_E-0.4.bfm --min -10000.2 --max 0.6 --bins 25002 -m 100000000000 -r 10 -b 10000 -f 1.1 --min-win -95.0 --max-win 0.4 --HGLnDOS LinearChain_N32_PerXYZ128_NNS_E-0.4_HGLnDOS_shifted_01.dat --dump 0 --overlap 0.5 --length-increase 0.125 --read-in-BFM 0 --flatness 0.2 > /dev/null 2>&1
````



## License

See the LICENSE in the root directory.
