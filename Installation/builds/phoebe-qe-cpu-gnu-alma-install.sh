#!/bin/bash

## CPU-only Phoebe / Quantum ESPRESSO environment
## Compiled with: GNU GCC/gfortran, OpenMPI 4.1.6, OpenBLAS 0.3.33,
##               ScaLAPACK 2.2.3, FFTW 3.3.10, Parallel HDF5 2.1.0
## Target: Zen 4 (x86_64), MPI+OpenMP, no GPU
##



conda deactivate # deactivate all conda env


##################################################
### System software 
##################################################


sudo dnf install binutils-devel
sudo dnf install libibverbs-devel rdma-core-devel
sudo dnf install hwloc-devel libevent-devel


sudo chown -R $USER:$USER /opt
sudo chmod -R 755 /opt/

export TEMPORARY=~/Downloads/PHOEBE


