#!/bin/bash

# nvhpc must be installed in /opt/nvidia/hpc_sdk
# Everything is compiled. Only the system compilers are needed to be installed
# Use Nvidia HPC SDK for CPU build
# openblas and lapack and scalapack is provided by SDK
# HDF5 is compiled 

# Skip to section C if you have all the librariers installed as in section A
##################################################
# Section A1
# Directory for user
##################################################

# Will stay in build directory
echo build_NVHPC_GPU


#################################################
# Section A2
# Set PATH and LIB
#################################################
export NVHPC=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3
export PATH=$NVHPC/compilers/bin:$NVHPC/comm_libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=$NVHPC/compilers/lib:$NVHPC/cuda/lib64:$NVHPC/math_libs/lib64:$LD_LIBRARY_PATH

export QE=/home/$USER/Downloads/qe-7.5
export QEX=$QE/build_NVHPC_GPU/bin
export EPW=$QE/EPW/bin
export PATH=$QEX:$PATH

export OMP_NUM_THREADS=2
NPROC=1
NPOOL=1
NIMG=1


##################################################
# Section A3
# Ryzen 7970X tuning (Zen 4)
##################################################

# 7970X = Zen 4 → use znver4, not znver2
# -Ofast usually performs better for QE than -O3
# keep fPIC for shared libs

# Correct flags for the NVHPC compiler targeting an RTX 4090
export CFLAGS="-fast -fPIC"
export FCFLAGS="-fast -fPIC -acc -mp"



##################################################
# Section B1
# Nvidia SDK
# Includes mpi compilers
##################################################


export CC=mpicc
export CXX=mpicxx
export FC=mpifort


##################################################
# Section B2
# OpenBLAS optimized for EPYC
##################################################

# Nvidia HPC SDK


##################################################
# Section B3
# FFTW3 MPI build
##################################################

# NVIDIA HPC SDK ships nvFFT and QE can use it


# ##################################################
# Section B4
# # ScaLAPACK (fixed version)
# ##################################################

# Nvidia HPC SDK


##################################################
# Section B5
# Build Parallel HDF5
##################################################
    cd /tmp
    git clone https://github.com/HDFGroup/hdf5.git
    cd hdf5

    git checkout hdf5_2.1.0
    mkdir build && cd build

    cmake .. \
    -DCMAKE_INSTALL_PREFIX=/opt/hdf5_2.1.0_nv \
    -DCMAKE_C_COMPILER=nvc \
    -DCMAKE_Fortran_COMPILER=nvfortran \
    -DHDF5_BUILD_FORTRAN=ON \
    -DHDF5_ENABLE_PARALLEL=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DHDF5_BUILD_CPP_LIB=ON \
    -DHDF5_BUILD_TOOLS=ON \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF

    make -j 32
    sudo make install

    # Check if you have the .mod files for the Fortran modules needed to build QE
    find /opt/hdf5_2.1.0_nv -name "*.mod"



##################################################
# Section C
# Quantum ESPRESSO build. 7.5
##################################################
    cd ~/Downloads

    git clone --depth 1 --branch qe-7.5 https://github.com/QEF/q-e.git
    cd q-e

    mkdir build_NVHPC_GPU
    cd build_NVHPC_GPU

    cmake .. \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DCMAKE_BUILD_TYPE=Release \
    -DQE_ENABLE_MPI=ON \
    -DQE_ENABLE_SCALAPACK=ON \
    -DQE_ENABLE_HDF5=ON  \
    -DQE_FFTW_VENDOR=Internal \
    -DHDF5_ROOT=/opt/hdf5_2.1.0_nv \
    -DQE_ENABLE_CUDA=ON \
    -DCMAKE_CUDA_ARCHITECTURES=89


##################################################
# IF : Quantum ESPRESSO build. from gitlab, latest version. Jun 8, 2026 
##################################################

    cmake .. \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpifort     \
    -DCMAKE_BUILD_TYPE=Release     \
    -DQE_ENABLE_MPI=ON     \
    -DQE_ENABLE_SCALAPACK=ON     \
    -DQE_ENABLE_HDF5=ON      \
    -DQE_FFTW_VENDOR=Internal     \
    -DHDF5_ROOT=/opt/hdf5_2.1.0_nv     \
    -DQE_GPU=cuda
    

    make all -j$(nproc)


# Check if pw.x can find everything
    ldd ~/Downloads/q-e/build_NVHPC_GPU/bin/pw.x


##################################################
# Section D
# Cleanup
##################################################

    rm -rf /tmp/openmpi*
    rm -rf /tmp/OpenBLAS*
    rm -rf /tmp/fftw*
    rm -rf /tmp/scalapack*
    rm -rf /tmp/hdf5*
    rm -rf /tmp/v*.tar.gz
    rm -rf /tmp/*.tar.gz


